#include "omp/transforms.hpp"
#include "omp/constants.hpp"
#include <cmath>

namespace omp {

double jday(int year, int month, int day, int hour, int minute, double second) {
    // Standard astronomical Julian date (Vallado eq. 3-26).
    return 367.0 * year
           - std::floor((7 * (year + std::floor((month + 9) / 12.0))) * 0.25)
           + std::floor(275.0 * month / 9.0)
           + day + 1721013.5
           + ((second / 60.0 + minute) / 60.0 + hour) / 24.0;
}

double gmst(double jd_ut1) {
    // IAU-82 GMST expansion. tut1 = Julian centuries from J2000.
    double tut1 = (jd_ut1 - 2451545.0) / 36525.0;
    double temp = -6.2e-6 * tut1 * tut1 * tut1
                  + 0.093104 * tut1 * tut1
                  + (876600.0 * 3600.0 + 8640184.812866) * tut1
                  + 67310.54841;  // seconds
    temp = std::fmod(temp * DEG2RAD / 240.0, TWO_PI);  // 360/86400 = 1/240 deg/s
    if (temp < 0.0) temp += TWO_PI;
    return temp;
}

Vec3 eci_to_ecef(const Vec3& r, double th) {
    double c = std::cos(th), s = std::sin(th);
    return {c * r.x + s * r.y, -s * r.x + c * r.y, r.z};
}

Geodetic ecef_to_geodetic(const Vec3& r) {
    // Closed-form-ish iterative solution on the WGS84 ellipsoid.
    const double a = R_EARTH;
    const double e2 = FLATTENING * (2.0 - FLATTENING);
    double lon = std::atan2(r.y, r.x);
    double p = std::sqrt(r.x * r.x + r.y * r.y);
    double lat = std::atan2(r.z, p * (1.0 - e2));
    double alt = 0.0;
    for (int it = 0; it < 8; ++it) {
        double sl = std::sin(lat);
        double N = a / std::sqrt(1.0 - e2 * sl * sl);
        alt = p / std::cos(lat) - N;
        lat = std::atan2(r.z, p * (1.0 - e2 * N / (N + alt)));
    }
    Geodetic g;
    g.lat = lat * RAD2DEG;
    g.lon = lon * RAD2DEG;
    if (g.lon > 180.0) g.lon -= 360.0;
    if (g.lon < -180.0) g.lon += 360.0;
    g.alt = alt;
    return g;
}

Geodetic eci_to_geodetic(const Vec3& r_eci, double th) {
    return ecef_to_geodetic(eci_to_ecef(r_eci, th));
}

StateVector elements_to_state(const OrbitalElements& el, double mu) {
    double p = el.a * (1.0 - el.e * el.e);     // semi-latus rectum
    double cn = std::cos(el.nu), sn = std::sin(el.nu);
    double r = p / (1.0 + el.e * cn);

    // Position & velocity in the perifocal (PQW) frame.
    Vec3 r_pqw{r * cn, r * sn, 0.0};
    double sq = std::sqrt(mu / p);
    Vec3 v_pqw{-sq * sn, sq * (el.e + cn), 0.0};

    double ci = std::cos(el.i), si = std::sin(el.i);
    double cr = std::cos(el.raan), sr = std::sin(el.raan);
    double cw = std::cos(el.argp), sw = std::sin(el.argp);

    // Rotation PQW -> ECI (3-1-3: Rz(raan) Rx(i) Rz(argp)).
    double m11 = cr * cw - sr * sw * ci;
    double m12 = -cr * sw - sr * cw * ci;
    double m21 = sr * cw + cr * sw * ci;
    double m22 = -sr * sw + cr * cw * ci;
    double m31 = sw * si;
    double m32 = cw * si;

    StateVector sv;
    sv.r = {m11 * r_pqw.x + m12 * r_pqw.y,
            m21 * r_pqw.x + m22 * r_pqw.y,
            m31 * r_pqw.x + m32 * r_pqw.y};
    sv.v = {m11 * v_pqw.x + m12 * v_pqw.y,
            m21 * v_pqw.x + m22 * v_pqw.y,
            m31 * v_pqw.x + m32 * v_pqw.y};
    return sv;
}

OrbitalElements state_to_elements(const StateVector& sv, double mu) {
    const Vec3& r = sv.r;
    const Vec3& v = sv.v;
    double rn = r.norm();
    double vn2 = v.norm2();

    Vec3 h = r.cross(v);
    double hn = h.norm();
    Vec3 n{-h.y, h.x, 0.0};  // node vector (k x h)
    double nn = n.norm();

    // Eccentricity vector.
    Vec3 ev = ((vn2 - mu / rn) * r - r.dot(v) * v) / mu;
    double e = ev.norm();

    double energy = vn2 / 2.0 - mu / rn;
    OrbitalElements el;
    el.a = (std::abs(energy) > 1e-12) ? -mu / (2.0 * energy) : 0.0;
    el.e = e;
    el.i = std::acos(std::max(-1.0, std::min(1.0, h.z / hn)));

    if (nn > 1e-9) {
        el.raan = std::acos(std::max(-1.0, std::min(1.0, n.x / nn)));
        if (n.y < 0.0) el.raan = TWO_PI - el.raan;
    }
    if (nn > 1e-9 && e > 1e-9) {
        el.argp = std::acos(std::max(-1.0, std::min(1.0, n.dot(ev) / (nn * e))));
        if (ev.z < 0.0) el.argp = TWO_PI - el.argp;
    }
    if (e > 1e-9) {
        el.nu = std::acos(std::max(-1.0, std::min(1.0, ev.dot(r) / (e * rn))));
        if (r.dot(v) < 0.0) el.nu = TWO_PI - el.nu;
    } else {
        // Circular: fall back to argument of latitude.
        el.nu = std::atan2(r.dot(h.cross(n)) , r.dot(n));
    }
    return el;
}

}  // namespace omp
