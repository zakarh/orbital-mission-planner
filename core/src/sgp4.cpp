// Near-Earth SGP4 propagator (WGS72), ported from the Hoots/Vallado reference
// formulation. Deep-space (SDP4) terms are intentionally omitted.
#include "omp/sgp4.hpp"
#include "omp/transforms.hpp"
#include "omp/constants.hpp"
#include <cmath>

namespace omp {

namespace {
const double XKE = 60.0 / std::sqrt(SGP4_RE * SGP4_RE * SGP4_RE / SGP4_MU);
const double JJ2 = SGP4_J2;
const double J3OJ2 = SGP4_J3 / SGP4_J2;
const double JJ4 = SGP4_J4;
const double X2O3 = 2.0 / 3.0;
const double VKMPERSEC = SGP4_RE * XKE / 60.0;
}  // namespace

Sgp4Sat sgp4_init(const TLE& tle) {
    Sgp4Sat s;
    s.epoch_jd = tle.epoch_jd;
    s.ecco = tle.ecco;
    s.inclo = tle.inclo;
    s.nodeo = tle.nodeo;
    s.argpo = tle.argpo;
    s.mo = tle.mo;
    s.bstar = tle.bstar;

    const double no_kozai = tle.no_kozai;
    const double cosio = std::cos(s.inclo);
    const double cosio2 = cosio * cosio;
    const double sinio = std::sin(s.inclo);
    const double eccsq = s.ecco * s.ecco;
    const double omeosq = 1.0 - eccsq;
    const double rteosq = std::sqrt(omeosq);

    // Recover original (un-Kozai'd) mean motion and semi-major axis.
    const double ak = std::pow(XKE / no_kozai, X2O3);
    const double d1 = 0.75 * JJ2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    double del = d1 / (ak * ak);
    const double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del = d1 / (adel * adel);
    s.no_unkozai = no_kozai / (1.0 + del);

    const double ao = std::pow(XKE / s.no_unkozai, X2O3);
    const double po = ao * omeosq;
    const double con42 = 1.0 - 5.0 * cosio2;
    s.con41 = -con42 - cosio2 - cosio2;  // 3*cosio2 - 1
    const double rp = ao * (1.0 - s.ecco);  // perigee, earth radii

    s.deep_space = (TWO_PI / s.no_unkozai) >= 225.0;
    s.isimp = (rp < (220.0 / SGP4_RE + 1.0)) ? 1 : 0;

    // Drag / atmospheric coefficients.
    const double ss = 78.0 / SGP4_RE + 1.0;
    double sfour = ss;
    double qzms24 = std::pow((120.0 - 78.0) / SGP4_RE, 4.0);
    const double perige = (rp - 1.0) * SGP4_RE;  // km
    if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) sfour = 20.0;
        qzms24 = std::pow((120.0 - sfour) / SGP4_RE, 4.0);
        sfour = sfour / SGP4_RE + 1.0;
    }

    const double pinvsq = 1.0 / (po * po);
    const double tsi = 1.0 / (ao - sfour);
    s.eta = ao * s.ecco * tsi;
    const double etasq = s.eta * s.eta;
    const double eeta = s.ecco * s.eta;
    const double psisq = std::abs(1.0 - etasq);
    const double coef = qzms24 * std::pow(tsi, 4.0);
    const double coef1 = coef / std::pow(psisq, 3.5);

    const double cc2 = coef1 * s.no_unkozai *
        (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
         0.375 * JJ2 * tsi / psisq * s.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    s.cc1 = s.bstar * cc2;

    double cc3 = 0.0;
    if (s.ecco > 1.0e-4)
        cc3 = -2.0 * coef * tsi * J3OJ2 * s.no_unkozai * sinio / s.ecco;
    s.x1mth2 = 1.0 - cosio2;
    s.cc4 = 2.0 * s.no_unkozai * coef1 * ao * omeosq *
        (s.eta * (2.0 + 0.5 * etasq) + s.ecco * (0.5 + 2.0 * etasq) -
         JJ2 * tsi / (ao * psisq) *
             (-3.0 * s.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
              0.75 * s.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
                  std::cos(2.0 * s.argpo)));
    s.cc5 = 2.0 * coef1 * ao * omeosq *
        (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    const double cosio4 = cosio2 * cosio2;
    const double temp1 = 1.5 * JJ2 * pinvsq * s.no_unkozai;
    const double temp2 = 0.5 * temp1 * JJ2 * pinvsq;
    const double temp3 = -0.46875 * JJ4 * pinvsq * pinvsq * s.no_unkozai;
    s.mdot = s.no_unkozai + 0.5 * temp1 * rteosq * s.con41 +
             0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    s.argpdot = -0.5 * temp1 * con42 +
                0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    const double xhdot1 = -temp1 * cosio;
    s.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                          2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;

    s.omgcof = s.bstar * cc3 * std::cos(s.argpo);
    s.xmcof = 0.0;
    if (s.ecco > 1.0e-4)
        s.xmcof = -X2O3 * coef * s.bstar / eeta;
    s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1;
    s.t2cof = 1.5 * s.cc1;

    if (std::abs(cosio + 1.0) > 1.5e-12)
        s.xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    else
        s.xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
    s.aycof = -0.5 * J3OJ2 * sinio;

    const double delmotemp = 1.0 + s.eta * std::cos(s.mo);
    s.delmo = delmotemp * delmotemp * delmotemp;
    s.sinmao = std::sin(s.mo);
    s.x7thm1 = 7.0 * cosio2 - 1.0;

    s.d2 = s.d3 = s.d4 = s.t3cof = s.t4cof = s.t5cof = 0.0;
    if (s.isimp != 1) {
        const double cc1sq = s.cc1 * s.cc1;
        s.d2 = 4.0 * ao * tsi * cc1sq;
        const double temp = s.d2 * tsi * s.cc1 / 3.0;
        s.d3 = (17.0 * ao + sfour) * temp;
        s.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s.cc1;
        s.t3cof = s.d2 + 2.0 * cc1sq;
        s.t4cof = 0.25 * (3.0 * s.d3 + s.cc1 * (12.0 * s.d2 + 10.0 * cc1sq));
        s.t5cof = 0.2 * (3.0 * s.d4 + 12.0 * s.cc1 * s.d3 + 6.0 * s.d2 * s.d2 +
                         15.0 * cc1sq * (2.0 * s.d2 + cc1sq));
    }

    s.gsto = gmst(s.epoch_jd);
    s.valid = true;
    return s;
}

StateVector sgp4_propagate(const Sgp4Sat& s, double t) {
    StateVector out;  // zeroed by default

    // Secular gravity and atmospheric drag.
    const double xmdf = s.mo + s.mdot * t;
    const double argpdf = s.argpo + s.argpdot * t;
    const double nodedf = s.nodeo + s.nodedot * t;
    double argpm = argpdf;
    double mm = xmdf;
    const double t2 = t * t;
    double nodem = nodedf + s.nodecf * t2;
    double tempa = 1.0 - s.cc1 * t;
    double tempe = s.bstar * s.cc4 * t;
    double templ = s.t2cof * t2;

    if (s.isimp != 1) {
        const double delomg = s.omgcof * t;
        const double delmt = 1.0 + s.eta * std::cos(xmdf);
        const double delm = s.xmcof * (delmt * delmt * delmt - s.delmo);
        const double temp = delomg + delm;
        mm = xmdf + temp;
        argpm = argpdf - temp;
        const double t3 = t2 * t;
        const double t4 = t3 * t;
        tempa = tempa - s.d2 * t2 - s.d3 * t3 - s.d4 * t4;
        tempe = tempe + s.bstar * s.cc5 * (std::sin(mm) - s.sinmao);
        templ = templ + s.t3cof * t3 + t4 * (s.t4cof + t * s.t5cof);
    }

    double nm = s.no_unkozai;
    double em = s.ecco;
    const double inclm = s.inclo;

    const double am = std::pow(XKE / nm, X2O3) * tempa * tempa;
    nm = XKE / std::pow(am, 1.5);
    em = em - tempe;
    if (em >= 1.0 || em < -0.001 || am < 0.95) return out;  // decayed / error
    if (em < 1.0e-6) em = 1.0e-6;

    mm = mm + s.no_unkozai * templ;
    double xlm = mm + argpm + nodem;
    nodem = std::fmod(nodem, TWO_PI);
    argpm = std::fmod(argpm, TWO_PI);
    xlm = std::fmod(xlm, TWO_PI);
    mm = std::fmod(xlm - argpm - nodem, TWO_PI);

    const double sinim = std::sin(inclm);
    const double cosim = std::cos(inclm);

    // Long-period periodics.
    const double axnl = em * std::cos(argpm);
    double temp = 1.0 / (am * (1.0 - em * em));
    const double aynl = em * std::sin(argpm) + temp * s.aycof;
    const double xl = mm + argpm + nodem + temp * s.xlcof * axnl;

    // Solve Kepler's equation for E + omega.
    const double u = std::fmod(xl - nodem, TWO_PI);
    double eo1 = u;
    double tem5 = 9999.9;
    double sineo1 = 0.0, coseo1 = 0.0;
    for (int ktr = 1; std::abs(tem5) >= 1.0e-12 && ktr <= 10; ++ktr) {
        sineo1 = std::sin(eo1);
        coseo1 = std::cos(eo1);
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (std::abs(tem5) >= 0.95) tem5 = (tem5 > 0.0) ? 0.95 : -0.95;
        eo1 += tem5;
    }

    // Short-period preliminary quantities.
    const double ecose = axnl * coseo1 + aynl * sineo1;
    const double esine = axnl * sineo1 - aynl * coseo1;
    const double el2 = axnl * axnl + aynl * aynl;
    const double pl = am * (1.0 - el2);
    if (pl < 0.0) return out;

    const double rl = am * (1.0 - ecose);
    const double rdotl = std::sqrt(am) * esine / rl;
    const double rvdotl = std::sqrt(pl) / rl;
    const double betal = std::sqrt(1.0 - el2);
    temp = esine / (1.0 + betal);
    const double sinu = am / rl * (sineo1 - aynl - axnl * temp);
    const double cosu = am / rl * (coseo1 - axnl + aynl * temp);
    double su = std::atan2(sinu, cosu);
    const double sin2u = (cosu + cosu) * sinu;
    const double cos2u = 1.0 - 2.0 * sinu * sinu;
    temp = 1.0 / pl;
    const double temp1 = 0.5 * JJ2 * temp;
    const double temp2 = temp1 * temp;

    // Short-period periodics applied to osculating elements.
    const double mrt = rl * (1.0 - 1.5 * temp2 * betal * s.con41) +
                       0.5 * temp1 * s.x1mth2 * cos2u;
    su = su - 0.25 * temp2 * s.x7thm1 * sin2u;
    const double xnode = nodem + 1.5 * temp2 * cosim * sin2u;
    const double xinc = inclm + 1.5 * temp2 * cosim * sinim * cos2u;
    const double mvt = rdotl - nm * temp1 * s.x1mth2 * sin2u / XKE;
    const double rvdot = rvdotl + nm * temp1 * (s.x1mth2 * cos2u + 1.5 * s.con41) / XKE;

    // Orientation vectors and final position/velocity (TEME).
    const double sinsu = std::sin(su), cossu = std::cos(su);
    const double snod = std::sin(xnode), cnod = std::cos(xnode);
    const double sini = std::sin(xinc), cosi = std::cos(xinc);
    const double xmx = -snod * cosi;
    const double xmy = cnod * cosi;
    const double ux = xmx * sinsu + cnod * cossu;
    const double uy = xmy * sinsu + snod * cossu;
    const double uz = sini * sinsu;
    const double vx = xmx * cossu - cnod * sinsu;
    const double vy = xmy * cossu - snod * sinsu;
    const double vz = sini * cossu;

    out.r = {mrt * ux * SGP4_RE, mrt * uy * SGP4_RE, mrt * uz * SGP4_RE};
    out.v = {(mvt * ux + rvdot * vx) * VKMPERSEC,
             (mvt * uy + rvdot * vy) * VKMPERSEC,
             (mvt * uz + rvdot * vz) * VKMPERSEC};
    return out;
}

}  // namespace omp
