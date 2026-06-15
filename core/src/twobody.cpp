#include "omp/twobody.hpp"
#include "omp/constants.hpp"
#include <cmath>

namespace omp {

Vec3 gravity_accel(const Vec3& r, bool with_j2) {
    double rn = r.norm();
    double rn3 = rn * rn * rn;
    Vec3 a = r * (-MU_EARTH / rn3);  // point-mass term

    if (with_j2) {
        double k = 1.5 * J2 * MU_EARTH * R_EARTH * R_EARTH / (rn3 * rn * rn);
        double z2r2 = (r.z * r.z) / (rn * rn);
        a.x += k * r.x * (5.0 * z2r2 - 1.0);
        a.y += k * r.y * (5.0 * z2r2 - 1.0);
        a.z += k * r.z * (5.0 * z2r2 - 3.0);
    }
    return a;
}

// Derivative of the state [r, v] -> [v, a].
static StateVector deriv(const StateVector& s, bool with_j2) {
    return {s.v, gravity_accel(s.r, with_j2)};
}

static StateVector add(const StateVector& s, const StateVector& k, double h) {
    return {s.r + k.r * h, s.v + k.v * h};
}

StateVector rk4_step(const StateVector& s, double dt, bool with_j2) {
    StateVector k1 = deriv(s, with_j2);
    StateVector k2 = deriv(add(s, k1, dt / 2.0), with_j2);
    StateVector k3 = deriv(add(s, k2, dt / 2.0), with_j2);
    StateVector k4 = deriv(add(s, k3, dt), with_j2);
    StateVector out;
    out.r = s.r + (k1.r + k2.r * 2.0 + k3.r * 2.0 + k4.r) * (dt / 6.0);
    out.v = s.v + (k1.v + k2.v * 2.0 + k3.v * 2.0 + k4.v) * (dt / 6.0);
    return out;
}

std::vector<StateVector> propagate_twobody(const StateVector& s0, double duration,
                                           double step, bool with_j2) {
    std::vector<StateVector> out;
    if (step <= 0.0) return out;

    // Cap the internal integration step for accuracy regardless of sampling.
    const double max_h = 10.0;  // seconds
    int substeps = std::max(1, (int)std::ceil(step / max_h));
    double h = step / substeps;

    StateVector s = s0;
    out.push_back(s);
    int nsamples = (int)std::floor(duration / step + 1e-9);
    for (int i = 0; i < nsamples; ++i) {
        for (int j = 0; j < substeps; ++j) s = rk4_step(s, h, with_j2);
        out.push_back(s);
    }
    return out;
}

}  // namespace omp
