#include "omp/transfer.hpp"
#include "omp/constants.hpp"
#include <cmath>

namespace omp {

HohmannTransfer hohmann(double r1, double r2, double mu) {
    HohmannTransfer t;
    t.r1 = r1;
    t.r2 = r2;
    t.transfer_a = 0.5 * (r1 + r2);

    t.v1 = std::sqrt(mu / r1);
    t.v2 = std::sqrt(mu / r2);
    t.v_transfer_peri = std::sqrt(mu * (2.0 / r1 - 1.0 / t.transfer_a));
    t.v_transfer_apo = std::sqrt(mu * (2.0 / r2 - 1.0 / t.transfer_a));

    t.dv1 = std::abs(t.v_transfer_peri - t.v1);
    t.dv2 = std::abs(t.v2 - t.v_transfer_apo);
    t.dv_total = t.dv1 + t.dv2;

    t.transfer_time = PI * std::sqrt(std::pow(t.transfer_a, 3) / mu);

    // Lead angle: how far ahead (rad) the target must be when we depart so it
    // arrives at the rendezvous point. Target sweeps n2 * transfer_time during
    // the coast; we must arrive 180 deg from departure.
    double n2 = std::sqrt(mu / std::pow(r2, 3));
    t.phase_angle = PI - n2 * t.transfer_time;
    // Normalize to (-pi, pi].
    while (t.phase_angle <= -PI) t.phase_angle += TWO_PI;
    while (t.phase_angle > PI) t.phase_angle -= TWO_PI;
    return t;
}

}  // namespace omp
