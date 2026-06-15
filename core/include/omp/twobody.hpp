#pragma once
#include "omp/types.hpp"
#include <vector>

namespace omp {

// Acceleration from Earth gravity. If with_j2 is true the dominant J2
// oblateness perturbation is included (important for realistic ground tracks).
Vec3 gravity_accel(const Vec3& r, bool with_j2);

// Propagate a Cartesian state forward by dt seconds using a single RK4 step.
StateVector rk4_step(const StateVector& s, double dt, bool with_j2);

// Propagate for `duration` seconds, sampling every `step` seconds.
// Returns one StateVector per sample including t=0. Integration substeps are
// capped so accuracy does not depend on the (possibly coarse) sample spacing.
std::vector<StateVector> propagate_twobody(const StateVector& s0,
                                           double duration,
                                           double step,
                                           bool with_j2);

}  // namespace omp
