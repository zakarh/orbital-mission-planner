#pragma once

namespace omp {

struct HohmannTransfer {
    double r1 = 0.0;          // initial circular radius, km
    double r2 = 0.0;          // final circular radius, km
    double v1 = 0.0;          // initial circular speed, km/s
    double v2 = 0.0;          // final circular speed, km/s
    double v_transfer_peri = 0.0;  // speed at transfer perigee, km/s
    double v_transfer_apo = 0.0;   // speed at transfer apogee, km/s
    double dv1 = 0.0;         // first burn magnitude, km/s
    double dv2 = 0.0;         // second burn magnitude, km/s
    double dv_total = 0.0;    // total delta-v, km/s
    double transfer_time = 0.0;    // half-ellipse coast time, s
    double transfer_a = 0.0;       // transfer ellipse semi-major axis, km
    double phase_angle = 0.0;      // required lead angle of target at departure, rad
};

// Classic two-burn Hohmann transfer between coplanar circular orbits.
// Radii are measured from Earth's center (km). mu in km^3/s^2.
HohmannTransfer hohmann(double r1, double r2, double mu);

}  // namespace omp
