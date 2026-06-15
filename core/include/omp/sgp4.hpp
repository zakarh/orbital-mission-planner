#pragma once
#include "omp/types.hpp"
#include "omp/tle.hpp"

namespace omp {

// SGP4 propagator state, initialized from a TLE. This is a near-Earth
// implementation (the SGP4 secular + periodic model). Deep-space (SDP4)
// lunar/solar resonance terms are NOT modeled; `deep_space` flags sats whose
// orbital period is >= 225 min, for which results are only approximate.
struct Sgp4Sat {
    // identity / epoch
    double epoch_jd = 0.0;
    bool deep_space = false;
    bool valid = false;

    // base mean elements (radians, rad/min)
    double ecco, inclo, nodeo, argpo, mo, bstar, no_unkozai;
    double gsto;  // Greenwich sidereal time at epoch (rad)

    // precomputed coefficients
    int isimp;
    double aycof, con41, cc1, cc4, cc5, d2, d3, d4, delmo, eta, argpdot,
        omgcof, sinmao, t2cof, t3cof, t4cof, t5cof, x1mth2, x7thm1, xlcof,
        xmcof, nodecf, mdot, nodedot;
};

// Initialize from a parsed TLE.
Sgp4Sat sgp4_init(const TLE& tle);

// Propagate to `minutes_since_epoch`. Returns TEME position (km) and
// velocity (km/s). On numerical failure, StateVector is left zeroed.
StateVector sgp4_propagate(const Sgp4Sat& s, double minutes_since_epoch);

}  // namespace omp
