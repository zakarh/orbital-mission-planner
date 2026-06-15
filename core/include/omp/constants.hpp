#pragma once

namespace omp {

// --- General math ---
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;

// --- WGS84 (used by the two-body propagator, transforms, transfers) ---
constexpr double MU_EARTH = 398600.4418;       // km^3 / s^2
constexpr double R_EARTH = 6378.137;           // km, equatorial radius
constexpr double R_EARTH_POLAR = 6356.752314245;  // km
constexpr double FLATTENING = 1.0 / 298.257223563;
constexpr double J2 = 1.08262668e-3;           // second zonal harmonic
constexpr double OMEGA_EARTH = 7.2921150e-5;   // rad/s, sidereal rotation rate

// --- WGS72 (SGP4 is defined against these — do not change) ---
constexpr double SGP4_MU = 398600.8;           // km^3 / s^2
constexpr double SGP4_RE = 6378.135;           // km
constexpr double SGP4_J2 = 0.001082616;
constexpr double SGP4_J3 = -0.00000253881;
constexpr double SGP4_J4 = -0.00000165597;

constexpr double MINUTES_PER_DAY = 1440.0;
constexpr double SECONDS_PER_DAY = 86400.0;

}  // namespace omp
