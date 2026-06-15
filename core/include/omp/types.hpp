#pragma once
#include "omp/vec3.hpp"

namespace omp {

// Cartesian state, inertial frame (km, km/s). For SGP4 this is TEME;
// for the two-body propagator it is whatever frame the input was given in.
struct StateVector {
    Vec3 r;  // position, km
    Vec3 v;  // velocity, km/s
};

// Classical (Keplerian) orbital elements.
struct OrbitalElements {
    double a = 0.0;      // semi-major axis, km
    double e = 0.0;      // eccentricity
    double i = 0.0;      // inclination, rad
    double raan = 0.0;   // right ascension of ascending node, rad
    double argp = 0.0;   // argument of perigee, rad
    double nu = 0.0;     // true anomaly, rad
};

// Geodetic position on/above the WGS84 ellipsoid.
struct Geodetic {
    double lat = 0.0;   // geodetic latitude, deg
    double lon = 0.0;   // longitude, deg [-180, 180]
    double alt = 0.0;   // altitude above ellipsoid, km
};

}  // namespace omp
