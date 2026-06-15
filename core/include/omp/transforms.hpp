#pragma once
#include "omp/types.hpp"

namespace omp {

// Julian date (UTC) from a calendar date. Valid for years 1900-2100.
double jday(int year, int month, int day, int hour, int minute, double second);

// Greenwich Mean Sidereal Time (radians, [0, 2pi)) for a given UT1 Julian date.
double gmst(double jd_ut1);

// Rotate an inertial (ECI/TEME) vector into the Earth-fixed (ECEF) frame
// using a simple rotation about the Z axis by the GMST angle.
Vec3 eci_to_ecef(const Vec3& r_eci, double theta_gmst);

// WGS84 ellipsoid: ECEF (km) -> geodetic latitude/longitude/altitude.
Geodetic ecef_to_geodetic(const Vec3& r_ecef);

// Convenience: inertial position + GMST -> geodetic sub-satellite point.
Geodetic eci_to_geodetic(const Vec3& r_eci, double theta_gmst);

// Classical elements <-> Cartesian state (gravitational parameter mu, km^3/s^2).
StateVector elements_to_state(const OrbitalElements& el, double mu);
OrbitalElements state_to_elements(const StateVector& sv, double mu);

}  // namespace omp
