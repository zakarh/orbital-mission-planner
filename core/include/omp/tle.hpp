#pragma once
#include <string>

namespace omp {

// Parsed Two-Line Element set, normalized into SGP4-ready units.
struct TLE {
    std::string name;
    int satnum = 0;

    double epoch_jd = 0.0;   // Julian date (UTC) of the element epoch
    int epoch_year = 0;      // full 4-digit year
    double epoch_doy = 0.0;  // day of year with fraction

    double bstar = 0.0;      // drag term, 1/earth-radii
    double ndot = 0.0;       // mean motion 1st derivative, rev/day^2 / 2
    double nddot = 0.0;      // mean motion 2nd derivative, rev/day^3 / 6

    // SGP4 input elements (angles in radians, mean motion in rad/min).
    double inclo = 0.0;
    double nodeo = 0.0;
    double ecco = 0.0;
    double argpo = 0.0;
    double mo = 0.0;
    double no_kozai = 0.0;
};

// Parse a TLE from its two data lines (and optional name line).
// Throws std::invalid_argument on malformed input.
TLE parse_tle(const std::string& line1, const std::string& line2,
              const std::string& name = "");

}  // namespace omp
