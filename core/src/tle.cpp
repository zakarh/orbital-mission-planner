#include "omp/tle.hpp"
#include "omp/transforms.hpp"
#include "omp/constants.hpp"
#include <cstdlib>
#include <stdexcept>
#include <cctype>

namespace omp {

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static double parse_d(const std::string& s) {
    std::string t = trim(s);
    return t.empty() ? 0.0 : std::atof(t.c_str());
}

// Parse the implied-decimal exponential fields (ndot'', bstar).
// signpos: column of the mantissa sign; mantissa is 5 digits; exp is 2 chars.
static double parse_exp(const std::string& line, int signpos, int mstart, int estart) {
    char sgn = line[signpos];
    std::string mant = line.substr(mstart, 5);
    std::string ex = trim(line.substr(estart, 2));
    if (ex.empty()) ex = "0";
    if (ex[0] != '+' && ex[0] != '-') ex = "+" + ex;  // implicit + exponent
    std::string num = std::string(1, (sgn == '-' ? '-' : '+')) + "0." + mant + "e" + ex;
    return std::atof(num.c_str());
}

TLE parse_tle(const std::string& line1, const std::string& line2,
              const std::string& name) {
    if (line1.size() < 64 || line2.size() < 63)
        throw std::invalid_argument("TLE lines too short");
    if (line1[0] != '1' || line2[0] != '2')
        throw std::invalid_argument("TLE line numbers invalid");

    TLE t;
    t.name = trim(name);
    t.satnum = std::atoi(line1.substr(2, 5).c_str());

    int yy = std::atoi(line1.substr(18, 2).c_str());
    t.epoch_year = (yy < 57) ? 2000 + yy : 1900 + yy;
    t.epoch_doy = parse_d(line1.substr(20, 12));
    // Jan 1 00:00 of the epoch year corresponds to day-of-year 1.0.
    t.epoch_jd = jday(t.epoch_year, 1, 1, 0, 0, 0.0) + (t.epoch_doy - 1.0);

    t.ndot = parse_d(line1.substr(33, 10));
    t.nddot = parse_exp(line1, 44, 45, 50);
    t.bstar = parse_exp(line1, 53, 54, 59);

    t.inclo = parse_d(line2.substr(8, 8)) * DEG2RAD;
    t.nodeo = parse_d(line2.substr(17, 8)) * DEG2RAD;
    t.ecco = std::atof(("0." + trim(line2.substr(26, 7))).c_str());
    t.argpo = parse_d(line2.substr(34, 8)) * DEG2RAD;
    t.mo = parse_d(line2.substr(43, 8)) * DEG2RAD;

    double no_revday = parse_d(line2.substr(52, 11));
    t.no_kozai = no_revday * TWO_PI / MINUTES_PER_DAY;  // rad/min
    return t;
}

}  // namespace omp
