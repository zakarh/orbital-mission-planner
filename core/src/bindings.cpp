#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "omp/constants.hpp"
#include "omp/vec3.hpp"
#include "omp/types.hpp"
#include "omp/tle.hpp"
#include "omp/sgp4.hpp"
#include "omp/twobody.hpp"
#include "omp/transforms.hpp"
#include "omp/transfer.hpp"

namespace py = pybind11;
using namespace omp;

namespace {
using Tuple3 = std::array<double, 3>;
Tuple3 to_tuple(const Vec3& v) { return {v.x, v.y, v.z}; }
Vec3 from_tuple(const Tuple3& t) { return {t[0], t[1], t[2]}; }
}  // namespace

PYBIND11_MODULE(omp_core, m) {
    m.doc() = "Orbital Mission Planner C++ core (two-body RK4 + SGP4)";

    // --- constants ---
    m.attr("MU_EARTH") = MU_EARTH;
    m.attr("R_EARTH") = R_EARTH;
    m.attr("OMEGA_EARTH") = OMEGA_EARTH;

    // --- time / frames ---
    m.def("jday", &jday, py::arg("year"), py::arg("month"), py::arg("day"),
          py::arg("hour"), py::arg("minute"), py::arg("second"),
          "Julian date (UTC) from a calendar date.");
    m.def("gmst", &gmst, py::arg("jd_ut1"),
          "Greenwich Mean Sidereal Time (rad) for a UT1 Julian date.");
    m.def("eci_to_geodetic",
          [](Tuple3 r, double theta) {
              Geodetic g = eci_to_geodetic(from_tuple(r), theta);
              return std::array<double, 3>{g.lat, g.lon, g.alt};
          },
          py::arg("r_eci"), py::arg("gmst"),
          "ECI/TEME position (km) + GMST -> (lat_deg, lon_deg, alt_km).");

    // --- TLE ---
    py::class_<TLE>(m, "TLE")
        .def_readonly("name", &TLE::name)
        .def_readonly("satnum", &TLE::satnum)
        .def_readonly("epoch_jd", &TLE::epoch_jd)
        .def_readonly("epoch_year", &TLE::epoch_year)
        .def_readonly("epoch_doy", &TLE::epoch_doy)
        .def_readonly("bstar", &TLE::bstar)
        .def_readonly("inclo", &TLE::inclo)
        .def_readonly("nodeo", &TLE::nodeo)
        .def_readonly("ecco", &TLE::ecco)
        .def_readonly("argpo", &TLE::argpo)
        .def_readonly("mo", &TLE::mo)
        .def_readonly("no_kozai", &TLE::no_kozai);
    m.def("parse_tle", &parse_tle, py::arg("line1"), py::arg("line2"),
          py::arg("name") = "");

    // --- SGP4 ---
    py::class_<Sgp4Sat>(m, "Sgp4")
        .def(py::init([](const TLE& t) { return sgp4_init(t); }), py::arg("tle"))
        .def_readonly("epoch_jd", &Sgp4Sat::epoch_jd)
        .def_readonly("deep_space", &Sgp4Sat::deep_space)
        .def_readonly("valid", &Sgp4Sat::valid)
        .def("propagate",
             [](const Sgp4Sat& s, double minutes) {
                 StateVector sv = sgp4_propagate(s, minutes);
                 return std::make_pair(to_tuple(sv.r), to_tuple(sv.v));
             },
             py::arg("minutes_since_epoch"),
             "Propagate to minutes since epoch -> ((rx,ry,rz), (vx,vy,vz)) TEME.")
        .def("geodetic_at_jd",
             [](const Sgp4Sat& s, double jd) {
                 double minutes = (jd - s.epoch_jd) * 1440.0;
                 StateVector sv = sgp4_propagate(s, minutes);
                 Geodetic g = eci_to_geodetic(sv.r, gmst(jd));
                 return std::array<double, 3>{g.lat, g.lon, g.alt};
             },
             py::arg("jd"),
             "Sub-satellite point at a given UTC Julian date -> (lat,lon,alt).");

    // --- two-body RK4 ---
    m.def("propagate_twobody",
          [](Tuple3 r0, Tuple3 v0, double duration, double step, bool with_j2) {
              StateVector s0{from_tuple(r0), from_tuple(v0)};
              auto traj = propagate_twobody(s0, duration, step, with_j2);
              std::vector<std::pair<Tuple3, Tuple3>> out;
              out.reserve(traj.size());
              for (const auto& s : traj)
                  out.emplace_back(to_tuple(s.r), to_tuple(s.v));
              return out;
          },
          py::arg("r0"), py::arg("v0"), py::arg("duration"), py::arg("step"),
          py::arg("with_j2") = true,
          "RK4-integrate a state for `duration` s, sampling every `step` s.");

    // --- elements <-> state ---
    py::class_<OrbitalElements>(m, "OrbitalElements")
        .def_readonly("a", &OrbitalElements::a)
        .def_readonly("e", &OrbitalElements::e)
        .def_readonly("i", &OrbitalElements::i)
        .def_readonly("raan", &OrbitalElements::raan)
        .def_readonly("argp", &OrbitalElements::argp)
        .def_readonly("nu", &OrbitalElements::nu);
    m.def("state_to_elements",
          [](Tuple3 r, Tuple3 v, double mu) {
              return state_to_elements({from_tuple(r), from_tuple(v)}, mu);
          },
          py::arg("r"), py::arg("v"), py::arg("mu") = MU_EARTH);
    m.def("elements_to_state",
          [](double a, double e, double i, double raan, double argp, double nu,
             double mu) {
              OrbitalElements el{a, e, i, raan, argp, nu};
              StateVector sv = elements_to_state(el, mu);
              return std::make_pair(to_tuple(sv.r), to_tuple(sv.v));
          },
          py::arg("a"), py::arg("e"), py::arg("i"), py::arg("raan"),
          py::arg("argp"), py::arg("nu"), py::arg("mu") = MU_EARTH);

    // --- Hohmann transfer ---
    py::class_<HohmannTransfer>(m, "HohmannTransfer")
        .def_readonly("r1", &HohmannTransfer::r1)
        .def_readonly("r2", &HohmannTransfer::r2)
        .def_readonly("v1", &HohmannTransfer::v1)
        .def_readonly("v2", &HohmannTransfer::v2)
        .def_readonly("v_transfer_peri", &HohmannTransfer::v_transfer_peri)
        .def_readonly("v_transfer_apo", &HohmannTransfer::v_transfer_apo)
        .def_readonly("dv1", &HohmannTransfer::dv1)
        .def_readonly("dv2", &HohmannTransfer::dv2)
        .def_readonly("dv_total", &HohmannTransfer::dv_total)
        .def_readonly("transfer_time", &HohmannTransfer::transfer_time)
        .def_readonly("transfer_a", &HohmannTransfer::transfer_a)
        .def_readonly("phase_angle", &HohmannTransfer::phase_angle);
    m.def("hohmann", &hohmann, py::arg("r1"), py::arg("r2"),
          py::arg("mu") = MU_EARTH);
}
