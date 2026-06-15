"""Validate the C++ core against known references.

SGP4 is checked against the official Hoots/Vallado verification satellite
(catalog 88888). Two-body RK4 is checked for energy/period conservation.
"""
import math
import omp_core as omp

L1 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    8"
L2 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  105"


def approx(a, b, tol):
    return abs(a - b) <= tol


def test_sgp4_reference():
    tle = omp.parse_tle(L1, L2, "VANGUARD-TEST")
    sat = omp.Sgp4(tle)

    # tsince = 0 min
    r, v = sat.propagate(0.0)
    exp_r = (2328.97048951, -5995.22076416, 1719.97067261)
    exp_v = (2.91207230, -0.98341546, -7.09081703)
    for i in range(3):
        assert approx(r[i], exp_r[i], 0.1), f"r[{i}] {r[i]} != {exp_r[i]}"
        assert approx(v[i], exp_v[i], 1e-3), f"v[{i}] {v[i]} != {exp_v[i]}"

    # tsince = 360 min (Spacetrack Report #3 verification output for 88888)
    r, v = sat.propagate(360.0)
    exp_r = (2456.10705566, -6071.93853760, 1222.89727783)
    exp_v = (2.67939004, -0.44829081, -7.22879215)
    for i in range(3):
        assert approx(r[i], exp_r[i], 0.2), f"r360[{i}] {r[i]} != {exp_r[i]}"
        assert approx(v[i], exp_v[i], 1e-3), f"v360[{i}] {v[i]} != {exp_v[i]}"
    print("SGP4 reference OK")


def test_twobody_energy():
    # Circular LEO at 500 km altitude.
    r0 = (omp.R_EARTH + 500.0, 0.0, 0.0)
    vc = math.sqrt(omp.MU_EARTH / r0[0])
    v0 = (0.0, vc, 0.0)
    period = 2 * math.pi * math.sqrt(r0[0] ** 3 / omp.MU_EARTH)

    traj = omp.propagate_twobody(r0, v0, period, period / 100, False)
    # After one full period (no J2) we should return near the start.
    rf, vf = traj[-1]
    err = math.dist(rf, r0)
    assert err < 1.0, f"closure error {err} km"
    print(f"Two-body 1-orbit closure error: {err*1000:.1f} m")


def test_hohmann_leo_geo():
    # LEO (300 km) -> GEO. Textbook total dv ~ 3.9 km/s.
    r1 = omp.R_EARTH + 300.0
    r2 = 42164.0
    h = omp.hohmann(r1, r2)
    assert 3.8 < h.dv_total < 4.0, h.dv_total
    print(f"Hohmann LEO->GEO dv1={h.dv1:.4f} dv2={h.dv2:.4f} "
          f"total={h.dv_total:.4f} km/s, t={h.transfer_time/3600:.2f} h")


def test_elements_roundtrip():
    a, e, i = 7000.0, 0.01, math.radians(51.6)
    raan, argp, nu = math.radians(40), math.radians(60), math.radians(120)
    r, v = omp.elements_to_state(a, e, i, raan, argp, nu)
    el = omp.state_to_elements(r, v)
    assert approx(el.a, a, 1e-3) and approx(el.e, e, 1e-6)
    assert approx(el.i, i, 1e-9) and approx(el.raan, raan, 1e-9)
    print("Elements round-trip OK")


if __name__ == "__main__":
    test_sgp4_reference()
    test_twobody_energy()
    test_hohmann_leo_geo()
    test_elements_roundtrip()
    print("\nAll core tests passed.")
