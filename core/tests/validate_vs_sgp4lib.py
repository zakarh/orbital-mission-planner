"""Cross-validate the C++ SGP4 against the reference `sgp4` Python library."""
import omp_core as omp
from sgp4.api import Satrec

L1 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    8"
L2 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  105"

sat = omp.Sgp4(omp.parse_tle(L1, L2))
ref = Satrec.twoline2rv(L1, L2)

maxr = maxv = 0.0
for t in [0, 90, 180, 360, 720, 1440]:
    r, v = sat.propagate(t)
    e, rr, vv = ref.sgp4_tsince(t)
    dr = max(abs(r[i] - rr[i]) for i in range(3))
    dv = max(abs(v[i] - vv[i]) for i in range(3))
    maxr, maxv = max(maxr, dr), max(maxv, dv)
    print(f"t={t:5.0f} min  dpos={dr*1000:8.2f} m   dvel={dv*1e6:8.3f} mm/s")

print(f"\nMAX pos err = {maxr*1000:.3f} m,  vel err = {maxv*1e6:.3f} mm/s")
assert maxr < 0.01, "position mismatch vs reference sgp4"
assert maxv < 1e-5, "velocity mismatch vs reference sgp4"
print("C++ SGP4 matches reference library.")
