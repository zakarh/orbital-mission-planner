"""Pairwise conjunction (close-approach) screening over a time window."""
from __future__ import annotations

import numpy as np

from propagation import Satellite


def screen(sats: list[Satellite], jd0: float, window_min: float = 180.0,
           step_s: float = 30.0, threshold_km: float = 10.0,
           max_sats: int = 400):
    """Find satellite pairs approaching within `threshold_km`.

    Coarse grid screen: propagate every satellite onto a common time grid,
    then test pairwise distances at each step and keep the minimum approach
    per flagged pair. Returns events sorted by closest distance.
    """
    sats = sats[:max_sats]
    n = len(sats)
    if n < 2:
        return []

    steps = int(window_min * 60 / step_s) + 1
    times_min = np.arange(steps) * (step_s / 60.0)

    # positions[t] -> (n, 3) ECI km. Mark decayed/failed sats with NaN.
    pos = np.full((steps, n, 3), np.nan)
    for j, sat in enumerate(sats):
        m0 = sat.minutes_since_epoch(jd0)
        for ti, tm in enumerate(times_min):
            r, _ = sat.prop.propagate(m0 + tm)
            if r != [0.0, 0.0, 0.0]:
                pos[ti, j] = r

    iu, ju = np.triu_indices(n, k=1)
    best = np.full(len(iu), np.inf)
    best_t = np.zeros(len(iu))

    for ti in range(steps):
        p = pos[ti]
        d = p[iu] - p[ju]
        dist = np.sqrt(np.einsum("ij,ij->i", d, d))
        closer = dist < best
        best = np.where(closer, dist, best)
        best_t = np.where(closer, times_min[ti], best_t)

    events = []
    for k in range(len(iu)):
        if best[k] <= threshold_km:
            a, b = sats[iu[k]], sats[ju[k]]
            events.append({
                "sat_a": {"satnum": a.satnum, "name": a.name},
                "sat_b": {"satnum": b.satnum, "name": b.name},
                "min_distance_km": round(float(best[k]), 3),
                "time_offset_min": round(float(best_t[k]), 2),
            })
    events.sort(key=lambda e: e["min_distance_km"])
    return events
