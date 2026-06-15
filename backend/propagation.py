"""Propagation helpers built on the omp_core C++ extension."""
from __future__ import annotations

import math
from datetime import datetime, timezone

import omp_core as omp

from tle_source import Sat


def jd_from_datetime(dt: datetime) -> float:
    """UTC datetime -> Julian date (via the C++ jday)."""
    dt = dt.astimezone(timezone.utc)
    sec = dt.second + dt.microsecond / 1e6
    return omp.jday(dt.year, dt.month, dt.day, dt.hour, dt.minute, sec)


class Satellite:
    """A TLE-backed satellite with an initialized SGP4 propagator."""

    def __init__(self, sat: Sat):
        self.name = sat.name
        self.satnum = sat.satnum
        self.line1 = sat.line1
        self.line2 = sat.line2
        self.tle = omp.parse_tle(sat.line1, sat.line2, sat.name)
        self.prop = omp.Sgp4(self.tle)

    @property
    def epoch_jd(self) -> float:
        return self.prop.epoch_jd

    @property
    def deep_space(self) -> bool:
        return self.prop.deep_space

    def minutes_since_epoch(self, jd: float) -> float:
        return (jd - self.prop.epoch_jd) * 1440.0

    def state_eci(self, jd: float):
        """TEME/ECI position & velocity (km, km/s) at a Julian date."""
        return self.prop.propagate(self.minutes_since_epoch(jd))

    def geodetic(self, jd: float):
        """(lat_deg, lon_deg, alt_km) sub-satellite point at a Julian date."""
        return self.prop.geodetic_at_jd(jd)

    def period_minutes(self) -> float:
        # Mean motion (rad/min) -> orbital period.
        n = self.tle.no_kozai
        return (2 * math.pi / n) if n > 0 else 0.0

    def orbit_eci(self, jd: float, samples: int = 120):
        """Sample one orbital period of ECI positions starting at `jd`."""
        period = self.period_minutes()
        if period <= 0:
            return []
        m0 = self.minutes_since_epoch(jd)
        out = []
        for k in range(samples + 1):
            t = m0 + period * k / samples
            r, _ = self.prop.propagate(t)
            out.append(list(r))
        return out

    def groundtrack(self, jd: float, minutes: float, step_s: float = 30.0):
        """Sub-satellite lat/lon track over `minutes` from `jd`."""
        out = []
        n = int(minutes * 60 / step_s)
        for k in range(n + 1):
            j = jd + (k * step_s) / 86400.0
            lat, lon, alt = self.prop.geodetic_at_jd(j)
            out.append({"t": k * step_s, "lat": lat, "lon": lon, "alt": alt})
        return out


class Catalog:
    """Loads a group of TLEs and builds propagators, skipping bad entries."""

    def __init__(self, sats: list[Sat]):
        self.sats: dict[int, Satellite] = {}
        self.errors = 0
        for s in sats:
            try:
                self.sats[s.satnum] = Satellite(s)
            except Exception:
                self.errors += 1

    def get(self, satnum: int) -> Satellite | None:
        return self.sats.get(satnum)

    def positions(self, jd: float):
        """Current ECI position + sub-point for every satellite."""
        out = []
        for sat in self.sats.values():
            try:
                r, v = sat.state_eci(jd)
                if r == [0.0, 0.0, 0.0]:
                    continue  # decayed / propagation error
                lat, lon, alt = sat.geodetic(jd)
                out.append({
                    "satnum": sat.satnum,
                    "name": sat.name,
                    "r": list(r),
                    "lat": lat, "lon": lon, "alt": alt,
                    "deep_space": sat.deep_space,
                })
            except Exception:
                continue
        return out
