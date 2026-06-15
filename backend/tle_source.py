"""Fetch and cache TLE sets from CelesTrak, with an offline fallback bundle."""
from __future__ import annotations

import os
import time
import logging
from dataclasses import dataclass

import requests

log = logging.getLogger("omp.tle")

CACHE_DIR = os.path.join(os.path.dirname(__file__), "cache")
CACHE_TTL = 2 * 3600  # seconds; CelesTrak asks clients not to refetch faster
CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php"

# Curated CelesTrak groups exposed to the UI.
GROUPS = {
    "stations": "Space stations",
    "active": "All active satellites",
    "starlink": "Starlink",
    "gps-ops": "GPS operational",
    "geo": "Geostationary",
    "visual": "Brightest (visual)",
}

# Minimal offline fallback so the app still runs without network access.
_FALLBACK = """ISS (ZARYA)
1 25544U 98067A   24168.51782528  .00016717  00000-0  30774-3 0  9994
2 25544  51.6394 211.6463 0004018  76.1736 283.9698 15.50108600456789
HUBBLE
1 20580U 90037B   24168.20366875  .00001234  00000-0  62000-4 0  9990
2 20580  28.4694 288.8421 0002628 321.7771  38.2772 15.09299865123456
GOES-16
1 41866U 16071A   24168.45000000 -.00000267  00000-0  00000+0 0  9991
2 41866   0.0412  88.1234 0000651 264.1234  95.1234  1.00271234 12345
"""


@dataclass
class Sat:
    name: str
    satnum: int
    line1: str
    line2: str


def _cache_path(group: str) -> str:
    return os.path.join(CACHE_DIR, f"{group}.tle")


def _parse_tle_text(text: str) -> list[Sat]:
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    sats: list[Sat] = []
    i = 0
    while i + 2 < len(lines) + 1:
        # Accept 3-line (name + L1 + L2) blocks.
        if i + 2 >= len(lines):
            break
        name, l1, l2 = lines[i], lines[i + 1], lines[i + 2]
        if l1.startswith("1 ") and l2.startswith("2 "):
            try:
                satnum = int(l1[2:7])
            except ValueError:
                satnum = 0
            sats.append(Sat(name.strip(), satnum, l1, l2))
            i += 3
        else:
            i += 1
    return sats


def get_tles(group: str = "stations", force: bool = False) -> list[Sat]:
    """Return TLEs for a group, using a fresh disk cache when available."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = _cache_path(group)

    if not force and os.path.exists(path):
        age = time.time() - os.path.getmtime(path)
        if age < CACHE_TTL:
            with open(path, "r", encoding="utf-8") as f:
                return _parse_tle_text(f.read())

    try:
        resp = requests.get(
            CELESTRAK_URL,
            params={"GROUP": group, "FORMAT": "tle"},
            timeout=15,
        )
        resp.raise_for_status()
        text = resp.text
        if "Invalid query" in text or len(text.strip()) < 10:
            raise ValueError("empty/invalid CelesTrak response")
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        log.info("Fetched %d bytes for group=%s", len(text), group)
        return _parse_tle_text(text)
    except Exception as exc:  # network down, rate-limited, etc.
        log.warning("CelesTrak fetch failed (%s); using cache/fallback", exc)
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                return _parse_tle_text(f.read())
        return _parse_tle_text(_FALLBACK)
