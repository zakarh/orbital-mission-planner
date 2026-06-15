"""Orbital Mission Planner REST API."""
from __future__ import annotations

import math
import time
import logging
from datetime import datetime, timezone

from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS

import omp_core as omp
from tle_source import get_tles, GROUPS
from propagation import Catalog, jd_from_datetime
from conjunction import screen

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("omp.api")

FRONTEND_DIR = "../frontend"
app = Flask(__name__, static_folder=None)
CORS(app)

# --- catalog cache (per group) -------------------------------------------------
_catalogs: dict[str, tuple[float, Catalog]] = {}
_CATALOG_TTL = 300  # seconds


def get_catalog(group: str) -> Catalog:
    now = time.time()
    cached = _catalogs.get(group)
    if cached and now - cached[0] < _CATALOG_TTL:
        return cached[1]
    cat = Catalog(get_tles(group))
    _catalogs[group] = (now, cat)
    log.info("Catalog group=%s: %d sats (%d errors)", group, len(cat.sats), cat.errors)
    return cat


def now_jd() -> float:
    return jd_from_datetime(datetime.now(timezone.utc))


# --- routes --------------------------------------------------------------------
@app.get("/api/health")
def health():
    return jsonify(status="ok", core="omp_core", mu_earth=omp.MU_EARTH)


@app.get("/api/groups")
def groups():
    return jsonify(groups=[{"id": k, "label": v} for k, v in GROUPS.items()])


@app.get("/api/satellites")
def satellites():
    group = request.args.get("group", "stations")
    cat = get_catalog(group)
    jd = now_jd()
    return jsonify(
        group=group,
        epoch=datetime.now(timezone.utc).isoformat(),
        gmst=omp.gmst(jd),
        count=len(cat.sats),
        satellites=cat.positions(jd),
    )


@app.get("/api/satellite/<int:satnum>/orbit")
def orbit(satnum: int):
    group = request.args.get("group", "stations")
    cat = get_catalog(group)
    sat = cat.get(satnum)
    if not sat:
        return jsonify(error="satellite not found"), 404
    jd = now_jd()
    r, v = sat.state_eci(jd)
    el = omp.state_to_elements(list(r), list(v))
    return jsonify(
        satnum=satnum,
        name=sat.name,
        gmst=omp.gmst(jd),
        period_min=sat.period_minutes(),
        deep_space=sat.deep_space,
        state={"r": list(r), "v": list(v)},
        elements={
            "a_km": el.a, "e": el.e,
            "i_deg": math.degrees(el.i),
            "raan_deg": math.degrees(el.raan),
            "argp_deg": math.degrees(el.argp),
            "nu_deg": math.degrees(el.nu),
            "apogee_km": el.a * (1 + el.e) - omp.R_EARTH,
            "perigee_km": el.a * (1 - el.e) - omp.R_EARTH,
        },
        orbit=sat.orbit_eci(jd, samples=160),
    )


@app.get("/api/satellite/<int:satnum>/groundtrack")
def groundtrack(satnum: int):
    group = request.args.get("group", "stations")
    minutes = float(request.args.get("minutes", 95))
    cat = get_catalog(group)
    sat = cat.get(satnum)
    if not sat:
        return jsonify(error="satellite not found"), 404
    return jsonify(satnum=satnum, name=sat.name,
                   track=sat.groundtrack(now_jd(), minutes))


@app.post("/api/conjunctions")
def conjunctions():
    body = request.get_json(silent=True) or {}
    group = body.get("group", "stations")
    window = float(body.get("window_min", 180))
    step = float(body.get("step_s", 30))
    threshold = float(body.get("threshold_km", 10))
    cat = get_catalog(group)
    sats = list(cat.sats.values())
    events = screen(sats, now_jd(), window_min=window, step_s=step,
                    threshold_km=threshold)
    return jsonify(
        group=group, screened=len(sats),
        window_min=window, threshold_km=threshold,
        count=len(events), conjunctions=events,
    )


@app.post("/api/transfer/hohmann")
def hohmann_transfer():
    body = request.get_json(silent=True) or {}
    alt1 = float(body.get("alt1_km", 400))
    alt2 = float(body.get("alt2_km", 35786))
    r1 = omp.R_EARTH + alt1
    r2 = omp.R_EARTH + alt2
    h = omp.hohmann(r1, r2)

    # Transfer-ellipse geometry in the equatorial plane (perigee at +X = r1).
    e_t = abs(r2 - r1) / (r2 + r1)
    a_t = h.transfer_a
    sign = 1.0 if r2 >= r1 else -1.0
    half = []
    for k in range(81):
        nu = math.pi * k / 80.0  # perigee -> apogee
        p = a_t * (1 - e_t * e_t)
        r = p / (1 + e_t * math.cos(nu))
        half.append([r * math.cos(nu), sign * r * math.sin(nu), 0.0])

    def circle(radius):
        return [[radius * math.cos(2 * math.pi * k / 128),
                 radius * math.sin(2 * math.pi * k / 128), 0.0]
                for k in range(129)]

    return jsonify(
        inputs={"alt1_km": alt1, "alt2_km": alt2, "r1_km": r1, "r2_km": r2},
        result={
            "dv1_kms": h.dv1, "dv2_kms": h.dv2, "dv_total_kms": h.dv_total,
            "transfer_time_s": h.transfer_time,
            "transfer_time_min": h.transfer_time / 60.0,
            "v1_kms": h.v1, "v2_kms": h.v2,
            "phase_angle_deg": math.degrees(h.phase_angle),
        },
        geometry={
            "orbit1": circle(r1),
            "orbit2": circle(r2),
            "transfer": half,
        },
    )


# --- static frontend -----------------------------------------------------------
@app.get("/")
def index():
    return send_from_directory(FRONTEND_DIR, "index.html")


@app.get("/<path:path>")
def static_files(path: str):
    return send_from_directory(FRONTEND_DIR, path)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
