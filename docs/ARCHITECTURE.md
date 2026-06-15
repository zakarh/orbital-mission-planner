# Architecture

This document describes how the Orbital Mission Planner is structured, how data
flows through it, and the conventions and design decisions that hold it
together. For setup and usage, see [README.md](README.md).

## 1. System overview

The application is a three-tier stack. Each tier has a single, well-defined
responsibility and talks to the next through a narrow interface.

```
                         pybind11                 REST / JSON over HTTP
   ┌──────────────────┐   (in-process)   ┌──────────────────┐   (localhost)   ┌────────────────┐
   │  C++ core         │ ───────────────▶ │  Flask backend    │ ──────────────▶ │  Three.js UI    │
   │  omp_core (.pyd)  │                  │  app.py           │                 │  browser        │
   │                   │ ◀─────────────── │                   │ ◀────────────── │                 │
   │ numeric engine    │   function calls │ orchestration     │   fetch()       │ presentation    │
   └──────────────────┘                  └──────────────────┘                 └────────────────┘
        no I/O, no state                   I/O, caching, glue                   rendering, input
```

**Guiding principle:** push all the numerically intensive, correctness-critical
math down into the C++ core; keep the backend a thin orchestration layer; keep
the frontend purely presentational. The backend never does orbital mechanics by
hand — it calls `omp_core`. The frontend never computes physics — it asks the
backend.

## 2. Layer responsibilities

| Layer | Owns | Does **not** own |
|-------|------|------------------|
| **C++ core** (`core/`) | Propagation, coordinate transforms, transfer math | Networking, caching, time-of-day, HTTP |
| **Backend** (`backend/`) | TLE fetch/cache, request handling, batching, screening orchestration | Orbital algorithms (delegates to core) |
| **Frontend** (`frontend/`) | 3D rendering, interaction, UI state | Physics, propagation, data persistence |

This separation means each layer is independently testable: the core is
validated numerically (`core/tests/`), the backend is checked by hitting its
endpoints, and the frontend is verified by rendering.

## 3. The C++ core (`omp_core`)

A single pybind11 extension module compiled from `core/src/*.cpp`. It is
**stateless and pure**: every function is a deterministic transformation of its
inputs, with no globals, no I/O, and no notion of "now". The caller (backend)
supplies absolute time as a Julian date.

### Module map

| File | Responsibility |
|------|----------------|
| `include/omp/vec3.hpp` | 3-vector value type (header-only) |
| `include/omp/constants.hpp` | WGS84 constants (general) **and** WGS72 constants (SGP4-only) |
| `include/omp/types.hpp` | `StateVector`, `OrbitalElements`, `Geodetic` |
| `src/twobody.cpp` | Gravity model + RK4 integrator |
| `src/sgp4.cpp` | SGP4 init + propagation (near-Earth) |
| `src/tle.cpp` | Fixed-column TLE parser |
| `src/transforms.cpp` | Julian date, GMST, ECI↔ECEF↔geodetic, elements↔state |
| `src/transfer.cpp` | Hohmann transfer solver |
| `src/bindings.cpp` | pybind11 surface — the only file that knows about Python |

### Key design decisions

- **Two gravity models, deliberately.** WGS84 constants drive the two-body
  propagator and geodetic transforms; WGS72 constants drive SGP4. SGP4 is
  *defined* against WGS72 — using WGS84 there would silently degrade accuracy.
  They are namespaced separately in `constants.hpp` to prevent cross-use. (This
  is also why the SGP4 source has a local `JJ2` constant — to avoid colliding
  with the global `omp::J2`.)

- **SGP4 is the near-Earth model only.** Deep-space (SDP4) lunar/solar resonance
  terms are not implemented. `Sgp4Sat::deep_space` flags objects with period
  ≥ 225 min so the upper layers can warn rather than silently mislead. This is
  the single largest scope cut and it is surfaced all the way to the UI.

- **The integrator decouples accuracy from sampling.** `propagate_twobody` caps
  the internal RK4 step at 10 s regardless of the requested sample spacing, so a
  caller asking for coarse output samples still gets an accurate trajectory.

- **The binding layer is the only Python-aware code.** `bindings.cpp` converts
  `Vec3` ↔ Python tuples and exposes classes; everything behind it is plain
  C++ that could be reused in a non-Python host.

### Validation

The core is the part most likely to be subtly wrong, so it is checked against
external references rather than itself:

- `tests/validate_vs_sgp4lib.py` cross-checks SGP4 against the independent
  `sgp4` Python package (matches to **< 1 mm** over 24 h).
- `tests/test_core.py` checks the Spacetrack #3 reference vector, two-body
  energy/period closure, a known LEO→GEO Δv, and element round-tripping.

## 4. The backend (`backend/`)

A Flask app that turns the stateless core into a stateful service: it knows what
"now" is, where to get TLEs, and how to batch many satellites.

### Module map

| File | Responsibility |
|------|----------------|
| `app.py` | Flask routes, catalog cache, serves the frontend |
| `tle_source.py` | CelesTrak fetch, disk cache (2 h TTL), offline fallback |
| `propagation.py` | `Satellite` / `Catalog` wrappers over `omp_core` + time helpers |
| `conjunction.py` | Vectorized pairwise close-approach screening |

### Request lifecycle

```
GET /api/satellites?group=stations
        │
        ▼
  get_catalog("stations")  ──cache miss──▶  get_tles()  ──▶  CelesTrak (or disk/fallback)
        │ (5-min in-memory cache)                                      │
        ▼                                                               ▼
  Catalog(sats)  ──builds──▶  Satellite → omp.parse_tle → omp.Sgp4   (one propagator per object)
        │
        ▼
  cat.positions(now_jd())  ──per sat──▶  Sgp4.propagate(minutes_since_epoch)
        │
        ▼
  JSON: { gmst, satellites:[{satnum,name,r,lat,lon,alt}, …] }
```

### Two layers of caching

1. **Disk cache** (`tle_source.py`) — raw TLE text per group, 2 h TTL. Respects
   CelesTrak's request to not refetch aggressively, and provides resilience: if
   the network is down it serves stale cache, then a bundled fallback set.
2. **In-memory catalog cache** (`app.py`) — parsed `Catalog` objects per group,
   5 min TTL. Avoids re-parsing thousands of TLEs and re-initializing SGP4 on
   every request.

### Conjunction screening

`conjunction.py` is the one place the backend does heavy numeric work itself
(in numpy, not the core), because it is inherently a batch operation over the
whole catalog:

1. Propagate every satellite onto a common time grid → `pos[time, sat, xyz]`.
2. At each time step, compute all pairwise distances via the upper-triangular
   index pair `(iu, ju)` and `einsum`.
3. Track the minimum approach per pair; report pairs under the threshold.

Decayed/failed propagations are marked `NaN` so they drop out of the distance
comparison. `max_sats` bounds the O(n²) cost.

## 5. The frontend (`frontend/`)

Plain ES modules loaded directly by the browser via an import map — **no build
step, no bundler, no npm**. Three.js comes from a CDN.

### Module map

| File | Responsibility |
|------|----------------|
| `index.html` | Layout + import map for `three` |
| `js/api.js` | Typed wrapper over the REST endpoints |
| `js/scene.js` | `SpaceScene` — all Three.js: Earth, satellites, orbits, transfer, picking |
| `js/main.js` | App glue — polling, selection, UI panels, event wiring |

`scene.js` owns everything WebGL and exposes intent-level methods
(`updateSatellites`, `drawOrbit`, `drawTransfer`, `pick`). `main.js` never
touches Three.js primitives directly — it speaks in domain terms. This keeps the
rendering swappable and the app logic readable.

### Live update loop

`main.js` polls `/api/satellites` every 4 s, updating the satellite point cloud
in place (reusing the buffer when the count is unchanged) and rotating the Earth
to the returned GMST. Selection and transfer overlays are drawn on demand and
persist across polls.

## 6. Coordinate frames and conventions

This is the contract that lets the three tiers agree on where things are.

| Concept | Convention |
|---------|-----------|
| **Inertial frame** | TEME/ECI throughout. SGP4 outputs TEME; the backend and frontend treat it as the inertial frame. |
| **Units (core/backend)** | Kilometres, kilometres/second, radians internally; degrees at the API boundary for human-facing angles. |
| **Time** | UTC Julian date is the absolute time currency. The core converts to "minutes since TLE epoch" for SGP4. |
| **Earth rotation** | Sub-satellite points come from rotating ECI by GMST into ECEF, then WGS84 geodetic. |
| **Frontend world** | World units = km / 1000. Frame is ECI with **+Z = north**; the camera up-vector is set to Z. The Earth mesh is rotated by GMST so geography sits under the right sub-points. |

Because satellites are rendered in true ECI, their positions are physically
exact. Only the Earth *texture's* absolute longitude offset is tuned by eye
(`TEXTURE_LON_OFFSET` in `scene.js`) — a cosmetic alignment, not a physics one.

## 7. End-to-end data flow examples

**Clicking a satellite:**
```
click → scene.pick() raycasts the point cloud → satnum
     → api.orbit(satnum) → GET /api/satellite/<id>/orbit
     → backend: Sgp4.propagate() sampled over one period + state_to_elements()
     → JSON {orbit:[…ECI…], elements:{…}}
     → scene.drawOrbit() + main.js renders the element panel
```

**Planning a transfer:**
```
"Plan & visualize" → api.hohmann({alt1,alt2}) → POST /api/transfer/hohmann
     → backend: omp.hohmann(r1,r2) for Δv/timing + builds ellipse & circle geometry
     → JSON {result:{dv1,dv2,dv_total,…}, geometry:{orbit1,orbit2,transfer}}
     → scene.drawTransfer() draws orbits + arc + burn markers; panel shows Δv
```

## 8. Build and packaging

- The core is built with **setuptools + pybind11** (`core/setup.py`), not CMake,
  and installed editable (`pip install -e core`). On Windows, setuptools locates
  MSVC automatically via vswhere.
- The backend serves the frontend directly (`send_from_directory`), so the whole
  app runs from one process on `127.0.0.1:5000` — no separate web server, no
  CORS hop in the default setup (CORS is enabled anyway for split deployments).

## 9. Extension points

The structure is designed so the obvious next features slot in cleanly:

| To add… | Touch… |
|---------|--------|
| SDP4 / deep-space accuracy | `src/sgp4.cpp` (add the resonance path); `deep_space` flag already plumbed |
| Drag/SRP in the numeric integrator | `gravity_accel()` in `twobody.cpp` |
| Plane-change or bi-elliptic transfers | new solver in `transfer.cpp` + a binding + an endpoint |
| Ground-station passes / visibility | new backend module over existing propagation |
| Client-side propagation (smoother motion) | port SGP4 to JS, or stream states; `scene.updateSatellites` already decoupled |
| New satellite catalogs | add a group id to `GROUPS` in `tle_source.py` |
```
