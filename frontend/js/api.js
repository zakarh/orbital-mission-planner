// Thin wrapper around the Flask REST API.
const BASE = location.origin.includes("5000")
  ? "" // served by Flask itself
  : "http://127.0.0.1:5000";

async function get(path) {
  const r = await fetch(BASE + path);
  if (!r.ok) throw new Error(`${path} -> ${r.status}`);
  return r.json();
}

async function post(path, body) {
  const r = await fetch(BASE + path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!r.ok) throw new Error(`${path} -> ${r.status}`);
  return r.json();
}

export const api = {
  groups: () => get("/api/groups"),
  satellites: (group) => get(`/api/satellites?group=${group}`),
  orbit: (satnum, group) => get(`/api/satellite/${satnum}/orbit?group=${group}`),
  groundtrack: (satnum, group, minutes) =>
    get(`/api/satellite/${satnum}/groundtrack?group=${group}&minutes=${minutes}`),
  conjunctions: (body) => post("/api/conjunctions", body),
  hohmann: (body) => post("/api/transfer/hohmann", body),
};
