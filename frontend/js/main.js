import * as THREE from "three";
import { api } from "./api.js";
import { SpaceScene } from "./scene.js";

const container = document.getElementById("canvas-container");
const space = new SpaceScene(container);

const els = {
  group: document.getElementById("group-select"),
  status: document.getElementById("status"),
  satInfo: document.getElementById("sat-info"),
  alt1: document.getElementById("alt1"),
  alt2: document.getElementById("alt2"),
  planBtn: document.getElementById("plan-transfer"),
  transferResult: document.getElementById("transfer-result"),
  conjThreshold: document.getElementById("conj-threshold"),
  conjWindow: document.getElementById("conj-window"),
  conjBtn: document.getElementById("run-conjunctions"),
  conjResult: document.getElementById("conj-result"),
  tooltip: document.getElementById("tooltip"),
};

let currentGroup = "stations";
let pollTimer = null;

function kv(label, value) {
  return `<div class="kv"><span>${label}</span><span>${value}</span></div>`;
}

async function init() {
  try {
    const { groups } = await api.groups();
    els.group.innerHTML = groups
      .map(g => `<option value="${g.id}">${g.label}</option>`).join("");
    els.group.value = currentGroup;
  } catch (e) {
    els.status.textContent = "backend offline?";
    console.error(e);
  }
  els.group.addEventListener("change", () => {
    currentGroup = els.group.value;
    space.clearOrbit();
    refresh();
  });
  startPolling();
  animate();
}

async function refresh() {
  try {
    const data = await api.satellites(currentGroup);
    space.setGmst(data.gmst);
    space.updateSatellites(data.satellites);
    els.status.textContent =
      `${data.count} satellites · ${new Date(data.epoch).toUTCString().slice(17, 25)} UTC`;
  } catch (e) {
    els.status.textContent = "fetch failed";
    console.error(e);
  }
}

function startPolling() {
  refresh();
  if (pollTimer) clearInterval(pollTimer);
  pollTimer = setInterval(refresh, 4000); // live update
}

// --- selection -----------------------------------------------------------
const ndc = new THREE.Vector2();

function toNDC(ev) {
  ndc.x = (ev.clientX / window.innerWidth) * 2 - 1;
  ndc.y = -(ev.clientY / window.innerHeight) * 2 + 1;
  return ndc;
}

space.renderer.domElement.addEventListener("click", async (ev) => {
  const hit = space.pick(toNDC(ev));
  if (!hit) return;
  els.satInfo.innerHTML = `Loading orbit for <b>${hit.meta.name}</b>…`;
  try {
    const o = await api.orbit(hit.satnum, currentGroup);
    space.drawOrbit(o.orbit);
    const e = o.elements;
    els.satInfo.innerHTML =
      `<b>${o.name}</b> (#${o.satnum})${o.deep_space ? " · deep-space" : ""}` +
      kv("Period", `${o.period_min.toFixed(1)} min`) +
      kv("Apogee", `${e.apogee_km.toFixed(0)} km`) +
      kv("Perigee", `${e.perigee_km.toFixed(0)} km`) +
      kv("Inclination", `${e.i_deg.toFixed(2)}°`) +
      kv("Eccentricity", e.e.toFixed(4)) +
      kv("RAAN", `${e.raan_deg.toFixed(1)}°`);
  } catch (e) {
    els.satInfo.textContent = "could not load orbit";
    console.error(e);
  }
});

space.renderer.domElement.addEventListener("mousemove", (ev) => {
  const hit = space.pick(toNDC(ev));
  if (hit) {
    els.tooltip.classList.remove("hidden");
    els.tooltip.style.left = ev.clientX + 12 + "px";
    els.tooltip.style.top = ev.clientY + 12 + "px";
    els.tooltip.textContent = `${hit.meta.name} · ${hit.meta.alt.toFixed(0)} km`;
  } else {
    els.tooltip.classList.add("hidden");
  }
});

// --- transfer ------------------------------------------------------------
els.planBtn.addEventListener("click", async () => {
  els.transferResult.textContent = "computing…";
  try {
    const data = await api.hohmann({
      alt1_km: parseFloat(els.alt1.value),
      alt2_km: parseFloat(els.alt2.value),
    });
    space.drawTransfer(data.geometry);
    const r = data.result;
    els.transferResult.innerHTML =
      kv("Burn 1 (Δv₁)", `${r.dv1_kms.toFixed(3)} km/s`) +
      kv("Burn 2 (Δv₂)", `${r.dv2_kms.toFixed(3)} km/s`) +
      kv("Total Δv", `<b>${r.dv_total_kms.toFixed(3)} km/s</b>`) +
      kv("Transfer time", `${r.transfer_time_min.toFixed(1)} min`) +
      kv("Phase angle", `${r.phase_angle_deg.toFixed(1)}°`);
  } catch (e) {
    els.transferResult.textContent = "transfer failed";
    console.error(e);
  }
});

// --- conjunctions --------------------------------------------------------
els.conjBtn.addEventListener("click", async () => {
  els.conjResult.textContent = "screening…";
  try {
    const data = await api.conjunctions({
      group: currentGroup,
      threshold_km: parseFloat(els.conjThreshold.value),
      window_min: parseFloat(els.conjWindow.value),
      step_s: 30,
    });
    if (!data.count) {
      els.conjResult.innerHTML =
        `<span class="muted">No approaches under ${data.threshold_km} km ` +
        `among ${data.screened} satellites.</span>`;
      return;
    }
    const rows = data.conjunctions.slice(0, 12).map(c =>
      `<div class="conj-row${c.min_distance_km < 1 ? " warn" : ""}">` +
      `${c.sat_a.name} ↔ ${c.sat_b.name}<br>` +
      `<b>${c.min_distance_km} km</b> at +${c.time_offset_min} min</div>`).join("");
    els.conjResult.innerHTML =
      `<div class="muted">${data.count} pair(s) screened from ` +
      `${data.screened} sats:</div>${rows}`;
  } catch (e) {
    els.conjResult.textContent = "screening failed";
    console.error(e);
  }
});

// --- render loop ---------------------------------------------------------
function animate() {
  requestAnimationFrame(animate);
  space.render();
}

init();
