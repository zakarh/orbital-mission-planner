// Three.js scene: Earth, satellites, orbits and transfer geometry.
// World units are kilometres / 1000. The frame is ECI (Earth-centered
// inertial) with +Z = north celestial pole, matching the backend.
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

export const SCALE = 1 / 1000; // km -> world units
const R_EARTH_KM = 6378.137;
// Longitude alignment offset for the equirectangular texture (tuned visually).
const TEXTURE_LON_OFFSET = -Math.PI / 2;

export class SpaceScene {
  constructor(container) {
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(
      45, window.innerWidth / window.innerHeight, 0.1, 5000);
    this.camera.up.set(0, 0, 1); // Z is up (north)
    this.camera.position.set(18, -22, 12);

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setSize(window.innerWidth, window.innerHeight);
    this.renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
    container.appendChild(this.renderer.domElement);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.minDistance = R_EARTH_KM * SCALE * 1.1;
    this.controls.maxDistance = 400;

    this._addLights();
    this._addStars();
    this._addEarth();

    this.satPoints = null;
    this.satnums = [];
    this.orbitLine = null;
    this.transferGroup = new THREE.Group();
    this.scene.add(this.transferGroup);

    this.raycaster = new THREE.Raycaster();
    this.raycaster.params.Points.threshold = 0.2;

    window.addEventListener("resize", () => this._onResize());
  }

  _addLights() {
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.35));
    const sun = new THREE.DirectionalLight(0xffffff, 1.4);
    sun.position.set(50, 30, 20);
    this.scene.add(sun);
  }

  _addStars() {
    const N = 3000, pos = new Float32Array(N * 3);
    for (let i = 0; i < N; i++) {
      const r = 800 + Math.random() * 1000;
      const t = Math.random() * Math.PI * 2;
      const p = Math.acos(2 * Math.random() - 1);
      pos[i * 3] = r * Math.sin(p) * Math.cos(t);
      pos[i * 3 + 1] = r * Math.sin(p) * Math.sin(t);
      pos[i * 3 + 2] = r * Math.cos(p);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    this.scene.add(new THREE.Points(g,
      new THREE.PointsMaterial({ color: 0x99aacc, size: 1.2 })));
  }

  _addEarth() {
    const r = R_EARTH_KM * SCALE;
    const geo = new THREE.SphereGeometry(r, 64, 48);
    const mat = new THREE.MeshPhongMaterial({ color: 0x224466, shininess: 8 });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.rotation.x = Math.PI / 2; // spin axis Y -> Z (north up)

    // Try to load a blue-marble texture; fall back to flat colour.
    new THREE.TextureLoader().load(
      "https://unpkg.com/three-globe@2.31.0/example/img/earth-blue-marble.jpg",
      (tex) => { mat.map = tex; mat.color.set(0xffffff); mat.needsUpdate = true; },
      undefined,
      () => console.warn("Earth texture failed; using flat colour"));

    this.earthSpin = new THREE.Group(); // rotates by GMST
    this.earthSpin.add(mesh);
    this.scene.add(this.earthSpin);

    // Faint equatorial circle for orientation (in the z=0 plane).
    const ringPts = [];
    for (let k = 0; k <= 128; k++) {
      const a = (k / 128) * Math.PI * 2;
      ringPts.push(new THREE.Vector3(
        Math.cos(a) * r * 1.001, Math.sin(a) * r * 1.001, 0));
    }
    const grid = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(ringPts),
      new THREE.LineBasicMaterial({ color: 0x335577 }));
    this.scene.add(grid);
  }

  setGmst(gmst) {
    if (this.earthSpin) this.earthSpin.rotation.z = gmst + TEXTURE_LON_OFFSET;
  }

  updateSatellites(sats) {
    const N = sats.length;
    const pos = new Float32Array(N * 3);
    this.satnums = new Array(N);
    this.satMeta = sats;
    for (let i = 0; i < N; i++) {
      pos[i * 3] = sats[i].r[0] * SCALE;
      pos[i * 3 + 1] = sats[i].r[1] * SCALE;
      pos[i * 3 + 2] = sats[i].r[2] * SCALE;
      this.satnums[i] = sats[i].satnum;
    }
    if (!this.satPoints) {
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
      const m = new THREE.PointsMaterial({ color: 0x4cc4ff, size: 0.35,
        sizeAttenuation: true });
      this.satPoints = new THREE.Points(g, m);
      this.scene.add(this.satPoints);
    } else {
      const attr = this.satPoints.geometry.getAttribute("position");
      if (attr.count !== N) {
        this.satPoints.geometry.setAttribute(
          "position", new THREE.BufferAttribute(pos, 3));
      } else {
        attr.array.set(pos);
        attr.needsUpdate = true;
      }
    }
  }

  drawOrbit(points) {
    this.clearOrbit();
    const v = points.map(p => new THREE.Vector3(
      p[0] * SCALE, p[1] * SCALE, p[2] * SCALE));
    const g = new THREE.BufferGeometry().setFromPoints(v);
    this.orbitLine = new THREE.Line(g,
      new THREE.LineBasicMaterial({ color: 0xffd34c }));
    this.scene.add(this.orbitLine);
  }

  clearOrbit() {
    if (this.orbitLine) {
      this.scene.remove(this.orbitLine);
      this.orbitLine.geometry.dispose();
      this.orbitLine = null;
    }
  }

  drawTransfer(geometry) {
    this.clearTransfer();
    const mk = (pts, color, width) => {
      const v = pts.map(p => new THREE.Vector3(
        p[0] * SCALE, p[1] * SCALE, p[2] * SCALE));
      const g = new THREE.BufferGeometry().setFromPoints(v);
      return new THREE.Line(g, new THREE.LineBasicMaterial({ color }));
    };
    this.transferGroup.add(mk(geometry.orbit1, 0x4cc4ff));
    this.transferGroup.add(mk(geometry.orbit2, 0x66ff99));
    this.transferGroup.add(mk(geometry.transfer, 0xff6b9d));

    // Burn markers at the ends of the transfer arc.
    const t = geometry.transfer;
    for (const p of [t[0], t[t.length - 1]]) {
      const dot = new THREE.Mesh(
        new THREE.SphereGeometry(0.25, 12, 12),
        new THREE.MeshBasicMaterial({ color: 0xff6b9d }));
      dot.position.set(p[0] * SCALE, p[1] * SCALE, p[2] * SCALE);
      this.transferGroup.add(dot);
    }
  }

  clearTransfer() {
    while (this.transferGroup.children.length) {
      const c = this.transferGroup.children.pop();
      c.geometry?.dispose();
      this.transferGroup.remove(c);
    }
  }

  pick(ndc) {
    if (!this.satPoints) return null;
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObject(this.satPoints);
    if (!hits.length) return null;
    const i = hits[0].index;
    return { satnum: this.satnums[i], meta: this.satMeta[i] };
  }

  _onResize() {
    this.camera.aspect = window.innerWidth / window.innerHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(window.innerWidth, window.innerHeight);
  }

  render() {
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }
}
