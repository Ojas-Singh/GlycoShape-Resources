from __future__ import annotations

import json
import math
import re
import time
from copy import deepcopy
from html import escape
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
API_BASE = "https://glycoshape.org"

REPO_TARGET_PDB = ROOT / "inputs" / "REPO.pdb"
REPO_COMPLEX_PDB = ROOT / "inputs" / "REPO_complex.pdb"
DEFAULT_SCAN_GLYCAN_ID = "G00028MO"

REPO_DEFAULT_GLYCANS: List[Dict[str, Any]] = [
    {"bio_site": "N51", "chain": "A", "resi": 24, "glycan": "G52258NG", "source": "repo_default"},
    {"bio_site": "N65", "chain": "A", "resi": 38, "glycan": "G6878EM", "source": "repo_default"},
    {"bio_site": "N110", "chain": "A", "resi": 83, "glycan": "G63547QM", "source": "repo_default"},
    {"bio_site": "A126", "chain": "A", "resi": 126, "glycan": "G722008IY", "source": "repo_default"},
]


def clone_entries(entries: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [deepcopy(entry) for entry in entries]


def format_site_key(chain: str, resi: int) -> str:
    return f"{int(resi)}_{str(chain).strip() or 'A'}"


def split_site_key(site: str) -> Tuple[int, str]:
    resi, chain = site.split("_", 1)
    return int(resi), chain


def site_sort_key(site: str) -> Tuple[int, str]:
    try:
        return split_site_key(site)
    except Exception:
        return (10**9, site)


def entry_label(entry: Dict[str, Any]) -> str:
    chain = entry.get("chain", "A")
    resi = int(entry["resi"])
    bio = entry.get("bio_site")
    if bio and str(bio) != f"{chain}{resi}":
        return f"{chain}{resi} ({bio})"
    return f"{chain}{resi}"


def load_pdb_bytes(path_or_name: str) -> Tuple[str, bytes]:
    p = Path(path_or_name)
    candidates = [
        p,
        ROOT / path_or_name,
        ROOT / p.name,
        ROOT / "inputs" / p.name,
        Path("inputs") / p.name,
    ]
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate.name, candidate.read_bytes()
    raise FileNotFoundError(
        f"PDB not found: {path_or_name} (tried: " + ", ".join(str(x) for x in candidates) + ")"
    )


def colab_upload_one() -> Tuple[str, bytes]:
    from google.colab import files  # type: ignore

    uploaded = files.upload()
    name = next(iter(uploaded.keys()))
    return name, uploaded[name]


def get_pdb_input(default_filename: str) -> Tuple[str, bytes]:
    try:
        return load_pdb_bytes(default_filename)
    except Exception:
        pass

    try:
        import google.colab  # type: ignore

        _ = google.colab
        return colab_upload_one()
    except Exception:
        return load_pdb_bytes(default_filename)


def _raise_for_api(resp) -> None:
    if resp.ok:
        return
    try:
        payload = resp.json()
        msg = payload.get("error") or payload
    except Exception:
        msg = resp.text
    raise RuntimeError(f"API error {resp.status_code}: {msg}")


def create_session_from_pdb_bytes(filename: str, pdb_bytes: bytes, *, timeout_s: int = 600) -> Dict[str, Any]:
    import requests

    files = {"protFile": (filename, pdb_bytes, "chemical/x-pdb")}
    resp = requests.post(f"{API_BASE}/api/sessions", files=files, timeout=timeout_s)
    _raise_for_api(resp)
    return resp.json()


def submit_job(
    session_uuid: str,
    job_type: str,
    *,
    selected_glycans: Optional[Dict[str, str]] = None,
    parameters: Optional[Dict[str, Any]] = None,
    timeout_s: int = 60,
) -> Dict[str, Any]:
    import requests

    payload: Dict[str, Any] = {
        "session_uuid": session_uuid,
        "jobType": job_type,
    }
    if selected_glycans is not None:
        payload["selectedGlycans"] = selected_glycans
    if parameters is not None:
        payload["parameters"] = parameters

    resp = requests.post(f"{API_BASE}/api/jobs", json=payload, timeout=timeout_s)
    _raise_for_api(resp)
    return resp.json()


def get_job(job_uuid: str, *, timeout_s: int = 30) -> Dict[str, Any]:
    import requests

    resp = requests.get(f"{API_BASE}/api/jobs/{job_uuid}", timeout=timeout_s)
    _raise_for_api(resp)
    return resp.json()


def wait_job(job_uuid: str, *, poll_s: float = 2.0, timeout_s: int = 3600) -> Dict[str, Any]:
    t0 = time.time()
    last_status = None
    while True:
        meta = get_job(job_uuid)
        status = (
            meta.get("status")
            or (meta.get("queue") or {}).get("status")
            or ((meta.get("results") or {}).get("status"))
        )
        if status != last_status:
            print(f"[{job_uuid[:8]}] status: {status}")
            last_status = status

        if status in {"completed", "finished", "failed"}:
            if status == "failed":
                err = (meta.get("results") or {}).get("error") or meta.get("error")
                raise RuntimeError(f"Job failed: {job_uuid} ({err})")
            return meta

        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"Timed out waiting for job: {job_uuid}")

        time.sleep(poll_s)


def job_file_url(job_uuid: str, filename: str) -> str:
    return f"{API_BASE}/api/jobs/{job_uuid}/files/{filename}"


def scan_passed_sites(job_meta: Dict[str, Any]) -> List[str]:
    items = ((job_meta.get("results") or {}).get("results")) or []
    passed = [it.get("residue") for it in items if it.get("clash_solved") is True and it.get("residue")]
    seen = set()
    out: List[str] = []
    for residue in passed:
        if residue in seen:
            continue
        seen.add(residue)
        out.append(residue)
    return out


def glycan_entries_to_selected_map(entries: Iterable[Dict[str, Any]]) -> Dict[str, str]:
    selected: Dict[str, str] = {}
    for entry in entries:
        selected[format_site_key(entry["chain"], int(entry["resi"]))] = str(entry["glycan"])
    return selected


def glycan_entries_to_js(entries: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    js_items: List[Dict[str, Any]] = []
    seen = set()
    for entry in entries:
        chain = str(entry.get("chain", "A")).strip() or "A"
        resi = int(entry["resi"])
        key = (chain, resi)
        if key in seen:
            continue
        seen.add(key)
        item = {"chain": chain, "resi": resi}
        if entry.get("bio_site"):
            item["label"] = str(entry["bio_site"])
        js_items.append(item)
    return js_items


def entries_from_site_keys(
    sites: Iterable[str],
    *,
    glycan: str = DEFAULT_SCAN_GLYCAN_ID,
    source: str = "scan",
    site_name_prefix: str = "Site",
) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    for site in sites:
        try:
            resi, chain = split_site_key(site)
        except Exception:
            continue
        entries.append(
            {
                "chain": chain,
                "resi": int(resi),
                "glycan": glycan,
                "source": source,
                "bio_site": f"{site_name_prefix} {chain}{resi}",
            }
        )
    return entries


def parse_manual_patch(text: str, *, default_chain: str = "A") -> List[Dict[str, Any]]:
    tokens = [tok.strip() for tok in re.split(r"[\s,;]+", text or "") if tok.strip()]
    out: List[Dict[str, Any]] = []
    seen = set()
    for token in tokens:
        if "_" in token:
            try:
                resi, chain = split_site_key(token)
            except Exception as exc:
                raise ValueError(f"Could not parse manual patch token: {token}") from exc
        else:
            match = re.fullmatch(r"([A-Za-z]?)(\d+)", token)
            if not match:
                raise ValueError(f"Could not parse manual patch token: {token}")
            chain = match.group(1) or default_chain
            resi = int(match.group(2))
        key = (chain, resi)
        if key in seen:
            continue
        seen.add(key)
        out.append({"chain": chain, "resi": int(resi), "source": "manual"})
    return out


def parse_custom_glycan_json(text: str) -> List[Dict[str, Any]]:
    raw = (text or "").strip()
    if not raw:
        raise ValueError("CUSTOM_GLYCAN_JSON is empty.")
    data = json.loads(raw)
    if not isinstance(data, list) or not data:
        raise ValueError("CUSTOM_GLYCAN_JSON must be a non-empty JSON list.")
    out: List[Dict[str, Any]] = []
    seen = set()
    for idx, item in enumerate(data, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"Entry {idx} is not an object.")
        chain = str(item.get("chain") or "A").strip() or "A"
        if item.get("resi") is None:
            raise ValueError(f'Entry {idx} is missing "resi".')
        try:
            resi = int(item["resi"])
        except Exception as exc:
            raise ValueError(f'Entry {idx} has an invalid "resi": {item.get("resi")}') from exc
        glycan = str(item.get("glycan") or "").strip()
        if not glycan:
            raise ValueError(f'Entry {idx} is missing "glycan".')
        key = (chain, resi)
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "chain": chain,
                "resi": resi,
                "glycan": glycan,
                "bio_site": item.get("bio_site"),
                "source": "custom_json",
            }
        )
    return out


def _parse_pdb_atom_line(line: str) -> Optional[Dict[str, Any]]:
    if not line.startswith("ATOM"):
        return None
    try:
        return {
            "atom": line[12:16].strip(),
            "resn": line[17:20].strip(),
            "chain": (line[21] or "?").strip() or "?",
            "resi": int(line[22:26].strip()),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
        }
    except Exception:
        return None


def parse_hotspots_pdb_text(pdb_text: str) -> List[str]:
    hotspot_sites = set()
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        chain = (line[21] or "?").strip() or "?"
        resi = line[22:26].strip()
        try:
            b = float(line[60:66])
        except Exception:
            continue
        if b >= 99.9 and resi:
            hotspot_sites.add(f"{resi}_{chain}")
    return sorted(hotspot_sites, key=site_sort_key)


def _residue_centers_from_pdb_text(pdb_text: str) -> Dict[Tuple[str, int, str], Tuple[float, float, float]]:
    residues: Dict[Tuple[str, int, str], List[Tuple[str, float, float, float]]] = {}
    for line in pdb_text.splitlines():
        atom = _parse_pdb_atom_line(line)
        if not atom:
            continue
        key = (atom["chain"], atom["resi"], atom["resn"])
        residues.setdefault(key, []).append((atom["atom"], atom["x"], atom["y"], atom["z"]))

    centers: Dict[Tuple[str, int, str], Tuple[float, float, float]] = {}
    for key, atoms in residues.items():
        picked = None
        for atom_name, x, y, z in atoms:
            if atom_name == "CB":
                picked = (x, y, z)
                break
        if picked is None:
            for atom_name, x, y, z in atoms:
                if atom_name == "CA":
                    picked = (x, y, z)
                    break
        if picked is not None:
            centers[key] = picked
    return centers


def derive_preset_patch_from_complex(
    complex_path: str | Path,
    glycan_entries: Iterable[Dict[str, Any]],
    *,
    target_chain: str = "A",
    distance_cutoff: float = 8.0,
    cluster_gap: int = 2,
    max_residues: int = 8,
) -> List[Dict[str, Any]]:
    complex_text = Path(complex_path).read_text()
    centers = _residue_centers_from_pdb_text(complex_text)
    target = {k: v for k, v in centers.items() if k[0] == target_chain}
    binders = {k: v for k, v in centers.items() if k[0] != target_chain}
    if not target or not binders:
        return []

    glyco_residues = [int(entry["resi"]) for entry in glycan_entries if str(entry.get("chain", "A")) == target_chain]
    contacts = []
    for (chain, resi, resn), xyz in target.items():
        best = min(math.dist(xyz, binder_xyz) for binder_xyz in binders.values())
        if best <= float(distance_cutoff):
            contacts.append({"chain": chain, "resi": resi, "resn": resn, "best": best})
    contacts = sorted(contacts, key=lambda item: item["resi"])
    if not contacts:
        return []

    clusters: List[List[Dict[str, Any]]] = []
    current: List[Dict[str, Any]] = []
    for item in contacts:
        if not current or item["resi"] - current[-1]["resi"] <= int(cluster_gap):
            current.append(item)
        else:
            clusters.append(current)
            current = [item]
    if current:
        clusters.append(current)

    def cluster_score(cluster: List[Dict[str, Any]]) -> Tuple[Any, ...]:
        residues = [item["resi"] for item in cluster]
        if glyco_residues:
            min_glyco_distance = min(min(abs(resi - g) for g in glyco_residues) for resi in residues)
            near_glyco_count = sum(1 for resi in residues for g in glyco_residues if abs(resi - g) <= 6)
        else:
            min_glyco_distance = 10**6
            near_glyco_count = 0
        mean_contact = sum(item["best"] for item in cluster) / len(cluster)
        return (
            min_glyco_distance,
            -near_glyco_count,
            -len(cluster),
            mean_contact,
            residues[0],
        )

    best_cluster = sorted(clusters, key=cluster_score)[0]
    residues = [item["resi"] for item in best_cluster][: int(max_residues)]
    return [
        {
            "chain": target_chain,
            "resi": int(resi),
            "source": "repo_preset",
            "bio_site": f"Patch {target_chain}{int(resi)}",
        }
        for resi in residues
    ]


def render_entry_table(title: str, entries: Iterable[Dict[str, Any]], *, empty_text: str = "No entries selected yet.") -> str:
    rows = []
    for entry in entries:
        rows.append(
            "<tr>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{escape(entry_label(entry))}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{escape(str(entry.get('glycan', '')))}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{escape(str(entry.get('source', '')))}</td>"
            "</tr>"
        )
    body = "".join(rows) if rows else (
        "<tr><td colspan='3' style='padding:10px; color:#666; border-bottom:1px solid #eee;'>"
        + escape(empty_text)
        + "</td></tr>"
    )
    return (
        "<div style='font-family:system-ui, -apple-system, Segoe UI, Roboto, Arial; font-size:14px;'>"
        f"<div style='font-weight:700; margin-bottom:8px;'>{escape(title)}</div>"
        "<table style='border-collapse:collapse; width:100%;'>"
        "<thead><tr>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Residue</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Glycan</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Source</th>"
        "</tr></thead>"
        f"<tbody>{body}</tbody></table></div>"
    )


def build_3dmol_gallery_html(structures: List[Dict[str, Any]], *, height: int = 520) -> str:
    html = """
<style>
  .rg-viewer-shell {
    font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial;
    font-size: 14px;
  }
  .rg-viewer-nav {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 12px;
    margin: 8px 0;
  }
  .rg-viewer-title {
    font-weight: 700;
  }
  .rg-viewer-sub {
    color: #666;
    font-size: 12px;
    margin-top: 2px;
  }
  .rg-viewer-btn {
    appearance: none;
    border: 1px solid #ddd;
    background: #fafafa;
    border-radius: 10px;
    padding: 8px 12px;
    font-weight: 600;
    cursor: pointer;
  }
  .rg-viewer-btn:hover {
    background: #f2f2f2;
  }
  #rg-3dmol {
    width: 100%;
    height: __HEIGHT__px;
    position: relative;
    border: 1px solid #eee;
    border-radius: 12px;
    overflow: hidden;
  }
</style>

<script src="https://3dmol.org/build/3Dmol-min.js"></script>

<div class="rg-viewer-shell">
  <div class="rg-viewer-nav">
    <div>
      <div class="rg-viewer-title" id="rg-3dmol-title"></div>
      <div class="rg-viewer-sub" id="rg-3dmol-sub"></div>
    </div>
    <button class="rg-viewer-btn" id="rg-3dmol-next">Next view</button>
  </div>
  <div id="rg-3dmol"></div>
</div>

<script>
const rgStructures = __STRUCTURES__;
let rgViewer = null;
let rgIndex = 0;

function rgHeader(i) {
  const item = rgStructures[i] || {};
  document.getElementById('rg-3dmol-title').textContent = item.label || `View ${i + 1}`;
  document.getElementById('rg-3dmol-sub').textContent = item.subtitle || '';
}

function rgStyleResidues(viewer, sites, style) {
  if (!Array.isArray(sites)) return;
  for (const site of sites) {
    if (!site || !site.chain || !site.resi) continue;
    viewer.addStyle({ chain: site.chain, resi: site.resi }, style);
  }
}

function rgRender(i) {
  const container = document.getElementById('rg-3dmol');
  const item = rgStructures[i] || {};
  rgHeader(i);

  if (!rgViewer) {
    rgViewer = $3Dmol.createViewer(container, { backgroundColor: 'white' });
  }
  rgViewer.clear();
  rgViewer.addModel(item.pdb_text || '', 'pdb');

  rgViewer.setStyle({ chain: 'A' }, { cartoon: { color: '#5F7C8A' } });
  rgViewer.addStyle({ chain: 'B' }, { cartoon: { color: '#D17A6A' } });
  rgViewer.addStyle({ chain: 'C' }, { cartoon: { color: '#C8A24A' } });
  rgViewer.addStyle({ hetflag: true }, { stick: { radius: 0.18, colorscheme: 'greenCarbon' } });

  rgStyleResidues(rgViewer, item.patch_sites, {
    stick: { radius: 0.22, color: '#2E86AB' },
    cartoon: { color: '#2E86AB' },
  });

  rgStyleResidues(rgViewer, item.glyco_sites, {
    sphere: { radius: 0.55, color: '#F6C945', opacity: 0.95 },
    stick: { radius: 0.18, color: '#F6C945' },
  });

  if (item.surface_chain) {
    rgViewer.addSurface(
      $3Dmol.SurfaceType.VDW,
      { opacity: 0.13, color: '#8FB9CA' },
      { chain: item.surface_chain }
    );
  }

  rgViewer.zoomTo();
  rgViewer.render();
}

const rgNextButton = document.getElementById('rg-3dmol-next');
if (rgNextButton) {
  rgNextButton.addEventListener('click', () => {
    rgIndex = (rgIndex + 1) % rgStructures.length;
    rgRender(rgIndex);
  });
}

if (rgStructures.length) {
  rgRender(0);
}
</script>
"""
    html = html.replace("__STRUCTURES__", json.dumps(structures))
    html = html.replace("__HEIGHT__", str(int(height)))
    return html


def build_remote_molstar_gallery_html(structures: List[Dict[str, Any]], *, height: int = 640) -> str:
    html = """
<style>
  .rg-mol-wrap {
    font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial;
    font-size: 14px;
  }
  .rg-mol-nav {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 12px;
    margin: 8px 0;
  }
  .rg-mol-title {
    font-weight: 700;
  }
  .rg-mol-sub {
    color: #666;
    font-size: 12px;
    margin-top: 2px;
  }
  .rg-mol-btn {
    appearance: none;
    border: 1px solid #ddd;
    background: #fafafa;
    border-radius: 10px;
    padding: 8px 12px;
    font-weight: 600;
    cursor: pointer;
  }
  .rg-mol-btn:hover {
    background: #f2f2f2;
  }
  #rg-molstar {
    width: 100%;
    height: __HEIGHT__px;
    position: relative;
    border: 1px solid #eee;
    border-radius: 12px;
    overflow: hidden;
  }
  .rg-mol-loading {
    position: absolute;
    inset: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    background: white;
    color: #666;
    font-size: 12px;
    letter-spacing: 0.08em;
    z-index: 5;
  }
</style>

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/molstar@3/build/viewer/molstar.css" />
<script src="https://cdn.jsdelivr.net/npm/molstar@3/build/viewer/molstar.js"></script>

<div class="rg-mol-wrap">
  <div class="rg-mol-nav">
    <div>
      <div class="rg-mol-title" id="rg-molstar-title"></div>
      <div class="rg-mol-sub" id="rg-molstar-sub"></div>
    </div>
    <button class="rg-mol-btn" id="rg-molstar-next">Next binder</button>
  </div>
  <div id="rg-molstar">
    <div id="rg-molstar-loading" class="rg-mol-loading">LOADING...</div>
  </div>
</div>

<script>
const rgMolStructures = __STRUCTURES__;
let rgMolViewer = null;
let rgMolIndex = 0;

function rgMolHeader(i) {
  const item = rgMolStructures[i] || {};
  document.getElementById('rg-molstar-title').textContent = item.label || `Binder ${i + 1}`;
  document.getElementById('rg-molstar-sub').textContent = item.subtitle || '';
}

async function rgMolEnsureViewer() {
  if (rgMolViewer) return rgMolViewer;
  rgMolViewer = await molstar.Viewer.create('rg-molstar', {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowRemoteState: false,
    layoutShowSequence: true,
    layoutShowLog: false,
  });
  return rgMolViewer;
}

async function rgAddHighlight(plugin, structure, expr, name, colorValue, sizeValue) {
  const comp = await plugin.builders.structure.tryCreateComponentFromExpression(structure, expr, name);
  if (!comp) return;
  await plugin.builders.structure.representation.addRepresentation(comp, {
    type: 'ball-and-stick',
    color: 'uniform',
    colorParams: { value: colorValue },
    size: 'uniform',
    sizeParams: { value: sizeValue },
  });
}

async function rgMolLoad(i) {
  const loading = document.getElementById('rg-molstar-loading');
  if (loading) loading.style.display = 'flex';

  rgMolHeader(i);
  const item = rgMolStructures[i] || {};
  const viewer = await rgMolEnsureViewer();
  const plugin = viewer.plugin;

  try {
    if (plugin && plugin.clear) await plugin.clear();
  } catch (e) {
    console.warn('Mol* clear failed:', e);
  }

  try {
    const data = await plugin.builders.data.download({ url: item.url });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);

    let presetApplied = false;
    try {
      if (plugin.builders?.structure?.representation?.applyPreset) {
        await plugin.builders.structure.representation.applyPreset(structure, 'auto');
        presetApplied = true;
      }
    } catch (e) {
      console.warn('Preset load failed:', e);
    }

    if (!presetApplied) {
      const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
      if (polymer) {
        await plugin.builders.structure.representation.addRepresentation(polymer, {
          type: 'cartoon',
          color: 'chain-id',
        });
      }
      const ligands = await plugin.builders.structure.tryCreateComponentStatic(structure, 'ligand');
      if (ligands) {
        await plugin.builders.structure.representation.addRepresentation(ligands, {
          type: 'ball-and-stick',
          color: 'element-symbol',
        });
      }
    }

    const MS = molstar.MolScriptBuilder;
    const sites = Array.isArray(item.glyco_sites) ? item.glyco_sites : [];
    for (const site of sites) {
      const chain = String(site.chain || '').trim();
      const resi = Number(site.resi);
      if (!chain || !Number.isFinite(resi)) continue;

      if (MS.struct?.atomProperty?.macromolecular?.auth_asym_id && MS.struct?.atomProperty?.macromolecular?.auth_seq_id) {
        const expr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), chain]),
          'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resi]),
        });
        await rgAddHighlight(plugin, structure, expr, `glyco-auth-${chain}-${resi}`, 0xF6C945, 0.30);
      }
    }

  } catch (e) {
    console.error('Mol* load failed:', e);
  }

  if (loading) loading.style.display = 'none';
}

const rgMolNext = document.getElementById('rg-molstar-next');
if (rgMolNext) {
  rgMolNext.addEventListener('click', () => {
    rgMolIndex = (rgMolIndex + 1) % rgMolStructures.length;
    rgMolLoad(rgMolIndex);
  });
}

function rgWaitForMolstar(cb, tries = 0) {
  if (typeof molstar !== 'undefined' && molstar.Viewer) {
    cb();
  } else if (tries < 50) {
    setTimeout(() => rgWaitForMolstar(cb, tries + 1), 200);
  } else {
    const loading = document.getElementById('rg-molstar-loading');
    if (loading) loading.textContent = 'ERROR: Mol* viewer failed to load.';
  }
}

if (rgMolStructures.length) {
  rgWaitForMolstar(() => rgMolLoad(0));
}
</script>
"""
    html = html.replace("__STRUCTURES__", json.dumps(structures))
    html = html.replace("__HEIGHT__", str(int(height)))
    return html


def load_target_context(config: Dict[str, Any]) -> Dict[str, Any]:
    repo_default_glycans = clone_entries(REPO_DEFAULT_GLYCANS)
    repo_preset_patch_entries = derive_preset_patch_from_complex(REPO_COMPLEX_PDB, repo_default_glycans)

    if config["target_source"] == "repo":
        target_name, target_bytes = load_pdb_bytes(str(REPO_TARGET_PDB))
        target_is_repo = True
    else:
        fallback_name = config["custom_target_pdb_path"] or "target.pdb"
        target_name, target_bytes = get_pdb_input(fallback_name)
        target_is_repo = False

    if config["glycan_source"] == "repo_default":
        configured_glycan_entries = clone_entries(repo_default_glycans)
    elif config["glycan_source"] == "custom_json":
        configured_glycan_entries = parse_custom_glycan_json(config["custom_glycan_json"])
    else:
        configured_glycan_entries = []

    if config["patch_mode"] == "manual" and config["custom_patch_residues"]:
        preview_patch_entries = parse_manual_patch(config["custom_patch_residues"])
    elif target_is_repo:
        preview_patch_entries = clone_entries(repo_preset_patch_entries)
    else:
        preview_patch_entries = []

    return {
        "repo_default_glycans": repo_default_glycans,
        "repo_preset_patch_entries": repo_preset_patch_entries,
        "target_name": target_name,
        "target_bytes": target_bytes,
        "target_text": target_bytes.decode("utf-8"),
        "target_is_repo": target_is_repo,
        "configured_glycan_entries": configured_glycan_entries,
        "preview_patch_entries": preview_patch_entries,
    }


def prepare_patch(config: Dict[str, Any], target_state: Dict[str, Any]) -> Dict[str, Any]:
    scan_passed: List[str] = []
    hotspot_sites: List[str] = []
    hotspots_pdb_url = None

    effective_patch_mode = config["patch_mode"]
    if effective_patch_mode == "preset" and not target_state["target_is_repo"]:
        effective_patch_mode = "scan"
        print("Preset patch is only defined for the bundled rEPO target. Falling back to scan mode for this target.")

    effective_glycan_source = config["glycan_source"]
    configured_glycan_entries = clone_entries(target_state["configured_glycan_entries"])
    if effective_glycan_source == "repo_default" and not target_state["target_is_repo"]:
        effective_glycan_source = "custom_json" if configured_glycan_entries else "scan"
        print(f"Repo-default glycan map is only defined for the bundled rEPO target. Falling back to {effective_glycan_source}.")

    target_session_uuid = None

    def ensure_target_session_uuid() -> str:
        nonlocal target_session_uuid
        if target_session_uuid:
            return target_session_uuid
        target_session = create_session_from_pdb_bytes(target_state["target_name"], target_state["target_bytes"])
        target_session_uuid = target_session["session_uuid"]
        print("Target session_uuid:", target_session_uuid)
        return target_session_uuid

    if effective_patch_mode == "scan" or effective_glycan_source == "scan":
        ensure_target_session_uuid()
        scan_job_uuid = submit_job(target_session_uuid, "scan")["job_uuid"]
        scan_meta = wait_job(scan_job_uuid)
        scan_passed = scan_passed_sites(scan_meta)
        print("Scan-passed sites:", scan_passed)

    if effective_glycan_source == "scan":
        configured_glycan_entries = entries_from_site_keys(
            scan_passed,
            glycan=DEFAULT_SCAN_GLYCAN_ID,
            source="scan",
            site_name_prefix="Scan site",
        )

    if effective_patch_mode == "preset":
        design_patch_entries = clone_entries(target_state["repo_preset_patch_entries"])
        patch_source = "preset patch derived from inputs/REPO_complex.pdb"
    elif effective_patch_mode == "manual":
        design_patch_entries = parse_manual_patch(config["custom_patch_residues"])
        patch_source = "manual residue list"
    else:
        if not scan_passed:
            raise RuntimeError('PATCH_MODE="scan" requires scan results, but no scan-passed residues were returned.')
        import requests

        ensemble_params = {
            "ensembleSize": 50,
            "calculateSASA": True,
            "calculateHotspots": True,
            "outputFormat": "PDB",
        }
        selected_scan = {site: DEFAULT_SCAN_GLYCAN_ID for site in scan_passed}
        ensemble_job_uuid = submit_job(
            ensure_target_session_uuid(),
            "ensemble",
            selected_glycans=selected_scan,
            parameters=ensemble_params,
        )["job_uuid"]
        _ = wait_job(ensemble_job_uuid)
        hotspots_pdb_url = job_file_url(ensemble_job_uuid, "hotspots.pdb")
        hotspot_text = requests.get(hotspots_pdb_url, timeout=60).text
        hotspot_sites = parse_hotspots_pdb_text(hotspot_text)
        design_patch_entries = entries_from_site_keys(
            hotspot_sites or scan_passed[:8],
            glycan=DEFAULT_SCAN_GLYCAN_ID,
            source="scan_hotspot",
            site_name_prefix="Hotspot",
        )
        patch_source = "scan-derived hotspot patch"

    return {
        "configured_glycan_entries": clone_entries(configured_glycan_entries),
        "design_patch_entries": clone_entries(design_patch_entries),
        "patch_source": patch_source,
        "scan_passed": scan_passed,
        "hotspot_sites": hotspot_sites,
        "hotspots_pdb_url": hotspots_pdb_url,
        "target_session_uuid": target_session_uuid,
    }


def render_patch_summary(patch_state: Dict[str, Any]) -> str:
    rows = []
    for entry in patch_state["design_patch_entries"]:
        rows.append(
            "<tr>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{escape(entry_label(entry))}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{escape(patch_state['patch_source'])}</td>"
            "</tr>"
        )
    body = "".join(rows) if rows else (
        "<tr><td colspan='2' style='padding:10px; color:#666; border-bottom:1px solid #eee;'>No patch residues were selected.</td></tr>"
    )
    return (
        "<div style='font-family:system-ui, -apple-system, Segoe UI, Roboto, Arial; font-size:14px;'>"
        "<div style='font-weight:700; margin-bottom:8px;'>Resolved design patch</div>"
        "<table style='border-collapse:collapse; width:100%;'>"
        "<thead><tr>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Residue</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Patch source</th>"
        "</tr></thead>"
        f"<tbody>{body}</tbody></table></div>"
    )


def build_target_structures(target_state: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [
        {
            "label": "Protein-only rEPO target",
            "subtitle": "Bundled target used for binder generation",
            "pdb_text": target_state["target_text"],
            "glyco_sites": glycan_entries_to_js(target_state["configured_glycan_entries"]),
            "patch_sites": glycan_entries_to_js(target_state["preview_patch_entries"]),
            "surface_chain": "A",
        },
        {
            "label": "Bundled example complex",
            "subtitle": "Preset patch comes from inputs/REPO_complex.pdb",
            "pdb_text": REPO_COMPLEX_PDB.read_text(),
            "glyco_sites": glycan_entries_to_js(target_state["configured_glycan_entries"]),
            "patch_sites": glycan_entries_to_js(target_state["preview_patch_entries"]),
            "surface_chain": "A",
        },
    ]


def generate_candidate_binders(config: Dict[str, Any], target_state: Dict[str, Any], patch_state: Dict[str, Any]) -> Dict[str, Any]:
    import numpy as np
    import biotite.structure as bs
    from atomworks.io.utils.io_utils import to_cif_file
    from biotite.structure.io.pdb import PDBFile
    from biotite.structure.io.pdb import get_structure as get_pdb_structure
    from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine

    out_dir = Path("binder_design_out").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    target_pdb_path = out_dir / "target_input.pdb"
    target_pdb_path.write_bytes(target_state["target_bytes"])

    pdb = PDBFile.read(str(target_pdb_path))
    aa_all = get_pdb_structure(pdb, model=1)

    def infer_target_chain() -> str:
        chains = [
            str(entry.get("chain", "")).strip()
            for entry in patch_state["design_patch_entries"] + patch_state["configured_glycan_entries"]
            if str(entry.get("chain", "")).strip()
        ]
        if chains:
            ranked = sorted(set(chains), key=lambda chain: (-chains.count(chain), chain))
            return ranked[0]
        chain_ids = sorted(set(aa_all.chain_id.tolist()))
        return chain_ids[0] if chain_ids else "A"

    target_chain = infer_target_chain()
    chain_mask = aa_all.chain_id == target_chain
    aa_mask = bs.filter_amino_acids(aa_all)
    target_chain_atoms = aa_all[chain_mask & aa_mask]
    if len(target_chain_atoms) == 0:
        raise ValueError(f"No amino-acid atoms found for chain {target_chain!r} in target PDB.")

    old_res_ids = np.sort(np.unique(target_chain_atoms.res_id))
    resid_map = {int(old): int(i + 1) for i, old in enumerate(old_res_ids)}
    target_chain_atoms = target_chain_atoms.copy()
    target_chain_atoms.res_id = np.array(
        [resid_map[int(r)] for r in target_chain_atoms.res_id],
        dtype=target_chain_atoms.res_id.dtype,
    )

    target_seq_len = int(len(np.unique(target_chain_atoms.res_id)))
    print(f"Target chain: {target_chain} ({target_seq_len} residues after renumbering)")

    normalized_target_cif = out_dir / "target_normalized.cif"
    to_cif_file(target_chain_atoms, str(normalized_target_cif))
    print("Wrote:", normalized_target_cif)

    def evenly_spaced(items: List[int], k: int) -> List[int]:
        if not items:
            return []
        k = int(max(1, min(int(k), len(items))))
        idx = np.linspace(0, len(items) - 1, k, dtype=int)
        return [items[i] for i in idx]

    def atom_for_res(res_id: int) -> str:
        mask = target_chain_atoms.res_id == res_id
        if not np.any(mask):
            return "CB"
        res_name = str(target_chain_atoms.res_name[mask][0])
        return "CA" if res_name.upper() == "GLY" else "CB"

    raw_patch_resids = sorted(
        {
            int(entry["resi"])
            for entry in patch_state["design_patch_entries"]
            if str(entry.get("chain", target_chain)).strip() == target_chain and int(entry["resi"]) in resid_map
        }
    )
    mapped_patch_resids = sorted({resid_map[int(old)] for old in raw_patch_resids if int(old) in resid_map})
    if mapped_patch_resids:
        chosen_hotspots = evenly_spaced(mapped_patch_resids, 5)
        hotspot_source = "selected design patch"
    else:
        chosen_hotspots = evenly_spaced(list(range(1, target_seq_len + 1)), 5)
        hotspot_source = "fallback evenly spaced patch"

    hotspot_pool = list(mapped_patch_resids) if mapped_patch_resids else list(range(1, target_seq_len + 1))

    def pick_hotspots(pool: List[int], binder_idx: int) -> List[int]:
        if not pool:
            return []
        pool = sorted(set(int(x) for x in pool))
        k = int(max(1, min(5, len(pool))))
        step = max(1, len(pool) // max(1, int(config["n_binders"])))
        offset = (int(binder_idx) * step) % len(pool)
        doubled = pool + pool
        picked = doubled[offset: offset + k]
        return sorted(set(int(x) for x in picked))

    binder_min = max(1, int(config["binder_length"]) - 10)
    binder_max = int(config["binder_length"]) + 10

    base_binder_spec = {
        "dialect": 2,
        "infer_ori_strategy": "hotspots",
        "input": str(normalized_target_cif),
        "contig": f"{binder_min}-{binder_max},/0,{target_chain}1-{target_seq_len}",
        "is_non_loopy": True,
    }

    rfd3_input_spec = {}
    for i in range(1, int(config["n_binders"]) + 1):
        hs = pick_hotspots(hotspot_pool, binder_idx=i)
        sel = {f"{target_chain}{int(r)}": atom_for_res(int(r)) for r in hs}
        rfd3_input_spec[f"binder_design_{i:05d}"] = {**base_binder_spec, "select_hotspots": sel}
        print(f"Binder {i} select_hotspots:", sel)

    input_json_path = out_dir / "rfd3_binder_input.json"
    input_json_path.write_text(json.dumps(rfd3_input_spec, indent=2) + "\n")
    print("Patch source:", hotspot_source)
    print("Patch residues in input numbering:", raw_patch_resids)
    print("Wrote:", input_json_path)

    engine_config = RFD3InferenceConfig(diffusion_batch_size=1, seed=0, verbose=True)
    try:
        engine_config.inference_sampler["step_scale"] = 3.0
        engine_config.inference_sampler["gamma_0"] = 0.2
    except Exception:
        pass
    engine = RFD3InferenceEngine(**engine_config.__dict__)
    outputs = engine.run(inputs=str(input_json_path), out_dir=None, n_batches=1)

    named_results = []
    for key in sorted(outputs.keys()):
        values = list(outputs[key]) if isinstance(outputs[key], (list, tuple)) else [outputs[key]]
        values = [item for item in values if hasattr(item, "atom_array")]
        if values:
            named_results.append(values[0])
    flat_results = named_results[: int(config["n_binders"])]
    if not flat_results:
        raise RuntimeError("RFD3 did not return any results with atom_array.")

    def normalize_complex_target_chain_ids(atom_array: bs.AtomArray, *, target_len: int, desired_chain: str = "A"):
        try:
            aa = atom_array[bs.filter_amino_acids(atom_array)]
        except Exception:
            aa = atom_array
        chain_ids = sorted(set(getattr(aa, "chain_id", []).tolist())) if len(aa) else []
        if not chain_ids:
            return atom_array, {}
        chain_counts: Dict[str, int] = {}
        for chain in chain_ids:
            mask = aa.chain_id == chain
            chain_counts[chain] = int(len(np.unique(aa.res_id[mask])))

        def rank(chain: str):
            count = chain_counts.get(chain, 0)
            return (abs(int(count) - int(target_len)), -int(count), str(chain))

        target_old = sorted(chain_ids, key=rank)[0]
        letters = [chr(i) for i in range(ord("A"), ord("Z") + 1)]
        pool = [letter for letter in letters if letter != desired_chain]
        mapping: Dict[str, str] = {str(target_old): str(desired_chain)}
        idx = 0
        for chain in chain_ids:
            if str(chain) == str(target_old):
                continue
            mapping[str(chain)] = pool[idx] if idx < len(pool) else str(chain)
            idx += 1
        out = atom_array.copy()
        out.chain_id = np.array([mapping.get(str(chain), str(chain)) for chain in out.chain_id], dtype=out.chain_id.dtype)
        return out, mapping

    generated_complex_pdb_paths = []
    pdb_cls = PDBFile
    for i, result in enumerate(flat_results, start=1):
        complex_atoms = result.atom_array
        complex_atoms, chain_map = normalize_complex_target_chain_ids(complex_atoms, target_len=target_seq_len, desired_chain="A")
        if chain_map:
            print(f"[{i}] Chain renaming:", chain_map)
        complex_pdb_path = out_dir / f"complex_generated_{i:05d}.pdb"
        pdb_out = pdb_cls()
        pdb_out.set_structure(complex_atoms)
        pdb_out.write(str(complex_pdb_path))
        generated_complex_pdb_paths.append(str(complex_pdb_path))
        print(f"[{i}/{len(flat_results)}] Saved complex PDB:", complex_pdb_path)

    if generated_complex_pdb_paths:
        legacy_path = Path("complex_generated.pdb").resolve()
        legacy_path.write_bytes(Path(generated_complex_pdb_paths[0]).read_bytes())

    generated_complex_glycan_entries = []
    for entry in patch_state["configured_glycan_entries"]:
        old_resi = int(entry["resi"])
        if old_resi not in resid_map:
            continue
        new_entry = dict(entry)
        new_entry["chain"] = "A"
        new_entry["resi"] = int(resid_map[old_resi])
        generated_complex_glycan_entries.append(new_entry)

    generated_complex_patch_entries = []
    for entry in patch_state["design_patch_entries"]:
        old_resi = int(entry["resi"])
        if old_resi not in resid_map:
            continue
        new_entry = dict(entry)
        new_entry["chain"] = "A"
        new_entry["resi"] = int(resid_map[old_resi])
        generated_complex_patch_entries.append(new_entry)

    return {
        "target_chain": target_chain,
        "target_seq_len": target_seq_len,
        "resid_map": resid_map,
        "generated_complex_pdb_paths": generated_complex_pdb_paths,
        "generated_complex_glycan_entries": generated_complex_glycan_entries,
        "generated_complex_patch_entries": generated_complex_patch_entries,
    }


def build_protein_only_structures(generation_state: Dict[str, Any]) -> List[Dict[str, Any]]:
    structures = []
    for idx, path in enumerate(generation_state["generated_complex_pdb_paths"], start=1):
        structures.append(
            {
                "label": f"Candidate binder {idx}",
                "subtitle": "Protein-only complex before glyco-aware filtering",
                "pdb_text": Path(path).read_text(),
                "glyco_sites": [],
                "patch_sites": glycan_entries_to_js(generation_state["generated_complex_patch_entries"]),
                "surface_chain": "A",
            }
        )
    return structures


def run_compatibility_analysis(
    config: Dict[str, Any],
    patch_state: Dict[str, Any],
    generation_state: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    from pathlib import Path

    def load_many_complex_pdbs() -> List[Tuple[str, bytes]]:
        if generation_state and generation_state.get("generated_complex_pdb_paths"):
            paths = [Path(path) for path in generation_state["generated_complex_pdb_paths"]]
            paths = [path for path in paths if path.exists() and path.is_file()]
            if paths:
                return [(path.name, path.read_bytes()) for path in paths]
        fallback_name, fallback_bytes = get_pdb_input("complex_generated.pdb")
        return [(fallback_name, fallback_bytes)]

    def per_site_pass(job_meta: dict, expected_sites: List[str]) -> Dict[str, bool]:
        items = ((job_meta.get("results") or {}).get("results")) or []
        by_site: Dict[str, List[bool]] = {}
        for item in items:
            site = item.get("residue")
            if not site:
                continue
            by_site.setdefault(site, []).append(item.get("clash_solved") is True)
        return {site: any(by_site.get(site, [])) for site in expected_sites}

    def run_opt_one(session_uuid: str, selected_glycans: Dict[str, str], *, rotamer: bool) -> Tuple[str, dict]:
        params = {
            "populationSize": 64,
            "maxGenerations": 10,
            "outputFormat": "PDB",
            "enableRotamerScan": bool(rotamer),
        }
        job_uuid = submit_job(
            session_uuid,
            "optimization",
            selected_glycans=selected_glycans,
            parameters=params,
        )["job_uuid"]
        meta = wait_job(job_uuid)
        return job_uuid, meta

    complex_inputs = load_many_complex_pdbs()
    binder_results = []
    generated_names = {
        Path(path).name for path in (generation_state or {}).get("generated_complex_pdb_paths", [])
    }

    for binder_idx, (complex_name, complex_bytes) in enumerate(complex_inputs, start=1):
        print(f"\n=== Binder {binder_idx}/{len(complex_inputs)}: {complex_name} ===")
        complex_session = create_session_from_pdb_bytes(complex_name, complex_bytes)
        complex_session_uuid = complex_session["session_uuid"]

        is_generated = complex_name in generated_names
        active_entries = clone_entries(
            generation_state["generated_complex_glycan_entries"]
            if is_generated and generation_state and generation_state.get("generated_complex_glycan_entries")
            else patch_state["configured_glycan_entries"]
        )
        site_lookup = {format_site_key(entry["chain"], int(entry["resi"])): entry_label(entry) for entry in active_entries}

        scan_sites = []
        if config["use_complex_scan"]:
            complex_scan_uuid = submit_job(complex_session_uuid, "scan")["job_uuid"]
            complex_scan_meta = wait_job(complex_scan_uuid)
            scan_sites = scan_passed_sites(complex_scan_meta)
            print("Complex scan-passed sites:", scan_sites)

        if config["glycan_source"] == "scan" or not active_entries:
            sites_sorted = sorted(scan_sites, key=site_sort_key)
            selected_complex = {site: DEFAULT_SCAN_GLYCAN_ID for site in sites_sorted}
        else:
            sites_sorted = sorted(glycan_entries_to_selected_map(active_entries).keys(), key=site_sort_key)
            selected_complex = glycan_entries_to_selected_map(active_entries)

        if not sites_sorted:
            raise RuntimeError("No glycosylation sites available for compatibility analysis.")

        job_off, meta_off = run_opt_one(complex_session_uuid, selected_complex, rotamer=False)
        url_off = job_file_url(job_off, "output.pdb")
        pass_off = per_site_pass(meta_off, sites_sorted)
        clashed_off = [site for site in sites_sorted if not pass_off.get(site, False)]

        job_on = None
        url_on = None
        clashed_on: List[str] = []
        if clashed_off:
            job_on, meta_on = run_opt_one(complex_session_uuid, selected_complex, rotamer=True)
            url_on = job_file_url(job_on, "output.pdb")
            pass_on = per_site_pass(meta_on, sites_sorted)
            clashed_on = [site for site in sites_sorted if not pass_on.get(site, False)]

        if clashed_on:
            glyco_result = "Fail"
            reason = "Persistent glycan clash after rotamer-enabled optimization."
            recommendation = "Discard"
            final_url = url_on or url_off
            blocking_sites = clashed_on
        elif clashed_off:
            glyco_result = "Borderline"
            reason = "Resolvable only after rotamer-enabled local adjustment."
            recommendation = "Review"
            final_url = url_on or url_off
            blocking_sites = clashed_off
        else:
            glyco_result = "Pass"
            reason = "No severe glycan obstruction detected."
            recommendation = "Keep"
            final_url = url_off
            blocking_sites = []

        binder_results.append(
            {
                "binder_idx": binder_idx,
                "protein_only": "Yes",
                "glyco_result": glyco_result,
                "reason": reason,
                "recommendation": recommendation,
                "blocking_labels": [site_lookup.get(site, site) for site in blocking_sites],
                "final_url": final_url,
                "glyco_sites_js": glycan_entries_to_js(active_entries),
            }
        )

    remote_structures = [
        {
            "label": f"Binder {item['binder_idx']}",
            "subtitle": f"{item['glyco_result']} • {item['recommendation']} • {item['reason']}",
            "url": item["final_url"],
            "glyco_sites": item["glyco_sites_js"],
        }
        for item in binder_results
    ]

    rows = []
    for item in binder_results:
        blocking = ", ".join(item["blocking_labels"]) if item["blocking_labels"] else "None"
        rows.append(
            "<tr>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>B{item['binder_idx']}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{item['protein_only']}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{item['glyco_result']}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee;'>{escape(blocking)}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee; white-space:nowrap;'>{item['recommendation']}</td>"
            f"<td style='padding:8px; border-bottom:1px solid #eee;'><a href='{item['final_url']}' target='_blank'>final output.pdb</a></td>"
            "</tr>"
        )

    pass_count = sum(1 for item in binder_results if item["glyco_result"] == "Pass")
    borderline_count = sum(1 for item in binder_results if item["glyco_result"] == "Borderline")
    fail_count = sum(1 for item in binder_results if item["glyco_result"] == "Fail")
    summary_html = (
        "<div style='font-family:system-ui, -apple-system, Segoe UI, Roboto, Arial; font-size:14px;'>"
        "<div style='font-weight:700; margin-bottom:10px;'>Per-binder ReGlyco results</div>"
        "<table style='border-collapse:collapse; width:100%;'>"
        "<thead><tr>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Binder</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Protein-only looks OK</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Glyco-aware result</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Blocking site(s)</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Recommendation</th>"
        "<th style='text-align:left; border-bottom:1px solid #ddd; padding:8px;'>Final structure</th>"
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
        f"<div style='margin-top:14px; padding:12px 14px; border:1px solid #eee; border-radius:12px; background:#fafafa;'>"
        f"{len(binder_results)} candidate binders were evaluated. "
        f"{pass_count} passed immediately, {borderline_count} were borderline after local adjustment, and {fail_count} failed due to persistent glycan incompatibility."
        "</div></div>"
    )
    return {
        "binder_results": binder_results,
        "remote_structures": remote_structures,
        "summary_html": summary_html,
    }
