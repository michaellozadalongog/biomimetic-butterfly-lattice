#!/usr/bin/env python3
"""Build watertight STLs of the Morpho-inspired lattice coupons and demo tile.

Parts are a union of axis-aligned boxes so every wall is vertical. They are
meant to print flat on a Bambu Lab X1 Carbon Cool Plate in PLA, basal skin
down, with no supports.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import trimesh
from trimesh.creation import box

from geometry import (
    ALL_VARIANTS,
    BED_MARGIN_MM,
    BED_SIZE_MM,
    LatticeParams,
    analytical_volume_mm3,
    lattice_region_relative_density,
)


def _box(xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float) -> trimesh.Trimesh:
    extents = np.array([xmax - xmin, ymax - ymin, zmax - zmin], dtype=float)
    if np.any(extents <= 0.0):
        raise ValueError(f"non-positive box extents: {extents}")
    mesh = box(extents=extents)
    center = np.array([0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)])
    mesh.apply_translation(center)
    return mesh



def _solid_intervals(limit: float, blocks: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Merge [a,b] blocks and clip to [0, limit]."""
    items = []
    for a, b in blocks:
        a, b = max(0.0, a), min(limit, b)
        if b > a + 1e-9:
            items.append((a, b))
    if not items:
        return []
    items.sort()
    merged = [items[0]]
    for a, b in items[1:]:
        if a <= merged[-1][1] + 1e-9:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def _open_intervals(limit: float, solid: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Complement of merged solid intervals inside [0, limit]."""
    solid = _solid_intervals(limit, solid)
    open_ = []
    cursor = 0.0
    for a, b in solid:
        if a > cursor + 1e-9:
            open_.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < limit - 1e-9:
        open_.append((cursor, limit))
    return open_


def cavity_boxes(p: LatticeParams) -> list[trimesh.Trimesh]:
    """Interior window cavities to subtract from the solid envelope.

    Each cavity runs between the two skins, so the upper and lower laminae
    stay intact and the ridges/cross-ribs/pads remain as the core.
    """
    L, W, H, t = p.length, p.width, p.height, p.t_skin
    x_solid = [(0.0, p.end_pad), (L - p.end_pad, L)]
    if p.load_pad > 0.0:
        mid0 = L / 2.0 - p.load_pad / 2.0
        mid1 = L / 2.0 + p.load_pad / 2.0
        x_solid.append((mid0, mid1))
    for c in p.cross_centers():
        x_solid.append((c - p.cross_width / 2.0, c + p.cross_width / 2.0))

    y_solid = [(0.0, p.frame_width), (W - p.frame_width, W)]
    for c in p.ridge_centers():
        y_solid.append((c - p.ridge_width / 2.0, c + p.ridge_width / 2.0))

    cavities = []
    # Expand each cavity by a tiny amount in x,y so subtracted faces are not
    # exactly coincident with rib faces (helps the boolean stay watertight).
    eps = 0.0
    for x0, x1 in _open_intervals(L, x_solid):
        for y0, y1 in _open_intervals(W, y_solid):
            if (x1 - x0) < 0.4 or (y1 - y0) < 0.4:
                continue
            cavities.append(_box(x0 + eps, x1 - eps, y0 + eps, y1 - eps, t, H - t))
    return cavities


def build_mesh(p: LatticeParams) -> trimesh.Trimesh:
    L, W, H = p.length, p.width, p.height
    envelope = _box(0.0, L, 0.0, W, 0.0, H)
    if p.kind == "solid":
        mesh = envelope
    else:
        cavities = cavity_boxes(p)
        if not cavities:
            mesh = envelope
        else:
            # Subtract all windows from the solid plate of the same footprint.
            void = trimesh.boolean.union(cavities, engine="manifold") if len(cavities) > 1 else cavities[0]
            mesh = trimesh.boolean.difference([envelope, void], engine="manifold")
    if isinstance(mesh, list):
        raise RuntimeError(f"boolean returned {len(mesh)} bodies for {p.name}")
    # Do not run process(validate=True) / fix_normals: that flips the
    # inward cavity shells to positive volume and fills the windows.
    mesh.merge_vertices()
    mesh.update_faces(mesh.unique_faces())
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    if mesh.volume < 0:
        mesh.invert()
    mesh.metadata["name"] = p.name
    return mesh


def export_variant(p: LatticeParams, stl_dir: Path, meta_dir: Path) -> dict:
    print(f"Building {p.name} ...", flush=True)
    mesh = build_mesh(p)
    stl_path = stl_dir / f"{p.name}.stl"
    mesh.export(stl_path, file_type="stl")

    vol_mesh = float(mesh.volume)
    vol_anal = analytical_volume_mm3(p)
    rel_env = vol_anal / p.envelope_volume_mm3
    rel_lattice = lattice_region_relative_density(p)
    bounds = mesh.bounds
    extents = (bounds[1] - bounds[0]).tolist()
    watertight = bool(mesh.is_watertight)
    rec = {
        "name": p.name,
        "kind": p.kind,
        "params": p.to_dict(),
        "stl": str(stl_path.relative_to(stl_dir.parent)),
        "mesh_volume_mm3": vol_mesh,
        "analytical_volume_mm3": vol_anal,
        "volume_rel_error": (vol_mesh - vol_anal) / vol_anal if vol_anal else None,
        "relative_density_envelope": rel_env,
        "relative_density_lattice_region": rel_lattice,
        "extents_mm": extents,
        "bounds_mm": bounds.tolist(),
        "is_watertight": watertight,
        "n_faces": int(len(mesh.faces)),
        "n_vertices": int(len(mesh.vertices)),
        "fits_cool_plate": all(e <= BED_SIZE_MM - 2 * BED_MARGIN_MM for e in extents[:2]),
    }
    meta_path = meta_dir / f"{p.name}.json"
    meta_path.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(
        f"  STL {stl_path.name}: V={vol_anal:.1f} mm^3 (mesh {vol_mesh:.1f}), "
        f"ρ*/ρs envelope={rel_env:.3f}, lattice region={rel_lattice:.3f}, "
        f"watertight={watertight}, faces={rec['n_faces']}",
        flush=True,
    )
    if not watertight:
        print(f"  WARNING: {p.name} mesh is not watertight", flush=True)
    if rec["volume_rel_error"] is not None and abs(rec["volume_rel_error"]) > 0.03:
        print(f"  WARNING: mesh vs analytical volume differ by {rec['volume_rel_error']*100:.1f}%", flush=True)
    return rec


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--stl-dir", type=Path, default=root / "stl")
    parser.add_argument("--results-dir", type=Path, default=root / "results")
    args = parser.parse_args(argv)

    args.stl_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    records = [export_variant(p, args.stl_dir, args.results_dir) for p in ALL_VARIANTS]
    summary = {"variants": records, "bed_size_mm": BED_SIZE_MM, "bed_margin_mm": BED_MARGIN_MM}
    (args.results_dir / "geometry_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("Wrote", args.results_dir / "geometry_summary.json")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
