#!/usr/bin/env python3
"""Figures: software-rendered STLs, comparison plots, schematic, worksheet."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import trimesh
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle
from PIL import Image, ImageDraw, ImageFont

from geometry import COUPON_LOW, COUPON_MID, DEMO_MID, LatticeParams

# Honest, high-contrast figures on a light background (portfolio-friendly).
NAVY = "#1b3a4b"
TEAL = "#1d7a6a"
AMBER = "#c4811a"
SLATE = "#4a5560"
RUST = "#b4533a"
CREAM = "#f6f1e8"
INK = "#1a1a1a"


def _rot(az_deg: float, el_deg: float) -> np.ndarray:
    az = np.deg2rad(az_deg)
    el = np.deg2rad(el_deg)
    Rz = np.array(
        [[np.cos(az), -np.sin(az), 0.0], [np.sin(az), np.cos(az), 0.0], [0.0, 0.0, 1.0]]
    )
    Rx = np.array(
        [[1.0, 0.0, 0.0], [0.0, np.cos(el), -np.sin(el)], [0.0, np.sin(el), np.cos(el)]]
    )
    return Rx @ Rz


def render_mesh(
    mesh: trimesh.Trimesh,
    path: Path,
    size: int = 1600,
    az: float = 38.0,
    el: float = 28.0,
    color: tuple[float, float, float] = (0.18, 0.52, 0.46),
    bg: tuple[int, int, int] = (250, 247, 242),
    pad: float = 0.10,
    light: np.ndarray | None = None,
) -> Path:
    """Lambert + ambient software render with a numpy z-buffer."""
    verts = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    center = 0.5 * (verts.min(axis=0) + verts.max(axis=0))
    v = (verts - center) @ _rot(az, el).T
    tris = v[faces]
    # Face normals in camera space (z out of the page after this mapping:
    # we will use x, y for the image and -z as depth so nearer is smaller z).
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    area = np.linalg.norm(n, axis=1)
    keep = area > 1e-12
    tris = tris[keep]
    n = n[keep]
    n = n / np.linalg.norm(n, axis=1)[:, None]
    # Camera looks toward -Z in this rotated frame (x right, y up, z toward camera
    # after we flip): after Rx@Rz, z is "up-ish". Re-map to screen: X, Y, depth.
    # Our rotation leaves z as the original vertical tilted; we want:
    # screen_x = v_x, screen_y = v_z, depth = v_y (into the page).
    # After az around z and el around x, the "view" is looking roughly along -y
    # if we then take (x, z) as the screen. depth = -y (larger y is closer? )
    # Let's use a proper look: after R, use (x, y) as screen and z as depth
    # by first applying a basis change so camera looks along +Y then swap.

    # Simpler: after R, coordinates are already rotated. Treat them as
    # camera space with +Z toward the camera, +Y up, +X right — but Rx(el)@Rz(az)
    # applied to (x,y,z) world does not put the camera on +Z. Apply an extra
    # world-to-camera: map (x,y,z) -> (x, z, -y) so we look along +Y.
    screen = np.empty_like(tris)
    screen[:, :, 0] = tris[:, :, 0]
    screen[:, :, 1] = tris[:, :, 2]
    # Depth = rotated y. After Rx(el), larger world-z (the top skin) has
    # smaller y', so this keeps the lid closer than the cell floors.
    screen[:, :, 2] = tris[:, :, 1]
    # Normals in the same basis; camera looks toward -Y in this frame.
    n_s = np.column_stack([n[:, 0], n[:, 2], -n[:, 1]])

    # Back-face cull (normal toward camera, +z).
    front = n_s[:, 2] > 0
    # Keep some back faces slightly for thin features? Cull for cleanliness.
    screen = screen[front]
    n_s = n_s[front]
    if screen.shape[0] == 0:
        raise RuntimeError("no front-facing triangles")

    if light is None:
        light = np.array([0.35, 0.45, 0.82], dtype=float)
        light /= np.linalg.norm(light)
    ndotl = np.clip(n_s @ light, 0.0, 1.0)
    ambient = 0.28
    shade = ambient + 0.72 * ndotl  # (F,)

    # Bounds.
    xy = screen[:, :, :2].reshape(-1, 2)
    lo = xy.min(axis=0)
    hi = xy.max(axis=0)
    span = (hi - lo).max()
    mid = 0.5 * (lo + hi)
    lo = mid - 0.5 * span * (1.0 + pad)
    hi = mid + 0.5 * span * (1.0 + pad)
    # Image buffer.
    img = np.zeros((size, size, 3), dtype=np.float64)
    img[:] = np.array(bg, dtype=float) / 255.0
    zbuf = np.full((size, size), np.inf, dtype=np.float64)

    def to_px(pts):
        u = (pts[:, 0] - lo[0]) / (hi[0] - lo[0]) * (size - 1)
        v = (1.0 - (pts[:, 1] - lo[1]) / (hi[1] - lo[1])) * (size - 1)
        return u, v

    col = np.array(color, dtype=float)
    # Paint farthest-first is unnecessary with a z-buffer; loop faces.
    # For speed, skip tiny faces and use a simple barycentric fill.
    for tri, sh in zip(screen, shade):
        u, v = to_px(tri[:, :2])
        z = tri[:, 2]
        u0, u1, u2 = u
        v0, v1, v2 = v
        umin = max(int(np.floor(min(u0, u1, u2))), 0)
        umax = min(int(np.ceil(max(u0, u1, u2))), size - 1)
        vmin = max(int(np.floor(min(v0, v1, v2))), 0)
        vmax = min(int(np.ceil(max(v0, v1, v2))), size - 1)
        if umax < umin or vmax < vmin:
            continue
        denom = (v1 - v2) * (u0 - u2) + (u2 - u1) * (v0 - v2)
        if abs(denom) < 1e-12:
            continue
        us = np.arange(umin, umax + 1)
        vs = np.arange(vmin, vmax + 1)
        UU, VV = np.meshgrid(us, vs)
        w0 = ((v1 - v2) * (UU - u2) + (u2 - u1) * (VV - v2)) / denom
        w1 = ((v2 - v0) * (UU - u2) + (u0 - u2) * (VV - v2)) / denom
        w2 = 1.0 - w0 - w1
        mask = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
        if not mask.any():
            continue
        zpix = w0 * z[0] + w1 * z[1] + w2 * z[2]
        view = zbuf[vmin : vmax + 1, umin : umax + 1]
        closer = mask & (zpix < view)
        if not closer.any():
            continue
        view[closer] = zpix[closer]
        pix = img[vmin : vmax + 1, umin : umax + 1]
        rgb = np.clip(col * sh, 0.0, 1.0)
        pix[closer] = rgb
        img[vmin : vmax + 1, umin : umax + 1] = pix
        zbuf[vmin : vmax + 1, umin : umax + 1] = view

    # Soft drop shadow: darken below the occupied silhouette a few pixels.
    alpha = np.isfinite(zbuf)
    out = (img * 255.0).clip(0, 255).astype(np.uint8)
    Image.fromarray(out, mode="RGB").save(path)
    return path




def clip_mesh_ymax(mesh: trimesh.Trimesh, y_max: float) -> trimesh.Trimesh:
    """Keep the portion with y <= y_max. The cut is left open so cells show."""
    verts = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    y = verts[:, 1]
    new_faces = []
    extra_verts = []

    def vid(p):
        extra_verts.append(p)
        return len(verts) + len(extra_verts) - 1

    def interp(ia, ib):
        ya, yb = y[ia], y[ib]
        tpar = (y_max - ya) / (yb - ya)
        return verts[ia] + tpar * (verts[ib] - verts[ia])

    for a, b, c in faces:
        inside = [y[a] <= y_max, y[b] <= y_max, y[c] <= y_max]
        n_in = sum(inside)
        ids = [a, b, c]
        if n_in == 3:
            new_faces.append(ids)
            continue
        if n_in == 0:
            continue
        # Clip against the plane; keep the inside polygon (3 or 4 verts).
        pts_ids = []
        for i in range(3):
            i0, i1 = ids[i], ids[(i + 1) % 3]
            in0, in1 = y[i0] <= y_max, y[i1] <= y_max
            if in0:
                pts_ids.append(i0)
            if in0 != in1:
                pts_ids.append(vid(interp(i0, i1)))
        if len(pts_ids) == 3:
            new_faces.append(pts_ids)
        elif len(pts_ids) == 4:
            new_faces.append([pts_ids[0], pts_ids[1], pts_ids[2]])
            new_faces.append([pts_ids[0], pts_ids[2], pts_ids[3]])
    if not new_faces:
        raise RuntimeError("clip removed every face")
    all_v = np.vstack([verts, np.array(extra_verts)]) if extra_verts else verts
    out = trimesh.Trimesh(vertices=all_v, faces=np.array(new_faces, dtype=int), process=False)
    return out



def voxel_surface_mesh(params, pitch: float = 0.25, y_max: float | None = None) -> trimesh.Trimesh:
    """Exposed-face mesh of the solid occupancy. Optional y_max makes a cutaway."""
    from geometry import is_solid
    p = params
    nx = int(np.round(p.length / pitch))
    ny = int(np.round(p.width / pitch))
    nz = int(np.round(p.height / pitch))
    hx, hy, hz = p.length / nx, p.width / ny, p.height / nz
    ix, iy, iz = np.indices((nx, ny, nz))
    cx = (ix + 0.5) * hx
    cy = (iy + 0.5) * hy
    cz = (iz + 0.5) * hz
    solid = is_solid(cx, cy, cz, p)
    if y_max is not None:
        solid &= cy <= y_max
    # Pad for neighbor tests.
    pad = np.pad(solid, 1, mode="constant", constant_values=False)
    faces = []
    verts = []

    def add_face(corners):
        v0 = len(verts)
        verts.extend(corners)
        faces.append((v0, v0 + 1, v0 + 2))
        faces.append((v0, v0 + 2, v0 + 3))

    # For each solid voxel, emit a face if the neighbor is empty.
    # offsets: (axis, sign, corner order in the empty-facing plane)
    dirs = [
        (0, -1),
        (0, +1),
        (1, -1),
        (1, +1),
        (2, -1),
        (2, +1),
    ]
    # Corner templates for a unit cube, then scale.
    for i, j, k in zip(*np.nonzero(solid)):
        x0, y0, z0 = i * hx, j * hy, k * hz
        x1, y1, z1 = x0 + hx, y0 + hy, z0 + hz
        # pad index is +1
        neighbors = {
            (0, -1): (not pad[i, j + 1, k + 1], [(x0, y0, z0), (x0, y0, z1), (x0, y1, z1), (x0, y1, z0)]),
            (0, +1): (not pad[i + 2, j + 1, k + 1], [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)]),
            (1, -1): (not pad[i + 1, j, k + 1], [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)]),
            (1, +1): (not pad[i + 1, j + 2, k + 1], [(x0, y1, z0), (x0, y1, z1), (x1, y1, z1), (x1, y1, z0)]),
            (2, -1): (not pad[i + 1, j + 1, k], [(x0, y0, z0), (x0, y1, z0), (x1, y1, z0), (x1, y0, z0)]),
            (2, +1): (not pad[i + 1, j + 1, k + 2], [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]),
        }
        for exposed, corners in neighbors.values():
            if exposed:
                add_face(corners)
    return trimesh.Trimesh(vertices=np.array(verts, dtype=float), faces=np.array(faces, dtype=int), process=False)


def make_cutaway(mesh: trimesh.Trimesh, length: float, width: float, height: float) -> trimesh.Trimesh:
    # Back-compat wrapper; prefer voxel_surface_mesh with a named LatticeParams.
    return clip_mesh_ymax(mesh, y_max=width * 0.52)


def make_cutaway(mesh: trimesh.Trimesh, length: float, width: float, height: float) -> trimesh.Trimesh:
    """Open a longitudinal cut just past mid-width so the sandwich core is visible."""
    return clip_mesh_ymax(mesh, y_max=width * 0.52)


def load_comparison(results: Path) -> dict:

    return json.loads((results / "comparison.json").read_text(encoding="utf-8"))


def style_axes(ax, title: str, ylabel: str) -> None:
    ax.set_title(title, loc="left", fontsize=12, color=INK, pad=10, fontweight="semibold")
    ax.set_ylabel(ylabel, fontsize=10, color=SLATE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#c5c0b8")
    ax.spines["bottom"].set_color("#c5c0b8")
    ax.tick_params(colors=SLATE, labelsize=9)
    ax.set_facecolor("white")
    ax.grid(axis="y", color="#ece8e1", linewidth=0.8)
    ax.set_axisbelow(True)


def plot_bars(payload: dict, fig_dir: Path) -> None:
    rows = payload["coupons"]
    names = [r["name"].replace("coupon_", "").replace("_", "\n") for r in rows]
    colors = [SLATE, TEAL, "#2f6f8f", AMBER]
    masses = [r["mass"]["mass_g"] for r in rows]
    phis = [r["mass"]["relative_density_envelope"] for r in rows]
    ks = [r["beam_theory"]["k_timoshenko_N_per_mm"] for r in rows]
    kfea = [r.get("fea", {}).get("k_fea_N_per_mm") for r in rows]
    sk = [r["specific_stiffness_N_per_mm_per_g"] for r in rows]
    sp = [r["specific_first_crack_load_N_per_g"] for r in rows]

    def barfig(values, ylabel, title, fname, values2=None, label1=None, label2=None):
        fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=160)
        fig.patch.set_facecolor(CREAM)
        x = np.arange(len(names))
        ax.bar(x, values, color=colors, width=0.62, zorder=3, label=label1)
        if values2 is not None and all(v is not None for v in values2):
            ax.scatter(x, values2, color=INK, zorder=4, s=36, label=label2)
            ax.legend(frameon=False, fontsize=8)
        ax.set_xticks(x, names)
        style_axes(ax, title, ylabel)
        fig.text(
            0.01,
            0.015,
            "Computational estimate only — not physical test data. Assumed isotropic PLA, E = 2.7 GPa.",
            fontsize=7,
            color=SLATE,
        )
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        fig.savefig(fig_dir / fname, facecolor=fig.get_facecolor())
        plt.close(fig)

    barfig(masses, "Mass (g)", "Printed coupon mass (analytical volume × 1.24 g/cm³)", "fig_mass.png")
    barfig(phis, "Relative density  ρ*/ρs", "Envelope relative density (volume / bounding box)", "fig_relative_density.png")
    barfig(
        ks,
        "Flexural stiffness k (N/mm)",
        "3-point-bend stiffness, 100 mm span",
        "fig_stiffness.png",
        values2=kfea if any(v is not None for v in kfea) else None,
        label1="Beam theory (Timoshenko, I(x))",
        label2="Voxel FEA (linear-elastic)",
    )
    barfig(sk, "Specific stiffness  k / mass  (N/mm/g)", "Specific flexural stiffness vs the solid control", "fig_specific_stiffness.png")
    barfig(sp, "Specific first-crack load  P / mass  (N/g)", "Specific first-crack load estimate (assumed 60 MPa)", "fig_specific_strength.png")

    # E_eff / Es vs phi with Gibson-Ashby curves.
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=160)
    fig.patch.set_facecolor(CREAM)
    phi_line = np.linspace(0.05, 1.0, 200)
    ax.plot(phi_line, phi_line**2, color=RUST, lw=1.6, label="Gibson–Ashby open-cell  E*/Es = (ρ*/ρs)²")
    ax.plot(phi_line, phi_line, color=AMBER, lw=1.6, ls="--", label="Stretch-dominated  E*/Es = ρ*/ρs")
    for r, c in zip(rows, colors):
        ax.scatter(
            r["mass"]["relative_density_envelope"],
            r["beam_theory"]["E_eff_over_Es"],
            s=70,
            color=c,
            zorder=3,
            label=r["name"].replace("coupon_", "").replace("_", " "),
        )
    style_axes(ax, "Effective flexural modulus vs relative density", "E_eff / E_s")
    ax.set_xlabel("Envelope relative density  ρ*/ρs", fontsize=10, color=SLATE)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.15)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.text(
        0.01,
        0.015,
        "E_eff is back-calculated from beam-theory k using the solid-section I. "
        "Ribbed plates sit above open-cell foam scaling because the ridges are continuous spars. Computational estimate only.",
        fontsize=7,
        color=SLATE,
        wrap=True,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(fig_dir / "fig_scaling.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_schematic(fig_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 5.0), dpi=170)
    fig.patch.set_facecolor(CREAM)

    # Left: published Morpho-type scale architecture (schematic, not to scale).
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Morpho-type scale (schematic, µm-scale)", loc="left", fontsize=11, color=INK, fontweight="semibold")
    # Lower lamina
    ax.add_patch(Rectangle((1.0, 1.2), 8.0, 0.28, facecolor="#cfc6b6", edgecolor=NAVY, lw=0.8))
    ax.text(5, 0.7, "lower lamina (basal plate)", ha="center", fontsize=8, color=SLATE)
    # Trabeculae
    for x in (2.0, 3.6, 5.2, 6.8, 8.2):
        ax.add_patch(Rectangle((x - 0.08, 1.48), 0.16, 1.15, facecolor="#8aa39a", edgecolor=NAVY, lw=0.4))
    ax.text(1.15, 1.9, "trabeculae", fontsize=7, color=SLATE, rotation=90, va="center")
    # Upper lamina: ridges with stacked lamellae ("Christmas tree") and cross-ribs
    ridge_xs = [2.0, 3.6, 5.2, 6.8, 8.2]
    for x in ridge_xs:
        # ridge stem
        ax.add_patch(Rectangle((x - 0.09, 2.55), 0.18, 2.6, facecolor=TEAL, edgecolor=NAVY, lw=0.5))
        for i, z in enumerate((2.9, 3.35, 3.8, 4.25, 4.7)):
            w = 0.85 - i * 0.08
            ax.add_patch(Rectangle((x - w / 2, z), w, 0.16, facecolor="#4f9a8a", edgecolor=NAVY, lw=0.35))
    # cross-ribs between ridges at upper lamina level
    for a, b in zip(ridge_xs[:-1], ridge_xs[1:]):
        ax.plot([a + 0.1, b - 0.1], [2.7, 2.7], color=NAVY, lw=1.6, solid_capstyle="round")
    ax.text(5, 5.55, "ridges with stacked lamellae", ha="center", fontsize=8, color=SLATE)
    ax.text(5, 2.35, "cross-ribs  +  windows", ha="center", fontsize=8, color=SLATE)
    ax.text(
        5,
        6.55,
        "After Ghiradella; Giraldo & Stavenga 2016;\nDing, Xu & Wang 2009 (Morpho peleides ridge pitch ~1.7–1.8 µm)",
        ha="center",
        fontsize=6.5,
        color=SLATE,
    )

    # Right: printable abstraction.
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Printable coupon abstraction (mm-scale, this repo)", loc="left", fontsize=11, color=INK, fontweight="semibold")
    # basal skin (lower lamina)
    ax.add_patch(Rectangle((1.0, 1.2), 8.0, 0.7, facecolor="#cfc6b6", edgecolor=NAVY, lw=0.8))
    ax.text(5, 0.7, "basal skin  t = 1.0 mm  (lower lamina)", ha="center", fontsize=8, color=SLATE)
    # frames / ridges
    for x, w in ((1.5, 0.45), (3.15, 0.45), (4.85, 0.45), (6.55, 0.45), (8.2, 0.45)):
        ax.add_patch(Rectangle((x, 1.2), w, 3.6, facecolor=TEAL, edgecolor=NAVY, lw=0.6))
    # one cross-rib in front (drawn as a translucent slab)
    ax.add_patch(Rectangle((1.0, 1.2), 8.0, 0.55, facecolor="#2f6f8f", edgecolor=NAVY, lw=0.6, alpha=0.55))
    # upper lamina (top skin)
    ax.add_patch(Rectangle((1.0, 4.3), 8.0, 0.7, facecolor="#cfc6b6", edgecolor=NAVY, lw=0.8, alpha=0.95))
    ax.text(5, 5.25, "top skin  t = 1.0 mm  (upper lamina, bridged)", ha="center", fontsize=8, color=SLATE)
    ax.text(5, 2.05, "cross-rib / ridge core", ha="center", fontsize=7, color="white")
    ax.text(
        5,
        6.45,
        "Walls are vertical; the top skin is a short 3–8 mm bridge.\n"
        "Print flat, basal skin down, no supports (they would be trapped).\n"
        "Lamellae are not printed. Trabeculae = the walls between skins.",
        ha="center",
        fontsize=6.5,
        color=SLATE,
    )

    fig.suptitle(
        "Biology → design map  ·  a scaled concept study, not a replica of chitin armor",
        fontsize=12,
        color=INK,
        x=0.02,
        ha="left",
        y=0.98,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.92))
    fig.savefig(fig_dir / "fig_schematic.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_print_layout(fig_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=160)
    fig.patch.set_facecolor(CREAM)
    bed = 256.0
    ax.set_xlim(-8, bed + 8)
    ax.set_ylim(-8, bed + 8)
    ax.set_aspect("equal")
    ax.add_patch(Rectangle((0, 0), bed, bed, facecolor="white", edgecolor=NAVY, lw=1.4))
    ax.set_title("Suggested X1C Cool Plate layout (mm)", loc="left", fontsize=11, color=INK, fontweight="semibold")
    # 4 coupons in a row near the bottom, demo tile above.
    gap = 6.0
    margin = 8.0
    cw, cl = 25.0, 120.0
    xs = margin
    # Place coupons with length along X, width along Y, stacked in Y.
    y = margin
    labels = ["solid", "low", "mid", "high"]
    cols = [SLATE, TEAL, "#2f6f8f", AMBER]
    for lab, c in zip(labels, cols):
        ax.add_patch(Rectangle((xs, y), cl, cw, facecolor=c, alpha=0.55, edgecolor=NAVY, lw=0.7))
        ax.text(xs + cl / 2, y + cw / 2, lab, ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        y += cw + gap
    # Demo tile to the right
    tile = 100.0
    tx = xs + cl + gap + 4
    ty = margin
    ax.add_patch(Rectangle((tx, ty), tile, tile, facecolor=TEAL, alpha=0.35, edgecolor=NAVY, lw=0.7))
    ax.text(tx + tile / 2, ty + tile / 2, "demo tile\n100 × 100", ha="center", va="center", fontsize=8, color=NAVY)
    ax.text(
        bed / 2,
        bed - 14,
        "256 × 256 mm Cool Plate  ·  all parts flat, basal skin down, no supports\n"
        "Print PLA with the Bambu Studio Cool Plate + PLA preset, 0.20 mm layers, 100% infill, ≥3 walls.",
        ha="center",
        va="top",
        fontsize=7.5,
        color=SLATE,
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_print_layout.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_unit_cell(fig_dir: Path, stl_dir: Path) -> None:
    """Top-view plan of the mid-density lattice region (true geometry)."""
    p = COUPON_MID
    fig, ax = plt.subplots(figsize=(8.4, 3.2), dpi=170)
    fig.patch.set_facecolor(CREAM)
    ax.set_xlim(p.end_pad - 2, p.length - p.end_pad + 2)
    ax.set_ylim(-2, p.width + 2)
    ax.set_aspect("equal")
    ax.add_patch(Rectangle((0, 0), p.length, p.width, facecolor="#efe8dc", edgecolor=NAVY, lw=0.6, alpha=0.5))
    # skin already implied; draw full-height features
    def rib(x0, y0, w, h, c=TEAL):
        ax.add_patch(Rectangle((x0, y0), w, h, facecolor=c, edgecolor=NAVY, lw=0.3, alpha=0.9))

    rib(0, 0, p.length, p.frame_width, "#2f6f8f")
    rib(0, p.width - p.frame_width, p.length, p.frame_width, "#2f6f8f")
    rib(0, 0, p.end_pad, p.width, "#2f6f8f")
    rib(p.length - p.end_pad, 0, p.end_pad, p.width, "#2f6f8f")
    mid0 = p.length / 2 - p.load_pad / 2
    rib(mid0, 0, p.load_pad, p.width, AMBER)
    for c in p.ridge_centers():
        rib(0, c - p.ridge_width / 2, p.length, p.ridge_width, TEAL)
    for c in p.cross_centers():
        rib(c - p.cross_width / 2, 0, p.cross_width, p.width, "#3d8b7a")
    ax.annotate(
        "solid end pad\n(support land)",
        xy=(p.end_pad / 2.0, p.width * 0.72),
        xytext=(22, p.width + 7),
        fontsize=7,
        color=SLATE,
        arrowprops=dict(arrowstyle="->", color=SLATE, lw=0.7),
        ha="left",
    )
    ax.annotate(
        "midspan load land",
        xy=(p.length / 2, p.width),
        xytext=(p.length / 2 + 10, p.width + 6),
        fontsize=7,
        color=SLATE,
        arrowprops=dict(arrowstyle="->", color=SLATE, lw=0.7),
        ha="left",
    )
    ax.set_xlabel("length x (mm)", fontsize=8, color=SLATE)
    ax.set_ylabel("width y (mm)", fontsize=8, color=SLATE)
    style_axes(ax, "Mid-density coupon plan (full-height features in color, skin everywhere)", "width y (mm)")
    ax.set_xlim(0, p.length)
    ax.set_ylim(-8, p.width + 12)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_coupon_plan.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def compose_hero(fig_dir: Path, renders: list[Path], title: str, fname: str) -> Path:
    imgs = [Image.open(p).convert("RGB") for p in renders]
    w = imgs[0].width
    h = imgs[0].height
    n = len(imgs)
    gap = 16
    canvas_w = n * w + (n + 1) * gap
    header = 90
    footer = 54
    canvas = Image.new("RGB", (canvas_w, h + header + footer), (246, 241, 232))
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
        font_s = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16)
    except Exception:
        font = ImageFont.load_default()
        font_s = font
    draw.text((gap, 22), title, fill=(26, 26, 26), font=font)
    draw.text(
        (gap, header + h + 14),
        "Software render of the printable CAD  ·  PLA concept study, not ballistic armor  ·  Michael Lozada-Longog",
        fill=(74, 85, 96),
        font=font_s,
    )
    for i, im in enumerate(imgs):
        canvas.paste(im, (gap + i * (w + gap), header))
    out = fig_dir / fname
    canvas.save(out, quality=95)
    return out



def plot_sections(fig_dir: Path) -> None:
    """True-scale cross-sections through the mid-density coupon (computational geometry)."""
    from geometry import COUPON_LOW, COUPON_MID, COUPON_HIGH, COUPON_SOLID, is_solid

    def section_image(p, x_station, ny=500, nz=160):
        ys = np.linspace(0, p.width, ny)
        zs = np.linspace(0, p.height, nz)
        Y, Z = np.meshgrid(ys, zs, indexing="xy")
        X = np.full_like(Y, x_station)
        return is_solid(X, Y, Z, p), ys, zs

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 6.2), dpi=170, sharex=True)
    fig.patch.set_facecolor(CREAM)
    # A station in the open lattice, not on a cross-rib or pad.
    # Mid coupon: end_pad=15, load pad 56–64, pick x=40.
    variants = [
        (COUPON_LOW, TEAL, "low"),
        (COUPON_MID, "#2f6f8f", "mid"),
        (COUPON_HIGH, AMBER, "high"),
    ]
    x_station = 40.0
    for ax, (p, color, name) in zip(axes, variants):
        solid, ys, zs = section_image(p, x_station)
        ax.imshow(
            solid.astype(float),
            origin="lower",
            extent=(0, p.width, 0, p.height),
            aspect="equal",
            cmap="Blues",
            vmin=0,
            vmax=1.3,
        )
        ax.set_ylabel("z (mm)", fontsize=8, color=SLATE)
        ax.set_title(f"{name}  ·  cross-section at x = {x_station:.0f} mm (two 1 mm skins + ridge core)", loc="left", fontsize=10, color=INK)
        ax.tick_params(colors=SLATE, labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(-0.4, p.height + 0.4)
    axes[-1].set_xlabel("width y (mm)", fontsize=8, color=SLATE)
    fig.suptitle("Sandwich sections  ·  computational geometry, not a photograph of a print", fontsize=11, color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(fig_dir / "fig_sections.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_worksheet(fig_dir: Path) -> None:
    """One-page printable data sheet for the physical tests."""
    fig = plt.figure(figsize=(8.5, 11), dpi=160)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)

    def text(x, y, s, size=10, weight="normal", color=INK, ha="left", va="top"):
        ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, wrap=True)

    text(0.5, 10.55, "Biomimetic lattice  ·  physical test worksheet", 14, "semibold")
    text(0.5, 10.22, "Scaled PLA concept study  —  this is NOT ballistic armor. Computational numbers are estimates only.", 8, color=SLATE)
    text(0.5, 9.92, "Name: ____________________    Date: __________    Printer: Bambu X1C + Cool Plate    Filament: PLA", 9)

    text(0.5, 9.55, "1. Print record", 11, "semibold")
    headers = ["Variant", "Predicted mass (g)", "Measured mass (g)", "Print time", "Notes (warp, stringing, gaps)"]
    y = 9.25
    colx = [0.5, 1.9, 3.4, 4.9, 6.1]
    for x, h in zip(colx, headers):
        text(x, y, h, 7.5, "semibold", SLATE)
    for i, name in enumerate(["coupon_solid", "coupon_lattice_low", "coupon_lattice_mid", "coupon_lattice_high", "panel_demo_mid"]):
        yy = y - 0.32 * (i + 1)
        ax.plot([0.5, 8.0], [yy - 0.22, yy - 0.22], color="#ddd", lw=0.6)
        text(0.5, yy, name.replace("_", " "), 8)

    text(0.5, 7.45, "2. Three-point bend  (span S = 100 mm, supports 10 mm from each end, load on the midspan land)", 11, "semibold")
    text(
        0.5,
        7.18,
        "Supports must sit on the solid end pads. Load the amber midspan land with a dowel or printed roller. "
        "Measure midspan deflection from a fixed camera/ruler. Stop at first crack or when it feels unsafe. Eye protection on.",
        7.5,
        color=SLATE,
    )
    headers = ["Variant", "Mass used (g)", "Load P (N or g)", "Deflection δ (mm)", "k = P/δ (N/mm)", "Fail mode"]
    y = 6.78
    colx = [0.5, 1.9, 3.15, 4.5, 5.9, 7.2]
    for x, h in zip(colx, headers):
        text(x, y, h, 7.5, "semibold", SLATE)
    for i in range(8):
        yy = y - 0.28 * (i + 1)
        ax.plot([0.5, 8.0], [yy - 0.18, yy - 0.18], color="#ddd", lw=0.6)

    text(0.5, 4.25, "3. Drop / impact  (optional, low energy only)", 11, "semibold")
    text(
        0.5,
        3.98,
        "Simply support the coupon on the end pads. Drop a known blunt mass (e.g. 50–200 g steel or 3D-printed weight with a rounded face) "
        "from a measured height onto the midspan land. Start low (2–5 cm) and increase. Do not use sharp projectiles, BBs, or anything that can ricochet. "
        "Impact energy U = m g h. Record whether the part dents, cracks, or stays elastic. This is not an armor test and cannot be compared to NIJ or military ballistic standards.",
        7.5,
        color=SLATE,
    )
    headers = ["Variant", "Mass m (g)", "Height h (cm)", "U = mgh (J)", "Outcome (elastic / dent / crack / split)", "Photo?"]
    y = 3.48
    colx = [0.5, 1.7, 2.8, 4.0, 5.2, 7.5]
    for x, h in zip(colx, headers):
        text(x, y, h, 7.5, "semibold", SLATE)
    for i in range(6):
        yy = y - 0.28 * (i + 1)
        ax.plot([0.5, 8.0], [yy - 0.18, yy - 0.18], color="#ddd", lw=0.6)

    text(0.5, 1.5, "4. After the test", 11, "semibold")
    text(
        0.5,
        1.22,
        "• Compare measured mass to results/comparison.csv. A large miss usually means under-extrusion or a slicer that ignored thin walls.\n"
        "• Compare measured k to the Timoshenko and FEA estimates. As-printed PLA is anisotropic, so a 20–40% miss is not surprising.\n"
        "• Do not report these PLA results as armor performance. If you later test a tougher material, keep the same geometry and protocol.\n"
        "• Safety: PLA shards can be sharp. Deburr, wear eye protection, keep hands out from under the load, and stop if the fixture tips.",
        7.5,
        color=SLATE,
    )
    text(0.5, 0.35, "Repo: biomimetic-butterfly-lattice  ·  Michael Lozada-Longog  ·  HPU Biomedical Engineering", 7.5, color=SLATE)
    fig.savefig(fig_dir / "fig_test_worksheet.png", facecolor="white")
    fig.savefig(fig_dir / "test_worksheet.pdf", facecolor="white")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--stl-dir", type=Path, default=root / "stl")
    parser.add_argument("--results-dir", type=Path, default=root / "results")
    parser.add_argument("--fig-dir", type=Path, default=root / "figures")
    args = parser.parse_args(argv)
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    colors = {
        "coupon_solid": (0.35, 0.40, 0.44),
        "coupon_lattice_low": (0.16, 0.52, 0.45),
        "coupon_lattice_mid": (0.18, 0.44, 0.56),
        "coupon_lattice_high": (0.75, 0.50, 0.12),
        "panel_demo_mid": (0.16, 0.52, 0.45),
    }
    renders = []
    for name, color in colors.items():
        stl = args.stl_dir / f"{name}.stl"
        if not stl.exists():
            print("missing", stl)
            continue
        mesh = trimesh.load_mesh(stl, force="mesh")
        out = args.fig_dir / f"render_{name}.png"
        print("Rendering", name, "...", flush=True)
        # Slightly different camera for the square tile.
        az, el = (42, 30) if "panel" in name else (38, 28)
        render_mesh(mesh, out, size=1400, az=az, el=el, color=color)
        renders.append(out)
        from geometry import ALL_VARIANTS
        params_by_name = {q.name: q for q in ALL_VARIANTS}
        if name != "coupon_solid":
            q = params_by_name[name]
            cut = voxel_surface_mesh(q, pitch=0.25, y_max=q.width * 0.52)
            cut_path = args.fig_dir / f"render_{name}_cutaway.png"
            # Look at the cut face (camera on the +Y side).
            render_mesh(cut, cut_path, size=1400, az=az, el=el, color=color)
        if name == "coupon_lattice_mid":
            render_mesh(
                mesh,
                args.fig_dir / "render_coupon_lattice_mid_hero.png",
                size=1800,
                az=48,
                el=32,
                color=color,
            )
            q = params_by_name[name]
            cut = voxel_surface_mesh(q, pitch=0.20, y_max=q.width * 0.52)
            render_mesh(
                cut,
                args.fig_dir / "render_coupon_lattice_mid_cutaway_hero.png",
                size=1800,
                az=40,
                el=24,
                color=color,
            )

    plot_schematic(args.fig_dir)
    plot_print_layout(args.fig_dir)
    plot_unit_cell(args.fig_dir, args.stl_dir)
    plot_sections(args.fig_dir)
    plot_worksheet(args.fig_dir)
    if (args.results_dir / "comparison.json").exists():
        payload = load_comparison(args.results_dir)
        plot_bars(payload, args.fig_dir)

    # Hero strip of the four coupons if all exist.
    coupon_renders = [
        args.fig_dir / f"render_{n}.png"
        for n in ("coupon_solid", "coupon_lattice_low", "coupon_lattice_mid", "coupon_lattice_high")
    ]
    if all(p.exists() for p in coupon_renders):
        compose_hero(
            args.fig_dir,
            coupon_renders,
            "Morpho-inspired lattice coupons  ·  solid, low, mid, high density",
            "fig_coupon_renders.png",
        )
    hero = args.fig_dir / "render_coupon_lattice_mid_cutaway_hero.png"
    if not hero.exists():
        hero = args.fig_dir / "render_coupon_lattice_mid_hero.png"
    demo = args.fig_dir / "render_panel_demo_mid_cutaway.png"
    if not demo.exists():
        demo = args.fig_dir / "render_panel_demo_mid.png"
    if hero.exists() and demo.exists():
        compose_hero(
            args.fig_dir,
            [hero, demo],
            "Mid-density coupon and 100 mm demo tile (cutaway: two skins + ridge core)",
            "fig_hero.png",
        )
    print("Figures in", args.fig_dir)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
