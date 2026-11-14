#!/usr/bin/env python3
"""Linear-elastic voxel FEA of a 3-point bend on the lattice coupons.

Half-span model with midspan symmetry. 8-node hexahedral elements,
isotropic Hooke's law, small-strain. Results are *computational
estimates* under the same assumed PLA constants as analyze.py.

This is a coarse teaching-scale check on the beam-theory stiffness,
not a converged production FEA (no contact pads, no plasticity, no
print anisotropy, one voxel size).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import cg

from analyze import PLA
from geometry import COUPON_VARIANTS, LatticeParams, is_solid


# Canonical 8-node hex natural coordinates (Abaqus / standard ordering).
XI = np.array(
    [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ]
)
GAUSS = 1.0 / np.sqrt(3.0)
GP = np.array([-GAUSS, GAUSS])


def hooke_matrix(E: float, nu: float) -> np.ndarray:
    """6x6 isotropic stiffness, Voigt order xx,yy,zz,xy,yz,zx."""
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2.0 * mu
    C[0, 1] = C[0, 2] = C[1, 0] = C[1, 2] = C[2, 0] = C[2, 1] = lam
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    return C


def hex8_stiffness(hx: float, hy: float, hz: float, C: np.ndarray) -> np.ndarray:
    """24x24 element stiffness for a rectangular hex of size hx,hy,hz."""
    Ke = np.zeros((24, 24))
    # Physical coords of the 8 nodes, centered at the origin.
    xyz = np.column_stack(
        [
            0.5 * hx * XI[:, 0],
            0.5 * hy * XI[:, 1],
            0.5 * hz * XI[:, 2],
        ]
    )
    for gx in GP:
        for gy in GP:
            for gz in GP:
                xi, eta, zeta = gx, gy, gz
                dN_dxi = np.zeros((8, 3))
                for a, (xa, ya, za) in enumerate(XI):
                    dN_dxi[a, 0] = 0.125 * xa * (1.0 + ya * eta) * (1.0 + za * zeta)
                    dN_dxi[a, 1] = 0.125 * ya * (1.0 + xa * xi) * (1.0 + za * zeta)
                    dN_dxi[a, 2] = 0.125 * za * (1.0 + xa * xi) * (1.0 + ya * eta)
                J = dN_dxi.T @ xyz  # 3x3
                detJ = np.linalg.det(J)
                dN_dx = dN_dxi @ np.linalg.inv(J)
                B = np.zeros((6, 24))
                for a in range(8):
                    i = 3 * a
                    dNdx, dNdy, dNdz = dN_dx[a]
                    B[0, i] = dNdx
                    B[1, i + 1] = dNdy
                    B[2, i + 2] = dNdz
                    B[3, i] = dNdy
                    B[3, i + 1] = dNdx
                    B[4, i + 1] = dNdz
                    B[4, i + 2] = dNdy
                    B[5, i] = dNdz
                    B[5, i + 2] = dNdx
                Ke += (B.T @ C @ B) * detJ
    return Ke


def voxelize_half(p: LatticeParams, h: float) -> tuple[np.ndarray, np.ndarray, dict]:
    """Solid voxel index list for the left half, x in [0, L/2]."""
    Lh = p.length / 2.0
    nx = int(np.round(Lh / h))
    ny = int(np.round(p.width / h))
    nz = int(np.round(p.height / h))
    # Snap so the grid exactly covers the half coupon.
    hx, hy, hz = Lh / nx, p.width / ny, p.height / nz
    ix, iy, iz = np.indices((nx, ny, nz))
    cx = (ix + 0.5) * hx
    cy = (iy + 0.5) * hy
    cz = (iz + 0.5) * hz
    solid = is_solid(cx, cy, cz, p)
    vox = np.column_stack(np.nonzero(solid)).astype(np.int32)
    meta = {
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "hx": hx,
        "hy": hy,
        "hz": hz,
        "n_solid": int(vox.shape[0]),
        "n_vox_total": int(nx * ny * nz),
    }
    return vox, solid, meta


def node_id(i: int, j: int, k: int, ny_n: int, nz_n: int) -> int:
    return (i * ny_n + j) * nz_n + k


def assemble(p: LatticeParams, mat: dict, h: float) -> dict:
    t0 = time.time()
    C = hooke_matrix(mat["E_MPa"], mat["nu"])
    vox, solid, meta = voxelize_half(p, h)
    hx, hy, hz = meta["hx"], meta["hy"], meta["hz"]
    Ke = hex8_stiffness(hx, hy, hz, C)
    nx, ny, nz = meta["nx"], meta["ny"], meta["nz"]
    nx_n, ny_n, nz_n = nx + 1, ny + 1, nz + 1
    n_nodes = nx_n * ny_n * nz_n

    # Mark used nodes (only those belonging to a solid voxel).
    used = np.zeros(n_nodes, dtype=bool)
    corners = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=np.int32,
    )
    elem_nodes = np.zeros((vox.shape[0], 8), dtype=np.int32)
    for e, (i, j, k) in enumerate(vox):
        nd = []
        for di, dj, dk in corners:
            nid = node_id(i + di, j + dj, k + dk, ny_n, nz_n)
            nd.append(nid)
        elem_nodes[e] = nd
        used[nd] = True

    old_to_new = -np.ones(n_nodes, dtype=np.int32)
    used_ids = np.nonzero(used)[0]
    old_to_new[used_ids] = np.arange(used_ids.size, dtype=np.int32)
    n_used = int(used_ids.size)
    elem_nodes = old_to_new[elem_nodes]

    # COO assembly, 24x24 per element.
    n_e = elem_nodes.shape[0]
    dofs = np.stack(
        [3 * elem_nodes + 0, 3 * elem_nodes + 1, 3 * elem_nodes + 2],
        axis=2,
    ).reshape(n_e, 24)
    ii = np.repeat(dofs, 24, axis=1).ravel()
    jj = np.tile(dofs, 24).ravel()
    vv = np.broadcast_to(Ke.ravel(), (n_e, 24 * 24)).ravel()
    n_dof = 3 * n_used
    K = coo_matrix((vv, (ii, jj)), shape=(n_dof, n_dof)).tocsr()
    meta.update(
        {
            "n_used_nodes": n_used,
            "n_dof": n_dof,
            "assemble_s": time.time() - t0,
        }
    )
    return {
        "K": K,
        "used_ids": used_ids,
        "old_to_new": old_to_new,
        "meta": meta,
        "nx_n": nx_n,
        "ny_n": ny_n,
        "nz_n": nz_n,
        "hx": hx,
        "hy": hy,
        "hz": hz,
    }


def nodal_xyz(used_ids: np.ndarray, nx_n: int, ny_n: int, nz_n: int, hx: float, hy: float, hz: float) -> np.ndarray:
    # Inverse of node_id(i,j,k) = (i*ny_n + j)*nz_n + k
    k = used_ids % nz_n
    tmp = used_ids // nz_n
    j = tmp % ny_n
    i = tmp // ny_n
    xyz = np.column_stack([i * hx, j * hy, k * hz])
    return xyz


def apply_bend_bcs(asm: dict, p: LatticeParams, P_half: float) -> tuple[csr_matrix, np.ndarray, dict]:
    """Simply supported at support_inset (uz=0), midspan symmetry (ux=0),
    load P_half downward on the top midspan edge. One uy=0 to kill twist.
    """
    xyz = nodal_xyz(
        asm["used_ids"], asm["nx_n"], asm["ny_n"], asm["nz_n"], asm["hx"], asm["hy"], asm["hz"]
    )
    hx, hy, hz = asm["hx"], asm["hy"], asm["hz"]
    Lh = p.length / 2.0
    tol = 0.51 * max(hx, hy, hz)

    ux_fixed = np.nonzero(np.abs(xyz[:, 0] - Lh) <= tol)[0]
    uz_support = np.nonzero(
        (np.abs(xyz[:, 0] - p.support_inset) <= tol) & (xyz[:, 2] <= 0.51 * hz)
    )[0]
    # Midspan top nodes (load land): x ≈ L/2, z ≈ H
    load_nodes = np.nonzero(
        (np.abs(xyz[:, 0] - Lh) <= tol) & (np.abs(xyz[:, 2] - p.height) <= tol)
        & (xyz[:, 1] >= -tol)
        & (xyz[:, 1] <= p.width + tol)
    )[0]
    # Pin one mid-width support node in y (and that is enough with ux, uz).
    if uz_support.size:
        mid_y = np.argmin(np.abs(xyz[uz_support, 1] - 0.5 * p.width))
        uy_pin = np.array([uz_support[mid_y]], dtype=int)
    else:
        uy_pin = np.array([], dtype=int)

    fixed = set()
    for n in ux_fixed:
        fixed.add(3 * int(n) + 0)
    for n in uz_support:
        fixed.add(3 * int(n) + 2)
    for n in uy_pin:
        fixed.add(3 * int(n) + 1)
    if not ux_fixed.size or not uz_support.size or not load_nodes.size:
        raise RuntimeError(
            f"BC node sets empty: ux={ux_fixed.size}, uz={uz_support.size}, load={load_nodes.size}"
        )

    f = np.zeros(asm["meta"]["n_dof"])
    share = P_half / load_nodes.size
    for n in load_nodes:
        f[3 * int(n) + 2] -= share  # downward

    free = np.array(sorted(set(range(asm["meta"]["n_dof"])) - fixed), dtype=int)
    info = {
        "n_ux_sym": int(ux_fixed.size),
        "n_uz_support": int(uz_support.size),
        "n_load": int(load_nodes.size),
        "n_fixed_dof": len(fixed),
        "n_free_dof": int(free.size),
        "P_half_N": P_half,
        "load_node_share_N": share,
    }
    return free, f, info


def midspan_deflection(u: np.ndarray, asm: dict, p: LatticeParams) -> float:
    """Mean downward deflection of the midspan top load nodes (mm, positive down)."""
    xyz = nodal_xyz(
        asm["used_ids"], asm["nx_n"], asm["ny_n"], asm["nz_n"], asm["hx"], asm["hy"], asm["hz"]
    )
    hx, hz = asm["hx"], asm["hz"]
    Lh = p.length / 2.0
    tol = 0.51 * max(hx, hz)
    load_nodes = np.nonzero(
        (np.abs(xyz[:, 0] - Lh) <= tol) & (np.abs(xyz[:, 2] - p.height) <= tol)
    )[0]
    uz = u[3 * load_nodes + 2]
    return float(-np.mean(uz))


def solve_coupon(p: LatticeParams, mat: dict = PLA, h: float = 0.5, P: float = 10.0) -> dict:
    """Full-coupon load P (N); half-model is loaded with P/2."""
    t0 = time.time()
    asm = assemble(p, mat, h)
    free, f, bc = apply_bend_bcs(asm, p, P_half=P / 2.0)
    K = asm["K"]
    Kff = K[free][:, free]
    ff = f[free]
    # Diagonal preconditioner.
    diag = Kff.diagonal()
    diag = np.where(np.abs(diag) > 0, diag, 1.0)
    M = csr_matrix((1.0 / diag, (np.arange(free.size), np.arange(free.size))), shape=Kff.shape)
    t_sol0 = time.time()
    uf, info = cg(Kff, ff, M=M, rtol=1e-6, maxiter=4000)
    if info != 0:
        raise RuntimeError(f"CG failed for {p.name} with info={info}")
    u = np.zeros(asm["meta"]["n_dof"])
    u[free] = uf
    delta = midspan_deflection(u, asm, p)
    if delta <= 0.0:
        raise RuntimeError(f"non-positive midspan deflection for {p.name}: {delta}")
    k = P / delta  # N/mm, full coupon
    elapsed = time.time() - t0
    result = {
        "name": p.name,
        "label": "computational estimate (linear-elastic voxel FEA, assumed isotropic PLA)",
        "voxel_h_mm": h,
        "P_applied_N": P,
        "delta_midspan_mm": delta,
        "k_fea_N_per_mm": k,
        "cg_info": int(info),
        "solve_s": time.time() - t_sol0,
        "total_s": elapsed,
        "bc": bc,
        "mesh": {k2: (v if not isinstance(v, (np.floating, np.integer)) else float(v) if isinstance(v, np.floating) else int(v)) for k2, v in asm["meta"].items()},
        "material": {"E_MPa": mat["E_MPa"], "nu": mat["nu"]},
    }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--results-dir", type=Path, default=root / "results")
    parser.add_argument("--h", type=float, default=0.5, help="voxel size (mm)")
    parser.add_argument("--P", type=float, default=10.0, help="full-coupon load (N)")
    args = parser.parse_args(argv)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for p in COUPON_VARIANTS:
        print(f"FEA {p.name} (h={args.h} mm) ...", flush=True)
        rec = solve_coupon(p, PLA, h=args.h, P=args.P)
        results.append(rec)
        print(
            f"  k_fea={rec['k_fea_N_per_mm']:.2f} N/mm, "
            f"δ({rec['P_applied_N']:.0f} N)={rec['delta_midspan_mm']:.4f} mm, "
            f"free DOF={rec['bc']['n_free_dof']}, t={rec['total_s']:.1f}s",
            flush=True,
        )

    # Merge into comparison.json if present.
    out = {
        "disclaimer": "All FEA results are computational estimates, not physical test data.",
        "voxel_h_mm": args.h,
        "results": results,
    }
    path = args.results_dir / "fea_bend.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("Wrote", path)

    cmp_path = args.results_dir / "comparison.json"
    if cmp_path.exists():
        payload = json.loads(cmp_path.read_text(encoding="utf-8"))
        by_name = {r["name"]: r for r in results}
        for c in payload["coupons"]:
            if c["name"] in by_name:
                fea = by_name[c["name"]]
                c["fea"] = {
                    "k_fea_N_per_mm": fea["k_fea_N_per_mm"],
                    "delta_midspan_mm_at_P": fea["delta_midspan_mm"],
                    "P_applied_N": fea["P_applied_N"],
                    "voxel_h_mm": fea["voxel_h_mm"],
                    "k_fea_over_k_timoshenko": fea["k_fea_N_per_mm"]
                    / c["beam_theory"]["k_timoshenko_N_per_mm"],
                    "label": fea["label"],
                }
        cmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print("Updated", cmp_path)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
