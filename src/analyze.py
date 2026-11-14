#!/usr/bin/env python3
"""Computational mass, density, and flexural estimates for the lattice coupons.

Every stiffness and strength number in this module is a *computational
estimate*. No physical tests were performed. Material constants are
assumed literature-typical values for isotropic PLA, not measurements
of a printed coupon.

Primary model: Euler-Bernoulli (and Timoshenko-corrected) beam theory
with the actual ribbed cross-section I(x) along the span.

Comparison model: Gibson-Ashby cellular-solid scaling applied to the
same envelope, which is the right order-of-magnitude tool for foams
but is expected to under-predict stiffness for this extruded, nearly
full-height ribbed plate (the ridges carry bending stress like spars).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from geometry import (
    COUPON_VARIANTS,
    LatticeParams,
    analytical_volume_mm3,
    lattice_region_relative_density,
    section_properties,
)

# ---------------------------------------------------------------------------
# Assumed PLA properties (computational inputs, not measured).
# Density ~1.24 g/cm^3 is a typical filament datasheet value.
# As-printed PLA modulus is anisotropic and usually lower than injection-
# molded PLA; published tensile/flexural moduli commonly fall around
# 2–3.5 GPa depending on orientation, infill, and temperature. 2.7 GPa
# is a mid-range assumption used uniformly for every variant.
# Assumed flexural stress at first crack (60 MPa) is likewise a mid-range
# literature-typical value for printed PLA, used only to scale a
# first-crack *load estimate*. It is not a test result.
# ---------------------------------------------------------------------------
PLA = {
    "name": "PLA (assumed isotropic, literature-typical)",
    "density_g_cm3": 1.24,
    "E_MPa": 2700.0,
    "nu": 0.33,
    "flexural_stress_MPa": 60.0,
    "notes": (
        "Assumed computational inputs, not measured on these coupons. "
        "As-printed PLA is anisotropic; interlayer bonding can be weaker "
        "than the in-plane modulus used here."
    ),
}

# Gibson & Ashby, Cellular Solids, 2nd ed. (1997), Cambridge University Press.
# Open-cell, bending-dominated foams (Eq. 5.6 / Table 5.1 style):
#   E*/Es ≈ C1 (ρ*/ρs)^2,  C1 ~ 1
#   σ*/σs ≈ C2 (ρ*/ρs)^{3/2}, C2 ~ 0.3
# Stretch-dominated lattices (upper-bound / honeycomb out-of-plane style):
#   E*/Es ≈ C3 (ρ*/ρs), C3 taken as 1.0 as a simple rule of thumb
#   σ*/σs ≈ C4 (ρ*/ρs), C4 taken as 0.3
GIBSON_ASHBY = {
    "C1_open_cell_modulus": 1.0,
    "C2_open_cell_strength": 0.3,
    "C3_stretch_modulus": 1.0,
    "C4_stretch_strength": 0.3,
    "source": "Gibson, L. J., & Ashby, M. F. (1997). Cellular Solids: Structure and Properties (2nd ed.). Cambridge University Press. Constants are the conventional order-of-magnitude values, not fitted to these prints.",
}


def coupon_mass_g(p: LatticeParams, density_g_cm3: float) -> dict[str, float]:
    vol = analytical_volume_mm3(p)
    rel_env = vol / p.envelope_volume_mm3
    rel_lat = lattice_region_relative_density(p)
    mass = vol * density_g_cm3 / 1000.0  # mm^3 * g/cm^3 / 1000 = g
    return {
        "volume_mm3": vol,
        "envelope_volume_mm3": p.envelope_volume_mm3,
        "relative_density_envelope": rel_env,
        "relative_density_lattice_region": rel_lat,
        "mass_g": mass,
    }


def beam_theory(p: LatticeParams, mat: dict, n: int = 2400) -> dict:
    """3-point-bend stiffness and first-crack load from I(x).

    Span S is between supports at `support_inset` from each end.
    Moment M(s) = (P/2) s on the left half (standard simply-supported
    center-loaded beam). Units: mm, N, MPa (N/mm^2).
    """
    E = mat["E_MPa"]
    nu = mat["nu"]
    G = E / (2.0 * (1.0 + nu))
    sigma_f = mat["flexural_stress_MPa"]
    S = p.span_mm
    x0 = p.support_inset
    s = np.linspace(0.0, S, n)
    x = x0 + s
    props = [section_properties(xi, p) for xi in x]
    I = np.array([q["I_mm4"] for q in props])
    A = np.array([q["area_mm2"] for q in props])
    ymax = np.array([q["ymax_mm"] for q in props])
    Mp = np.where(s <= S / 2.0, 0.5 * s, 0.5 * (S - s))  # M/P, mm
    Vp = np.full_like(s, 0.5)
    Vp[s > S / 2.0] = -0.5

    # Castigliano / virtual work: δ/P = ∫ (M/P)^2 / (E I) ds
    flex_comp = np.trapezoid((Mp**2) / I, s) / E  # mm/N
    k_eb = 1.0 / flex_comp  # N/mm

    # Timoshenko shear, form factor 6/5 for a rectangular section
    # (used here as an approximate correction even for the ribbed section).
    f_s = 6.0 / 5.0
    shear_comp = f_s * np.trapezoid((Vp**2) / A, s) / G  # mm/N
    k_tim = 1.0 / (flex_comp + shear_comp)

    # Reference: constant-I formula with the midspan-adjacent lattice I
    # and with the solid-section I.
    I_solid = p.width * p.height**3 / 12.0
    # Section just outside the load pad (weakest high-moment station).
    s_edge = max(0.0, S / 2.0 - p.load_pad / 2.0 - 0.05)
    I_edge = section_properties(x0 + s_edge, p)["I_mm4"]
    k_const_edge = 48.0 * E * I_edge / S**3
    k_const_solid = 48.0 * E * I_solid / S**3

    # Effective modulus vs the solid envelope I, back-calculated from k_eb.
    E_eff = k_eb * S**3 / (48.0 * I_solid)

    # First-crack load estimate: max of |M| ymax / I along the span.
    # Avoid the exact midspan solid pad when the lattice section next to
    # it is more critical.
    sigma_per_P = np.abs(Mp) * ymax / I  # 1/mm^2 = 1/N * MPa? M/P is mm, ymax mm, I mm^4
    # σ (MPa) = (M_Nmm) * ymax / I = (P_N * Mp_mm) * ymax / I  → σ/P = Mp*ymax/I  (1/mm^2)
    # P (N) = sigma_f (N/mm^2) / (σ/P)
    crit = int(np.argmax(sigma_per_P))
    P_crack = sigma_f / sigma_per_P[crit]  # N

    return {
        "span_mm": S,
        "support_inset_mm": p.support_inset,
        "E_assumed_MPa": E,
        "G_assumed_MPa": G,
        "nu_assumed": nu,
        "flexural_stress_assumed_MPa": sigma_f,
        "k_euler_bernoulli_N_per_mm": k_eb,
        "k_timoshenko_N_per_mm": k_tim,
        "shear_share_of_deflection": shear_comp / (flex_comp + shear_comp),
        "k_constI_edge_N_per_mm": k_const_edge,
        "k_constI_solid_section_N_per_mm": k_const_solid,
        "I_solid_mm4": I_solid,
        "I_min_in_span_mm4": float(np.min(I)),
        "I_at_load_pad_edge_mm4": I_edge,
        "E_eff_from_k_and_Isolid_MPa": E_eff,
        "E_eff_over_Es": E_eff / E,
        "P_first_crack_est_N": P_crack,
        "critical_s_mm": float(s[crit]),
        "critical_x_mm": float(x[crit]),
        "sigma_per_P_max_per_mm2": float(sigma_per_P[crit]),
        "label": "computational estimate (beam theory, assumed isotropic PLA)",
    }


def gibson_ashby_estimates(p: LatticeParams, mat: dict, mass: dict) -> dict:
    """Naive foam-style scaling on the envelope, for comparison only."""
    Es = mat["E_MPa"]
    sigma_s = mat["flexural_stress_MPa"]
    # Use lattice-region relative density for the cellular-solid scaling;
    # the solid control is φ = 1.
    phi = mass["relative_density_lattice_region"] if p.kind != "solid" else 1.0
    C = GIBSON_ASHBY
    E_open = Es * C["C1_open_cell_modulus"] * phi**2
    E_stretch = Es * C["C3_stretch_modulus"] * phi
    sig_open = sigma_s * C["C2_open_cell_strength"] * phi**1.5
    sig_stretch = sigma_s * C["C4_stretch_strength"] * phi
    I_solid = p.width * p.height**3 / 12.0
    S = p.span_mm
    k_open = 48.0 * E_open * I_solid / S**3
    k_stretch = 48.0 * E_stretch * I_solid / S**3
    # First-crack load from σ = M c / I, c = H/2, M = P S / 4
    # P = σ I / (c S/4) = 4 σ I / (c S)
    c = p.height / 2.0
    P_open = 4.0 * sig_open * I_solid / (c * S)
    P_stretch = 4.0 * sig_stretch * I_solid / (c * S)
    return {
        "phi_used": phi,
        "E_open_cell_MPa": E_open,
        "E_stretch_MPa": E_stretch,
        "k_open_cell_N_per_mm": k_open,
        "k_stretch_N_per_mm": k_stretch,
        "sigma_open_cell_MPa": sig_open,
        "sigma_stretch_MPa": sig_stretch,
        "P_open_cell_N": P_open,
        "P_stretch_N": P_stretch,
        "constants": C,
        "label": (
            "computational estimate (Gibson-Ashby envelope scaling). "
            "This treats the coupon as a homogeneous foam of relative "
            "density φ and is shown as a comparison, not as the primary "
            "model for these extruded ribs."
        ),
    }


def analyze_coupon(p: LatticeParams, mat: dict = PLA) -> dict:
    mass = coupon_mass_g(p, mat["density_g_cm3"])
    beam = beam_theory(p, mat)
    ga = gibson_ashby_estimates(p, mat, mass)
    k = beam["k_timoshenko_N_per_mm"]
    specific_k = k / mass["mass_g"]  # N/mm per gram
    specific_P = beam["P_first_crack_est_N"] / mass["mass_g"]
    specific_k_solid_ref = None
    return {
        "name": p.name,
        "kind": p.kind,
        "material_assumptions": mat,
        "mass": mass,
        "beam_theory": beam,
        "gibson_ashby": ga,
        "specific_stiffness_N_per_mm_per_g": specific_k,
        "specific_first_crack_load_N_per_g": specific_P,
        "disclaimer": (
            "All stiffness, strength, and specific-property values are "
            "computational estimates under the assumed isotropic PLA "
            "constants. They are not physical test data."
        ),
    }


def add_ratios(rows: list[dict]) -> None:
    solid = next(r for r in rows if r["kind"] == "solid")
    k_s = solid["beam_theory"]["k_timoshenko_N_per_mm"]
    m_s = solid["mass"]["mass_g"]
    P_s = solid["beam_theory"]["P_first_crack_est_N"]
    sk_s = solid["specific_stiffness_N_per_mm_per_g"]
    sp_s = solid["specific_first_crack_load_N_per_g"]
    for r in rows:
        r["ratio_vs_solid"] = {
            "mass": r["mass"]["mass_g"] / m_s,
            "k_timoshenko": r["beam_theory"]["k_timoshenko_N_per_mm"] / k_s,
            "P_first_crack_est": r["beam_theory"]["P_first_crack_est_N"] / P_s,
            "specific_stiffness": r["specific_stiffness_N_per_mm_per_g"] / sk_s,
            "specific_first_crack_load": r["specific_first_crack_load_N_per_g"] / sp_s,
        }


def flatten_row(r: dict) -> dict:
    m, b, g, q = r["mass"], r["beam_theory"], r["gibson_ashby"], r["ratio_vs_solid"]
    return {
        "name": r["name"],
        "kind": r["kind"],
        "volume_mm3": m["volume_mm3"],
        "mass_g": m["mass_g"],
        "relative_density_envelope": m["relative_density_envelope"],
        "relative_density_lattice_region": m["relative_density_lattice_region"],
        "k_euler_bernoulli_N_per_mm": b["k_euler_bernoulli_N_per_mm"],
        "k_timoshenko_N_per_mm": b["k_timoshenko_N_per_mm"],
        "E_eff_MPa": b["E_eff_from_k_and_Isolid_MPa"],
        "E_eff_over_Es": b["E_eff_over_Es"],
        "P_first_crack_est_N": b["P_first_crack_est_N"],
        "specific_stiffness_N_per_mm_per_g": r["specific_stiffness_N_per_mm_per_g"],
        "specific_first_crack_load_N_per_g": r["specific_first_crack_load_N_per_g"],
        "k_GA_open_cell_N_per_mm": g["k_open_cell_N_per_mm"],
        "k_GA_stretch_N_per_mm": g["k_stretch_N_per_mm"],
        "mass_vs_solid": q["mass"],
        "k_vs_solid": q["k_timoshenko"],
        "specific_stiffness_vs_solid": q["specific_stiffness"],
        "P_vs_solid": q["P_first_crack_est"],
        "specific_P_vs_solid": q["specific_first_crack_load"],
        "label": "computational estimate; not physical test data",
    }


def write_markdown_table(rows: list[dict], path: Path) -> None:
    flat = [flatten_row(r) for r in rows]
    cols = [
        ("name", "Variant", "{:s}"),
        ("mass_g", "Mass (g)", "{:.2f}"),
        ("relative_density_envelope", "ρ*/ρs (envelope)", "{:.3f}"),
        ("relative_density_lattice_region", "ρ*/ρs (lattice region)", "{:.3f}"),
        ("k_timoshenko_N_per_mm", "k, Timoshenko (N/mm)", "{:.2f}"),
        ("E_eff_over_Es", "E_eff / E_s", "{:.3f}"),
        ("P_first_crack_est_N", "P first-crack est. (N)", "{:.1f}"),
        ("specific_stiffness_N_per_mm_per_g", "k / mass (N/mm/g)", "{:.2f}"),
        ("specific_first_crack_load_N_per_g", "P / mass (N/g)", "{:.1f}"),
        ("specific_stiffness_vs_solid", "specific k vs solid", "{:.2f}"),
    ]
    headers = [c[1] for c in cols]
    lines = [
        "# Computational comparison (not physical test data)",
        "",
        "Assumed isotropic PLA: E = 2.7 GPa, ν = 0.33, density = 1.24 g/cm³, "
        "flexural stress at first crack = 60 MPa. Three-point bend span = 100 mm. "
        "Primary stiffness is Timoshenko beam theory with the ribbed I(x). "
        "Gibson-Ashby foam scaling is included in the JSON/CSV as a comparison only.",
        "",
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in flat:
        vals = []
        for key, _, fmt in cols:
            vals.append(fmt.format(row[key]))
        lines.append("| " + " | ".join(vals) + " |")
    lines += [
        "",
        "These numbers will change if the print orientation, infill, or actual "
        "PLA modulus differs from the assumptions. Weigh each print and run the "
        "worksheet in `docs/test_protocol.md` before claiming an experimental result.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--results-dir", type=Path, default=root / "results")
    args = parser.parse_args(argv)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    rows = [analyze_coupon(p) for p in COUPON_VARIANTS]
    add_ratios(rows)
    payload = {
        "disclaimer": "All results are computational estimates, not physical test data.",
        "material_assumptions": PLA,
        "gibson_ashby_constants": GIBSON_ASHBY,
        "coupons": rows,
    }
    json_path = args.results_dir / "comparison.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    flat = [flatten_row(r) for r in rows]
    csv_path = args.results_dir / "comparison.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0].keys()))
        w.writeheader()
        w.writerows(flat)

    write_markdown_table(rows, args.results_dir / "comparison.md")
    print("Wrote", json_path)
    print("Wrote", csv_path)
    for r in flat:
        print(
            f"  {r['name']}: mass={r['mass_g']:.2f} g, "
            f"φ_env={r['relative_density_envelope']:.3f}, "
            f"k={r['k_timoshenko_N_per_mm']:.2f} N/mm, "
            f"k/m={r['specific_stiffness_N_per_mm_per_g']:.2f} N/mm/g, "
            f"P_est={r['P_first_crack_est_N']:.1f} N"
        )
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
