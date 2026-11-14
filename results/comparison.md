# Computational comparison (not physical test data)

Assumed isotropic PLA: **E = 2.7 GPa**, ν = 0.33, density = 1.24 g/cm³, flexural stress at first crack = 60 MPa. These are literature-typical *inputs*, not measurements of a printed coupon.

Three-point bend span **S = 100 mm** (supports 10 mm from each end). Primary stiffness is Timoshenko beam theory with the ribbed/sandwich I(x). A linear-elastic voxel FEA (0.5 mm hex elements, half-span model) is shown as a check; it agrees with beam theory to about 2%. Gibson–Ashby foam scaling is in `comparison.json` / `comparison.csv` as a comparison only — these extruded two-skin panels sit *above* open-cell (ρ*)² scaling, as expected for a sandwich with continuous ridges.

| Variant | Mass (g) | ρ*/ρs envelope | ρ*/ρs lattice region | k Timoshenko (N/mm) | k voxel FEA (N/mm) | E_eff / E_s | P first-crack est. (N) | k / mass (N/mm/g) | specific k vs solid |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| solid | 18.60 | 1.000 | 1.000 | 33.48 | 33.87 | 1.000 | 250.1 | 1.80 | 1.00 |
| lattice_low | 13.78 | 0.741 | 0.621 | 29.63 | 29.23 | 0.888 | 229.6 | 2.15 | 1.19 |
| lattice_mid | 15.21 | 0.818 | 0.733 | 30.83 | 30.78 | 0.922 | 237.6 | 2.03 | 1.13 |
| lattice_high | 16.48 | 0.886 | 0.833 | 31.81 | 31.23 | 0.951 | 245.1 | 1.93 | 1.07 |

The low-density coupon is the lightest (about 74% of the solid mass) and has the highest *estimated* specific flexural stiffness (about 1.19× the solid plate) and specific first-crack load (about 1.24×). That is a sandwich-theory result under the assumed isotropic PLA constants, **not** a physical test and **not** an armor rating.

These numbers will change if the print orientation, interlayer bonding, or actual PLA modulus differs from the assumptions. Weigh each print and run `docs/test_protocol.md` before claiming an experimental result.
