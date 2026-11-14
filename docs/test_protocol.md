# Physical test protocol and worksheet

**This is a scaled PLA concept study. It is not ballistic armor.** A print from these STLs cannot be compared to NIJ, military, or industry armor ratings, and it is not a stand-in for Oceanit’s materials work. The computational numbers in `results/` are estimates under assumed isotropic PLA constants, not lab data.

A printable one-page form is in `figures/fig_test_worksheet.png` and `figures/test_worksheet.pdf`.

## Safety

- Eye protection on for printing-finish, bend, and drop tests.
- PLA shards are sharp; deburr before you handle a coupon a lot.
- Keep fingers out from under the load. If the fixture tips or the coupon shatters, stop.
- Drop tests use a *blunt* known mass only. Do not use BBs, pellets, arrows, knives, or anything that can ricochet or penetrate.
- Work on a stable bench with a backstop (cardboard box or similar) for the drop test.

## 1. Print and weigh

Print settings: `print/bambu_x1c_coolplate.md`. Flat on the Cool Plate, basal skin down, PLA, no supports. The upper skin is a short bridge (about 3–8 mm). After printing, hold each lattice coupon up to a light and confirm the cells are hollow and the top skin actually bridged.

Predicted masses (analytical volume × 1.24 g/cm³; computational, not a promise about your filament):

| Variant | Predicted mass (g) |
|---|---:|
| coupon_solid | 18.60 |
| coupon_lattice_low | 13.78 |
| coupon_lattice_mid | 15.21 |
| coupon_lattice_high | 16.48 |
| panel_demo_mid | 45.23 |

Weigh each part (0.01 g resolution if you have it). A large miss usually means under-extrusion, a slicer that filled the cells, or a first-layer issue.

Record print time and any warp, stringing, or sag in the top skin.

## 2. Three-point bend

- Coupon envelope: 120 × 25 × 5 mm.
- Support span **S = 100 mm** (supports 10 mm in from each end), sitting on the **solid end pads**.
- Load on the midspan land (the slightly thicker-looking central strip) with a wooden dowel, metal rod, or a printed roller. Do not load a window.
- Measure midspan deflection from a fixed camera and a scale in the photo, or a dial indicator if you have one. Take the unloaded photo first.
- Add known masses (water bottles you have weighed, calibration weights, etc.). Convert grams to newtons with P (N) = m (kg) × 9.81.
- Stay in the small-deflection range first (a few tenths of a millimetre to maybe 1–2 mm) and compute k = P / δ from the linear portion.
- Then continue, if you want a failure point, until first crack. Stop when it feels unsafe.
- Record fail mode: crack at midspan land, crack in a window corner, delamination between layers, crushing under the load, or support indentation.

Compare measured k to the Timoshenko and voxel-FEA estimates in `results/comparison.csv`. As-printed PLA is anisotropic and weaker between layers, so a 20–40% difference from the isotropic estimate would not be surprising. Do not “correct” the spreadsheet to make it match.

## 3. Optional low-energy drop / impact

- Simply support the coupon on the end pads, same 100 mm span, midspan land facing up.
- Drop a known blunt mass (for example 50–200 g steel or a printed weight with a rounded face) from a measured height onto the midspan land.
- Start low (2–5 cm) and step up. One drop per height, or note if you repeat.
- Impact energy U = m g h. Example: 0.10 kg from 0.10 m is about 0.098 J.
- Record elastic (no visible damage), dent, crack, or split, and take a photo.
- This is a qualitative toughness / damage-mode test at very low energy. It is **not** an armor test.

## 4. What you can and cannot claim

You can claim: “I printed Morpho-inspired sandwich coupons, weighed them, and measured 3-point-bend stiffness on a 100 mm span.”

You cannot claim: ballistic performance, a strength-to-weight record, a result on behalf of Oceanit, or that the PLA part is “biomimetic armor.”
