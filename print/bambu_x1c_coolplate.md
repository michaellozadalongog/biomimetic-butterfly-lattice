# Bambu Lab X1 Carbon + Cool Plate print notes

These parts were designed to print **flat, basal (lower) skin down**, on the **Cool Plate**, in **PLA**, with **no supports**.

The official X1C build volume is 256 × 256 × 256 mm. Every STL in `stl/` fits inside a 4 mm margin of that plate.

## Recommended starting settings

Use the **Bambu Studio PLA + Cool Plate** preset as the baseline, then:

| Setting | Suggestion |
|---|---|
| Plate | Cool Plate (not the high-temperature plate) |
| Nozzle | 0.4 mm |
| Layer height | 0.20 mm |
| Walls / perimeters | 3 or more |
| Infill | 100% (the CAD *is* the structure; do not add a sparse infill on top of the windows) |
| Supports | Off |
| Brim | Optional 3–5 mm if the first layer looks shy |
| Orientation | Largest face down (the 1.0 mm lower skin on the plate) |

Bed temperature, nozzle temperature, cooling, and any glue-stick recommendation should follow the Cool Plate + PLA preset for the filament you actually use. Those values vary by filament and by Bambu Studio version, so they are not hard-coded here.

## Bridging the upper skin

The mid- and high-density windows span about 3–5 mm. The low-density windows span about 6–8 mm. The 1.0 mm upper skin (five 0.20 mm layers) is meant to be **bridged** by the slicer across those openings.

- Keep cooling on for the top skin.
- Do **not** enable supports inside the cells — they would be trapped.
- After the print, hold the coupon up to a light and check that the windows are hollow and the top skin actually bridged. If the low-density top skin sags, slow the bridge / top solid layers or add one more wall; do not invent a support scheme that fills the cells.

## Suggested plate layout

See `figures/fig_print_layout.png`. All four 120 × 25 mm coupons plus the 100 × 100 mm demo tile fit on one Cool Plate with gaps.

## What this is not

A PLA print is **not** ballistic armor. These parts are a scaled classroom / portfolio concept study of a butterfly-scale architecture. They will crack, they are temperature-sensitive, and they have nothing to do with NIJ or military armor ratings.
