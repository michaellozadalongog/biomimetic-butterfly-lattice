"""Parametric Morpho-inspired ridge / cross-rib / lower-lamina geometry.

This is a *scaled, printable abstraction* of published butterfly-scale
architecture, not a nanoscale replica. Natural ridge spacing is on the
order of 1–2 µm; these coupons are on the order of 100 mm.

Biological sources (public):
- Ghiradella, H. (1998). Hairs, bristles, and scales. In F. W. Harrison
  & M. Locke (Eds.), Microscopic Anatomy of Invertebrates, Vol. 11A:
  Insecta (pp. 257–287). Wiley-Liss.
- Giraldo, M. A., & Stavenga, D. G. (2016). Brilliant iridescence of
  Morpho butterfly wing scales is due to both a thin film lower lamina
  and a multilayered upper lamina. Journal of Comparative Physiology A,
  202, 381–388. https://doi.org/10.1007/s00359-016-1084-1
- Giraldo, M. A., & Stavenga, D. G. (2016). Coloration mechanisms and
  phylogeny of Morpho butterflies. Journal of Experimental Biology,
  219, 3936–3944. https://doi.org/10.1242/jeb.148726
- Ding, Y., Xu, S., & Wang, Z. L. (2009). Structural colors from
  Morpho peleides butterfly wing scales. Journal of Applied Physics,
  106, 074702. https://doi.org/10.1063/1.3239513
  (ridge separation ~1.7–1.8 µm; trabeculae connect shelves to a
  ~250 nm base plate; cross-ribs span neighboring ridges)
- Davis, A. L., et al. (2020). Diverse nanostructures underlie thin
  ultra-black scales in butterflies. Nature Communications, 11, 1294.
  https://doi.org/10.1038/s41467-020-15033-1
  (papilionid scales often show a honeycomb-like upper lamina of
  ridges + cross-ribs; trabeculae connect the laminae)

Mapped features (printable, all walls vertical except a short-span
top skin that the slicer can bridge):
- lower lamina  -> continuous basal skin of thickness t_skin
- upper lamina  -> continuous top skin of thickness t_skin (bridged
                   across the windows; typical spans are 3–8 mm)
- longitudinal ridges -> full-height ribs along the coupon length
- cross-ribs    -> full-height transverse ribs
- trabeculae    -> abstracted as the vertical walls that join the
                   two skins (not separate nanoscale pillars)
- lamellae / "Christmas-tree" ridge shelves are optical, ~10^1–10^2 nm,
  and would be unsupported overhangs if scaled up literally, so they
  are *not* printed. They are discussed only as biological context.

All dimensions are millimetres.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Iterable

import numpy as np

# Bambu Lab X1 Carbon usable area on the Cool Plate (mm). Official
# build volume is 256 x 256 x 256 mm; stay inside a small margin.
BED_SIZE_MM = 256.0
BED_MARGIN_MM = 4.0


@dataclass(frozen=True)
class LatticeParams:
    name: str
    length: float = 120.0
    width: float = 25.0
    height: float = 5.0
    t_skin: float = 1.0
    frame_width: float = 2.0
    end_pad: float = 15.0
    load_pad: float = 8.0
    ridge_pitch: float = 5.0
    ridge_width: float = 1.4
    cross_pitch: float = 6.0
    cross_width: float = 1.2
    # Recommended 3-point-bend support inset from each end (mm).
    support_inset: float = 10.0
    kind: str = "lattice"  # "lattice" or "solid"

    def __post_init__(self) -> None:
        if 2.0 * self.t_skin >= self.height:
            raise ValueError("two skins would fill the thickness; t_skin must be thinner")
        if self.end_pad * 2 + self.load_pad >= self.length:
            raise ValueError("pads leave no lattice span")
        if min(self.ridge_width, self.cross_width, self.frame_width, self.t_skin) < 0.8:
            raise ValueError("features below 0.8 mm are risky on a 0.4 mm nozzle")
        footprint = max(self.length, self.width)
        if footprint > BED_SIZE_MM - 2 * BED_MARGIN_MM:
            raise ValueError("footprint exceeds the X1C Cool Plate printable area")

    @property
    def envelope_volume_mm3(self) -> float:
        return self.length * self.width * self.height

    @property
    def span_mm(self) -> float:
        return self.length - 2.0 * self.support_inset

    def to_dict(self) -> dict:
        d = asdict(self)
        d["envelope_volume_mm3"] = self.envelope_volume_mm3
        d["span_mm"] = self.span_mm
        d["ridge_centers_mm"] = list(self.ridge_centers())
        d["cross_centers_mm"] = list(self.cross_centers())
        return d

    def ridge_centers(self) -> np.ndarray:
        if self.kind == "solid":
            return np.array([])
        # Overlap the side frames slightly so the boolean union cannot
        # leave a four-fold knife edge where a ridge meets a frame.
        overlap = 0.2
        inner0 = max(0.0, self.frame_width - overlap)
        inner1 = min(self.width, self.width - self.frame_width + overlap)
        return _centered_pitch(inner0, inner1, self.ridge_width, self.ridge_pitch)

    def cross_centers(self) -> np.ndarray:
        if self.kind == "solid":
            return np.array([])
        overlap = 0.2
        inner0 = max(0.0, self.end_pad - overlap)
        inner1 = min(self.length, self.length - self.end_pad + overlap)
        # Keep cross-ribs out of the midspan load pad so the load land stays a
        # clean rectangle (union would still work; this keeps the cell even).
        centers = _centered_pitch(inner0, inner1, self.cross_width, self.cross_pitch)
        mid0 = self.length / 2.0 - self.load_pad / 2.0
        mid1 = self.length / 2.0 + self.load_pad / 2.0
        keep = []
        for c in centers:
            a, b = c - self.cross_width / 2.0, c + self.cross_width / 2.0
            if b <= mid0 or a >= mid1:
                keep.append(c)
        return np.array(keep, dtype=float)


def _centered_pitch(start: float, end: float, width: float, pitch: float) -> np.ndarray:
    """Place as many `width` ribs as fit at approximately `pitch`, centered."""
    span = end - start
    if span < width - 1e-9:
        return np.array([])
    n = int(np.floor((span - width) / pitch + 1e-9)) + 1
    total = width + (n - 1) * pitch
    origin = start + 0.5 * (span - total)
    return origin + width / 2.0 + pitch * np.arange(n)


# Coupon variants: same envelope, three relative densities, plus solid control.
# Widths stay >= 1.0 mm for reliable Cool Plate printing with a 0.4 mm nozzle.
COUPON_SOLID = LatticeParams(name="coupon_solid", kind="solid", ridge_width=1.0, cross_width=1.0)
COUPON_LOW = LatticeParams(
    name="coupon_lattice_low",
    ridge_width=1.0,
    cross_width=1.0,
    ridge_pitch=7.0,
    cross_pitch=8.5,
)
COUPON_MID = LatticeParams(
    name="coupon_lattice_mid",
    ridge_width=1.6,
    cross_width=1.4,
    ridge_pitch=5.0,
    cross_pitch=6.0,
)
COUPON_HIGH = LatticeParams(
    name="coupon_lattice_high",
    ridge_width=2.4,
    cross_width=2.0,
    ridge_pitch=5.0,
    cross_pitch=5.5,
)

# Square demo tile (mid density) for photos; not the bend coupon.
DEMO_MID = LatticeParams(
    name="panel_demo_mid",
    length=100.0,
    width=100.0,
    height=5.0,
    t_skin=1.0,
    frame_width=2.5,
    end_pad=2.5,  # same as frame; a uniform border, no bend pads
    load_pad=0.0,
    ridge_pitch=5.0,
    ridge_width=1.6,
    cross_pitch=5.0,
    cross_width=1.4,
    support_inset=0.0,
    kind="lattice",
)

COUPON_VARIANTS: tuple[LatticeParams, ...] = (
    COUPON_SOLID,
    COUPON_LOW,
    COUPON_MID,
    COUPON_HIGH,
)
ALL_VARIANTS: tuple[LatticeParams, ...] = COUPON_VARIANTS + (DEMO_MID,)


def is_full_height(x: np.ndarray, y: np.ndarray, p: LatticeParams) -> np.ndarray:
    """True where the part is solid through the full coupon thickness."""
    x = np.asarray(x)
    y = np.asarray(y)
    inside = (x >= 0.0) & (x <= p.length) & (y >= 0.0) & (y <= p.width)
    if p.kind == "solid":
        return inside
    on = np.zeros(np.broadcast(x, y).shape, dtype=bool)
    on |= (x <= p.end_pad) | (x >= p.length - p.end_pad)
    if p.load_pad > 0.0:
        mid0 = p.length / 2.0 - p.load_pad / 2.0
        mid1 = p.length / 2.0 + p.load_pad / 2.0
        on |= (x >= mid0) & (x <= mid1)
    on |= (y <= p.frame_width) | (y >= p.width - p.frame_width)
    for c in p.ridge_centers():
        on |= np.abs(y - c) <= p.ridge_width / 2.0
    for c in p.cross_centers():
        on |= np.abs(x - c) <= p.cross_width / 2.0
    return on & inside


def is_solid(x: np.ndarray, y: np.ndarray, z: np.ndarray, p: LatticeParams) -> np.ndarray:
    """Point-in-solid test used by the voxel FEA and volume checks."""
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    inside = (x >= 0.0) & (x <= p.length) & (y >= 0.0) & (y <= p.width) & (z >= 0.0) & (z <= p.height)
    skin = (z <= p.t_skin) | (z >= p.height - p.t_skin)
    return inside & (skin | is_full_height(x, y, p))


def analytical_volume_mm3(p: LatticeParams, dx: float = 0.05, dy: float = 0.05) -> float:
    """Occupancy integral of the 2.5-D solid. Returns mm^3."""
    if p.kind == "solid":
        return p.envelope_volume_mm3
    nx = int(np.round(p.length / dx))
    ny = int(np.round(p.width / dy))
    xs = (np.arange(nx) + 0.5) * (p.length / nx)
    ys = (np.arange(ny) + 0.5) * (p.width / ny)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    full = is_full_height(X, Y, p)
    cell = (p.length / nx) * (p.width / ny)
    area_full = float(full.sum()) * cell
    area_skin_only = p.length * p.width - area_full
    # Two skins (upper + lower lamina) in the open cells.
    return area_full * p.height + area_skin_only * (2.0 * p.t_skin)


def section_rectangles(x: float, p: LatticeParams, dy: float = 0.05) -> list[tuple[float, float, float]]:
    """Cross-section at station x as rectangles (y_width, z0, z1) for I(x).

    The section is uniform along y in strips of width `dy`.
    """
    if x < 0.0 or x > p.length:
        return []
    if p.kind == "solid" or x <= p.end_pad or x >= p.length - p.end_pad:
        return [(p.width, 0.0, p.height)]
    if p.load_pad > 0.0:
        mid0 = p.length / 2.0 - p.load_pad / 2.0
        mid1 = p.length / 2.0 + p.load_pad / 2.0
        if mid0 <= x <= mid1:
            return [(p.width, 0.0, p.height)]
    ny = int(np.round(p.width / dy))
    ys = (np.arange(ny) + 0.5) * (p.width / ny)
    full = is_full_height(np.full_like(ys, x), ys, p)
    rects: list[tuple[float, float, float]] = []
    # Merge adjacent strips with the same height.
    i = 0
    width_step = p.width / ny
    while i < ny:
        if full[i]:
            # One full-height rectangle.
            kind = "full"
        else:
            kind = "skins"
        j = i + 1
        while j < ny:
            kind_j = "full" if full[j] else "skins"
            if kind_j != kind:
                break
            j += 1
        w = (j - i) * width_step
        if kind == "full":
            rects.append((w, 0.0, p.height))
        else:
            rects.append((w, 0.0, p.t_skin))
            rects.append((w, p.height - p.t_skin, p.height))
        i = j
    return rects


def section_properties(x: float, p: LatticeParams, dy: float = 0.05) -> dict[str, float]:
    """Area, neutral-axis height, and second moment about the NA (mm units)."""
    rects = section_rectangles(x, p, dy=dy)
    if not rects:
        return {"area_mm2": 0.0, "na_mm": 0.0, "I_mm4": 0.0, "ymax_mm": 0.0}
    areas = []
    cys = []
    i_own = []
    for w, z0, z1 in rects:
        h = z1 - z0
        a = w * h
        areas.append(a)
        cys.append(0.5 * (z0 + z1))
        i_own.append(w * h**3 / 12.0)
    A = float(sum(areas))
    na = float(sum(a * c for a, c in zip(areas, cys)) / A)
    I = float(sum(io + a * (c - na) ** 2 for io, a, c in zip(i_own, areas, cys)))
    z_top = max(z1 for _, _, z1 in rects)
    z_bot = min(z0 for _, z0, _ in rects)
    ymax = max(abs(z_top - na), abs(na - z_bot))
    return {"area_mm2": A, "na_mm": na, "I_mm4": I, "ymax_mm": ymax}


def lattice_region_relative_density(p: LatticeParams, dx: float = 0.05, dy: float = 0.05) -> float:
    """Relative density in the open lattice span, excluding solid end/load pads.

    For the solid control and the demo tile (no load pad, tiny end border)
    this falls back to the whole-envelope relative density.
    """
    if p.kind == "solid":
        return 1.0
    x0 = p.end_pad
    x1 = p.length - p.end_pad
    if p.load_pad > 0.0:
        # Exclude the midspan pad from the "open lattice" measure by
        # integrating only the open cells on either side of it.
        mid0 = p.length / 2.0 - p.load_pad / 2.0
        mid1 = p.length / 2.0 + p.load_pad / 2.0
        ranges = [(x0, mid0), (mid1, x1)]
    else:
        ranges = [(x0, x1)]
    vol = 0.0
    env = 0.0
    for a, b in ranges:
        if b <= a:
            continue
        nx = max(1, int(np.round((b - a) / dx)))
        ny = int(np.round(p.width / dy))
        xs = a + (np.arange(nx) + 0.5) * ((b - a) / nx)
        ys = (np.arange(ny) + 0.5) * (p.width / ny)
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        full = is_full_height(X, Y, p)
        cell = ((b - a) / nx) * (p.width / ny)
        area_full = float(full.sum()) * cell
        area_all = (b - a) * p.width
        vol += area_full * p.height + (area_all - area_full) * (2.0 * p.t_skin)
        env += area_all * p.height
    if env <= 0.0:
        return analytical_volume_mm3(p) / p.envelope_volume_mm3
    return vol / env
