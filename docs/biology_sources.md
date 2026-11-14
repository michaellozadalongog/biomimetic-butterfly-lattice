# Biological sources (public)

This project is a *scaled, printable abstraction* of published butterfly-scale architecture. It is not a replica of chitin, it is not at biological length scales, and it is not a claim about how a living wing behaves under impact.

## Architecture that was mapped

Typical nymphalid / Morpho-type scales (and many papilionid scales) have:

- a relatively smooth **lower lamina** (basal plate)
- an **upper lamina** of longitudinal **ridges** joined by **cross-ribs**, leaving **windows**
- pillar-like **trabeculae** connecting the two laminae
- in Morpho-type ridges, stacked **lamellae** (the “Christmas-tree” multilayer) that are primarily an *optical* device

Those features are described in:

- Ghiradella, H. (1998). Hairs, bristles, and scales. In F. W. Harrison & M. Locke (Eds.), *Microscopic Anatomy of Invertebrates, Vol. 11A: Insecta* (pp. 257–287). Wiley-Liss.
- Giraldo, M. A., & Stavenga, D. G. (2016). Brilliant iridescence of *Morpho* butterfly wing scales is due to both a thin film lower lamina and a multilayered upper lamina. *Journal of Comparative Physiology A, 202*, 381–388. https://doi.org/10.1007/s00359-016-1084-1
- Giraldo, M. A., & Stavenga, D. G. (2016). Coloration mechanisms and phylogeny of *Morpho* butterflies. *Journal of Experimental Biology, 219*, 3936–3944. https://doi.org/10.1242/jeb.148726 (open PDF: https://journals.biologists.com/jeb/article-pdf/219/24/3936/1893519/jeb148726.pdf)
- Ding, Y., Xu, S., & Wang, Z. L. (2009). Structural colors from *Morpho peleides* butterfly wing scales. *Journal of Applied Physics, 106*, 074702. https://doi.org/10.1063/1.3239513 — ridge separation about 1.7–1.8 µm; trabeculae connect the shelves to a ~250 nm base plate; cross-rib density differs between cover and ground scales. Author-posted PDF: https://xugroup.eng.ucsd.edu/wp-content/uploads/2015/09/09_j_appl_phys.pdf
- Davis, A. L., Nijhout, H. F., & Johnsen, S. (2020). Diverse nanostructures underlie thin ultra-black scales in butterflies. *Nature Communications, 11*, 1294. https://doi.org/10.1038/s41467-020-15033-1 — papilionid scales often show a honeycomb-like upper lamina of ridges and cross-ribs; robust trabeculae connect the laminae.

Some papilionids (for example *Parides sesostris*) also have internal **gyroid** photonic crystals. Those are optical, nanoscale, and were *not* used as the printable architecture here (a gyroid at this coupon size is harder to print support-free and is a different biological structure than the Morpho ridge / cross-rib scale).

## What was deliberately not printed

Ridge lamellae in Morpho-type scales are on the order of 10¹–10² nm (Ding, Xu & Wang report cuticle / air layers near 90 nm and 70 nm). Scaled up literally they would be unsupported shelves. They are discussed only as biological context.

Natural ridge pitch (~1.7–1.8 µm in *M. peleides*) is scaled by roughly 2,000–4,000× in these coupons (5–7 mm pitch). Nothing in this repo should be read as a statement that the printed PLA part recovers the mechanics of chitin at biological scale.

## Cellular-solid reference

- Gibson, L. J., & Ashby, M. F. (1997). *Cellular Solids: Structure and Properties* (2nd ed.). Cambridge University Press. Open-cell bending-dominated foams are conventionally described by E*/Es ≈ (ρ*/ρs)² and σ*/σs ≈ 0.3 (ρ*/ρs)^{3/2}. Those relations are used here only as a comparison curve. The printed architecture is a two-skin sandwich with continuous ridges, so section-based beam theory (and the voxel FEA) is the primary model.
