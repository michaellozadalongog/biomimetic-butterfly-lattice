#!/usr/bin/env python3
"""Generate STLs, run the computational comparison and FEA, and build figures."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import analyze
import fea_bend
import generate_lattice
import plot_figures


def main() -> int:
    rc = generate_lattice.main([])
    if rc:
        return rc
    rc = analyze.main([])
    if rc:
        return rc
    rc = fea_bend.main([])
    if rc:
        return rc
    return plot_figures.main([])


if __name__ == "__main__":
    raise SystemExit(main())
