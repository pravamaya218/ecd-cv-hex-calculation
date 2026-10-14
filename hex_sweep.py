"""Why is HEX low? Sweep the minimum wall length and compare side-count histograms.

Run from the project folder (score_deks.py must be in the same folder):

    python hex_sweep.py "SAMPLE - DEKS-EYEBANK/Reader1"
"""
import sys
from pathlib import Path

import numpy as np

import score_deks as sd

root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
sd.MIN_PEELED = 99          # the peel filter stays off, as in the best DEKS run

data = []
for img, k, a_path, c_path in sd.frames(root):
    areas, sides, _ = sd.read_A(a_path)
    pts, _, _ = sd.read_C(c_path)
    if len(areas) >= sd.MIN_AREAS and len(pts) >= sd.MIN_POINTS and sides is not None:
        data.append((pts, sides))


def pct(s, k):
    return 100 * float((s == k).mean())


gt_all = np.concatenate([s for _, s in data])
print(f"{len(data)} frames. Instrument sides (all cells): "
      + "  ".join(f"{k}:{pct(gt_all, k):.1f}%" for k in (4, 5, 6, 7, 8)))
print(f"\n{'MIN_EDGE_PX':>12}{'HEX bias %':>12}{'HEX MAPE %':>12}   ours 4 / 5 / 6 / 7 / 8 (%)")

for me in (0.0, 0.5, 1.0, 2.0, 3.0):
    sd.MIN_EDGE_PX = me
    err, ours = [], []
    for pts, sides in data:
        cl = sd.power_cells(pts, np.zeros(len(pts)))
        keep = sd.keep_cells(pts, cl)
        s = np.array([len(cl[i]) for i in keep])
        ours.append(s)
        err.append(((s == 6).mean() - (sides == 6).mean()) / (sides == 6).mean())
    err = np.array(err)
    o = np.concatenate(ours)
    print(f"{me:>12}{100 * err.mean():>+12.1f}{100 * np.abs(err).mean():>12.1f}   "
          + " / ".join(f"{pct(o, k):.1f}" for k in (4, 5, 6, 7, 8)))
