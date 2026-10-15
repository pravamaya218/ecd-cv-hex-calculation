"""Do the instrument's per-cell areas line up with our cells, and how do they differ?

Run from the project folder (score_deks.py must be in the same folder):

    python percell_check.py "SAMPLE - DEKS-EYEBANK/Reader1"

Only frames where we keep exactly as many cells as the instrument are used,
because the A file does not say which marker each area belongs to. The test
assumes the areas are listed in marker order (kept markers, ascending index).
If that is wrong the correlations will be near zero and the output says so.
"""
import sys
from pathlib import Path

import numpy as np

import score_deks as sd

UM = 0.768
root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
sd.MIN_PEELED = 99

rows = []
for img, k, a_path, c_path in sd.frames(root):
    areas, sides, _ = sd.read_A(a_path)
    pts, _, _ = sd.read_C(c_path)
    if len(areas) < sd.MIN_AREAS or len(pts) < sd.MIN_POINTS or sides is None:
        continue
    cl = sd.power_cells(pts, np.zeros(len(pts)))
    keep = sd.keep_cells(pts, cl)
    if len(keep) != len(areas):
        continue
    ours = np.array([sd.area(cl[i]) for i in keep]) * UM ** 2
    s_ours = np.array([len(cl[i]) for i in keep])
    r = float(np.corrcoef(ours, areas)[0, 1])
    slope = float(np.polyfit(ours, areas, 1)[0])
    rows.append((f"{img}.{k}", len(areas), r, slope, float((s_ours == sides).mean()),
                 ours.std() / ours.mean(), areas.std() / areas.mean()))

if not rows:
    print("No frame has matching cell counts, so nothing can be compared.")
    sys.exit()

print(f"{'frame':>6}{'cells':>7}{'r(area)':>9}{'slope':>7}{'sides ok':>10}{'CV ours':>9}{'CV instr':>10}")
for f, n, r, sl, so, c1, c2 in rows:
    print(f"{f:>6}{n:>7}{r:>9.2f}{sl:>7.2f}{so:>10.2f}{c1:>9.3f}{c2:>10.3f}")
R = np.array([x[2] for x in rows])
print(f"\n{len(rows)} frames with equal counts. median r = {np.median(R):.2f}, "
      f"median slope = {np.median([x[3] for x in rows]):.2f}, "
      f"median sides agreement = {np.median([x[4] for x in rows]):.2f}")
print("r near 0 -> the A file is NOT in marker order; r above 0.7 -> cell-by-cell comparison is possible.")