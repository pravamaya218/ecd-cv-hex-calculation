"""Does the Num error follow patch size, hull vertices, or hull perimeter?

Run from the project folder, after score_deks.py has written scores.csv:

    python boundary_check.py "SAMPLE - DEKS-EYEBANK/Reader1"
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
rows = list(csv.DictReader(open("scores.csv")))

feat = {"n_markers": [], "n_hull": [], "perim_pitch": [], "pitch_px": []}
err = []
for r in rows:
    img, k = r["frame"].split(".")
    c = next(p for p in root.rglob(f"C{k}.dat") if p.parent.name == img)
    t = np.array(c.read_text().split(), float)
    n = int(t[0])
    pts = t[1:1 + 2 * n].reshape(-1, 2)
    hull = ConvexHull(pts)
    pitch = np.median(cKDTree(pts).query(pts, k=2)[0][:, 1])
    feat["n_markers"].append(len(pts))
    feat["n_hull"].append(len(hull.vertices))
    feat["perim_pitch"].append(hull.area / pitch)   # 2-D hull.area is the perimeter
    feat["pitch_px"].append(pitch)
    err.append(int(r["num"]) - int(r["num_gt"]))

print(f"{len(rows)} frames; correlation of (Num pred - truth) with:")
for name, v in feat.items():
    print(f"  {name:12} {np.corrcoef(v, err)[0, 1]:+.2f}")


n_mk = np.array(feat["n_markers"], float)
perim = np.array(feat["perim_pitch"], float)
peel_gt = n_mk - np.array([int(r["num_gt"]) for r in rows])
peel_pred = n_mk - np.array([int(r["num"]) for r in rows])
for label, y in (("instrument", peel_gt), ("ours", peel_pred)):
    slope, icpt = np.polyfit(perim, y, 1)
    print(f"cells removed, {label}: {slope:.2f} x perimeter(pitches) {icpt:+.1f}")