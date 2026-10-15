"""Draw the P rectangles with the C markers, and check how the box relates to the marked patch.

Run from the project folder (score_deks.py must be in the same folder):

    python p_box_check.py "SAMPLE - DEKS-EYEBANK/Reader1"        # table for every field
    python p_box_check.py "SAMPLE - DEKS-EYEBANK/Reader1" 2      # also draw image 2

P<k>.dat is read as: left, top, width, height (pixels). The table tests that
reading: if it is right, nearly every marker falls inside its own box.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv
from scipy.spatial import ConvexHull, cKDTree

import score_deks as sd

root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
draw_image = sys.argv[2] if len(sys.argv) > 2 else None

rows, drawn, extra = [], [], []
print(f"{'frame':>6}{'box (x, y, w, h)':>26}{'markers':>9}{'inside %':>10}"
      f"{'min gap px':>12}{'hull/box %':>12}")
for img, k, a_path, c_path in sd.frames(root):
    p_path = c_path.with_name(f"P{k}.dat")
    if not p_path.exists():
        continue
    x, y, w, h = np.array(p_path.read_text().split(), float)[:4]
    pts, _, _ = sd.read_C(c_path)
    if len(pts) < 4:
        continue
    inside = ((pts[:, 0] >= x) & (pts[:, 0] <= x + w)
              & (pts[:, 1] >= y) & (pts[:, 1] <= y + h))
    hull = ConvexHull(pts)
    hv = pts[hull.vertices]
    gap = min(hv[:, 0].min() - x, x + w - hv[:, 0].max(),
              hv[:, 1].min() - y, y + h - hv[:, 1].max())
    rows.append((inside.mean(), gap, hull.volume / (w * h)))
    pitch = np.median(cKDTree(pts).query(pts, k=2)[0][:, 1])
    extra.append((f"{img}.{k}", gap / pitch, hull.volume / (w * h)))
    print(f"{img + '.' + k:>6}{f'({x:.0f}, {y:.0f}, {w:.0f}, {h:.0f})':>26}{len(pts):>9}"
          f"{100 * inside.mean():>10.0f}{gap:>12.0f}{100 * hull.volume / (w * h):>12.0f}")
    if draw_image is not None and img == draw_image:
        drawn.append((k, (x, y, w, h), pts, hv))

if not rows:
    sys.exit("No P files found next to the C files.")
R = np.array(rows)
print(f"\n{len(R)} fields. Median share of markers inside the box: {100 * np.median(R[:, 0]):.0f}%. "
      f"Fields with every marker inside: {int((R[:, 0] == 1).sum())}. "
      f"Median gap from marked patch to box edge: {np.median(R[:, 1]):.0f} px. "
      f"Median marked area as share of box: {100 * np.median(R[:, 2]):.0f}%.")

# Does the Num error depend on how close the marked patch comes to the box edge?
sc = Path("scores.csv")
if sc.exists():
    err = {r["frame"]: int(r["num"]) - int(r["num_gt"]) for r in csv.DictReader(open(sc))}
    use = [(f, g, rt, err[f]) for f, g, rt in extra if f in err]
    if len(use) > 5:
        g = np.array([u[1] for u in use]); rt = np.array([u[2] for u in use]); e = np.array([u[3] for u in use])
        print(f"\nNum error (pred - truth) on {len(use)} scored fields, from scores.csv:")
        print(f"  correlation with gap to box (in cell spacings): {np.corrcoef(g, e)[0, 1]:+.2f}")
        print(f"  correlation with marked area / box area:        {np.corrcoef(rt, e)[0, 1]:+.2f}")
        for name, m in (("patch within 1 cell spacing of box edge", g < 1.0),
                        ("patch farther than 1 cell spacing", g >= 1.0)):
            if m.any():
                print(f"  {name}: {int(m.sum())} fields, mean error {e[m].mean():+.2f}, "
                      f"mean |error| {np.abs(e[m]).mean():.2f}")
    else:
        print("\nscores.csv has too few matching fields; run score_deks.py first.")
else:
    print("\nNo scores.csv here; run score_deks.py first to compare with the Num error.")

if drawn:
    fig, ax = plt.subplots(figsize=(10, 8))
    for k, (x, y, w, h), pts, hv in drawn:
        ax.add_patch(patches.Rectangle((x, y), w, h, linewidth=2, edgecolor="blue", facecolor="none"))
        ax.text(x + 4, y + 16, f"P{k}  ({x:.0f}, {y:.0f})", fontsize=10, color="red")
        ax.scatter(pts[:, 0], pts[:, 1], s=8, color="black")
        ax.plot(*np.vstack([hv, hv[:1]]).T, color="green", linewidth=1)
    ax.autoscale()
    ax.set_aspect("equal")
    ax.invert_yaxis()                      # image coordinates: y grows downward
    ax.set_title(f"Image {draw_image}: P boxes (blue), markers (black), marker hull (green)")
    out = f"p_boxes_image{draw_image}.png"
    fig.savefig(out, dpi=120)
    print(f"saved {out}")
    plt.show()
