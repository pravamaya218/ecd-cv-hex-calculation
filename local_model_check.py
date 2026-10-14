"""One-parameter CV rule:  area_i = V_i * (V_i / mean area of neighbours)^gamma

V_i is our plain Voronoi area. gamma is fitted on frames where the A file lines
up with our cells (equal counts, area correlation >= 0.7), leaving out the whole
image being scored. The held-out image is then scored on EVERY frame, including
frames whose counts do not match, so this is a fair comparison with the current
weighted model over all frames.

Run from the project folder (score_deks.py must be in the same folder):

    python local_model_check.py "SAMPLE - DEKS-EYEBANK/Reader1"
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay

import new.score_deks as sd

UM = 0.768
root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
sd.MIN_PEELED = 99

fr = []
for img, k, a_path, c_path in sd.frames(root):
    areas, sides, _ = sd.read_A(a_path)
    pts, _, _ = sd.read_C(c_path)
    if len(areas) < sd.MIN_AREAS or len(pts) < sd.MIN_POINTS:
        continue
    cl = sd.power_cells(pts, np.zeros(len(pts)))
    keep = sd.keep_cells(pts, cl)
    if len(keep) < 5:
        continue
    V = np.array([sd.area(cl[i]) for i in keep])
    Vfull = np.full(len(pts), np.nan)
    Vfull[keep] = V
    ptr, ind = Delaunay(pts).vertex_neighbor_vertices
    loc = []
    for i, v in zip(keep, V):
        nb = Vfull[ind[ptr[i]:ptr[i + 1]]]
        nb = nb[np.isfinite(nb)]
        loc.append(nb.mean() if len(nb) else v)
    loc = np.array(loc)
    aligned = len(keep) == len(areas) and np.corrcoef(V, areas)[0, 1] >= 0.7
    fr.append(dict(img=img, name=f"{img}.{k}", V=V, loc=loc, aligned=aligned,
                   A=areas.astype(float), cv_gt=areas.std() / areas.mean(),
                   cv_now=sd.analyse(pts, UM)["cv"]))

al = [f for f in fr if f["aligned"]]
print(f"{len(fr)} frames scored, {len(al)} aligned frames used for fitting, "
      f"{len({f['img'] for f in fr})} images")


def fit_gamma(frames):
    x2, z = [], []
    for f in frames:
        x1 = np.log(f["V"] / f["V"].mean())
        x2.append(np.log(f["V"] / f["loc"]))
        z.append(np.log(f["A"] / f["A"].mean()) - x1)
    return float(np.polyfit(np.concatenate(x2), np.concatenate(z), 1)[0])


def cv_local(f, g):
    a = f["V"] * (f["V"] / f["loc"]) ** g
    return a.std() / a.mean()


print(f"gamma fitted on all aligned frames: {fit_gamma(al):.2f}\n")
err_now, err_loc = [], []
print(f"{'held-out image':>15}{'gamma':>8}")
for im in sorted({f["img"] for f in fr}, key=lambda s: (len(s), s)):
    train = [f for f in al if f["img"] != im]
    if not train:
        continue
    g = fit_gamma(train)
    print(f"{im:>15}{g:>8.2f}")
    for f in (f for f in fr if f["img"] == im):
        err_now.append((f["cv_now"] - f["cv_gt"]) / f["cv_gt"])
        err_loc.append((cv_local(f, g) - f["cv_gt"]) / f["cv_gt"])

print(f"\nCV error over all {len(err_now)} frames, leave-one-image-out")
print(f"{'':28}{'bias %':>8}{'MAPE %':>8}")
for name, e in (("current weighted model", err_now), ("local-contrast rule", err_loc)):
    e = np.array(e)
    print(f"{name:28}{100 * e.mean():>+8.1f}{100 * np.abs(e).mean():>8.1f}")