"""What simple rule turns our Voronoi cell areas into the instrument's areas?

Run from the project folder (score_deks.py must be in the same folder):

    python stretch_check.py "SAMPLE - DEKS-EYEBANK/Reader1"

Uses only frames where the A file lines up with our kept cells (equal counts and
area correlation >= 0.7, as found by percell_check.py). CV scores are
leave-one-image-out: the rule is fitted on the other images, then applied to the
held-out image, so fields from the same image never help predict each other.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay

import score_deks as sd

UM = 0.768
root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
sd.MIN_PEELED = 99

fr = []
for img, k, a_path, c_path in sd.frames(root):
    areas, sides, _ = sd.read_A(a_path)
    pts, _, _ = sd.read_C(c_path)
    if len(areas) < sd.MIN_AREAS or len(pts) < sd.MIN_POINTS or sides is None:
        continue
    cl = sd.power_cells(pts, np.zeros(len(pts)))
    keep = sd.keep_cells(pts, cl)
    if len(keep) != len(areas):
        continue
    V = np.array([sd.area(cl[i]) for i in keep]) * UM ** 2
    if np.corrcoef(V, areas)[0, 1] < 0.7:
        continue
    Vfull = np.full(len(pts), np.nan)
    Vfull[keep] = V / UM ** 2
    ptr, ind = Delaunay(pts).vertex_neighbor_vertices
    loc = []
    for j, i in enumerate(keep):
        nb = Vfull[ind[ptr[i]:ptr[i + 1]]]
        nb = nb[np.isfinite(nb)]
        loc.append(nb.mean() if len(nb) else Vfull[i])
    loc = np.array(loc) * UM ** 2
    s6 = np.array([len(cl[i]) for i in keep]) - 6.0
    fr.append(dict(img=img, name=f"{img}.{k}", V=V, A=areas.astype(float), s6=s6,
                   x1=np.log(V / V.mean()), y=np.log(areas / areas.mean()),
                   x2=np.log(V / loc), cv_now=sd.analyse(pts, UM)["cv"],
                   cv_gt=areas.std() / areas.mean()))

print(f"{len(fr)} aligned frames, {sum(len(f['V']) for f in fr)} cells, "
      f"{len({f['img'] for f in fr})} images")
if len(fr) < 3:
    sys.exit("Too few aligned frames to fit anything.")


def design(f, cols):
    c = {"1": np.ones(len(f["V"])), "x1": f["x1"], "s6": f["s6"], "x2": f["x2"]}
    return np.column_stack([c[n] for n in cols])


MODELS = {
    "own area only":              ["1", "x1"],
    "+ number of sides":          ["1", "x1", "s6"],
    "+ own / neighbour area":     ["1", "x1", "x2"],
    "+ both":                     ["1", "x1", "s6", "x2"],
}


def fit(frames, cols):
    X = np.vstack([design(f, cols) for f in frames])
    y = np.concatenate([f["y"] for f in frames])
    return np.linalg.lstsq(X, y, rcond=None)[0], X, y


print("\nCell level, all aligned frames: log(instrument area / frame mean) regressed on")
for name, cols in MODELS.items():
    b, X, y = fit(fr, cols)
    r2 = 1 - ((y - X @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    print(f"  {name:26} R2 = {r2:.2f}   coefficients "
          + ", ".join(f"{c}={v:+.2f}" for c, v in zip(cols, b)))


def cv_of(a):
    return a.std() / a.mean()


imgs = sorted({f["img"] for f in fr})
res = {"current weighted model": [], "linear stretch": []}
res.update({f"power law, {n}": [] for n in MODELS})
for im in imgs:
    train = [f for f in fr if f["img"] != im]
    test = [f for f in fr if f["img"] == im]
    if not train:
        continue
    num = sum(((np.exp(f["x1"]) - 1) * (np.exp(f["y"]) - 1)).sum() for f in train)
    den = sum(((np.exp(f["x1"]) - 1) ** 2).sum() for f in train)
    s = num / den
    coefs = {n: fit(train, c)[0] for n, c in MODELS.items()}
    for f in test:
        e = lambda cv: abs(cv - f["cv_gt"]) / f["cv_gt"]
        res["current weighted model"].append(e(f["cv_now"]))
        res["linear stretch"].append(e(cv_of(1 + s * (np.exp(f["x1"]) - 1))))
        for n, c in MODELS.items():
            res[f"power law, {n}"].append(e(cv_of(np.exp(design(f, c) @ coefs[n]))))

print("\nFrame-level CV error, leave-one-image-out (MAPE %):")
for n, v in res.items():
    print(f"  {n:34} {100 * np.mean(v):5.1f}")