"""Score the Center Method model (patent_model.py) against the instrument's A files.

Run from the project folder (score_deks.py and patent_model.py in the same folder):

    python patent_check.py "SAMPLE - DEKS-EYEBANK/Reader1"

Prints one row per setting, next to the current Voronoi model, then sweeps the outer-cell
angle for the best setting, and writes the best setting's per-field results to
patent_scores.csv.

Columns
    Num exact / +-2   fields where our count equals / is within 2 of the instrument's
    CD, CV, HEX       typical error (MAPE, %) against the A-file values; CV also shows bias
    aligned           fields where the counts are equal, so cells can be paired in marker order
    med r, med slope  median correlation and slope of instrument area against ours on those
                      fields (the Voronoi model gave about 0.90 and 1.76)
    cell err          median relative error of a single cell's area on those fields
"""
import csv
import sys
from pathlib import Path

import numpy as np

import score_deks as sd
import patent_model as pm

UM = 0.768
root = Path(sys.argv[1] if len(sys.argv) > 1 else "SAMPLE - DEKS-EYEBANK/Reader1")
sd.MIN_PEELED = 99

data = []
for img, k, a_path, c_path in sd.frames(root):
    areas, sides, _ = sd.read_A(a_path)
    pts, _, _ = sd.read_C(c_path)
    if len(areas) < sd.MIN_AREAS or len(pts) < sd.MIN_POINTS:
        continue
    data.append(dict(name=f"{img}.{k}", pts=pts, areas=areas, sides=sides, cache=pm.periph(pts)))
print(f"{len(data)} fields")


def run(fn):
    rows = []
    for d in data:
        r = fn(d)
        if r is None:
            continue
        a = d["areas"]
        row = dict(frame=d["name"], n_markers=len(d["pts"]), num=r["num"], num_gt=len(a),
                   cd=r["cd"], cd_gt=1e6 / a.mean(), cv=r["cv"], cv_gt=a.std(ddof=0) / a.mean(),
                   hex=r["hex"],
                   hex_gt=float((d["sides"] == 6).mean()) if d["sides"] is not None else np.nan,
                   r=np.nan, slope=np.nan, cell_err=np.nan)
        if "areas_um2" in r and r["num"] == len(a):
            ours = r["areas_um2"]
            row["r"] = float(np.corrcoef(ours, a)[0, 1])
            row["slope"] = float(np.polyfit(ours, a, 1)[0])
            row["cell_err"] = float(np.median(np.abs(ours - a) / a))
        rows.append(row)
    return rows


def stats(rows, k):
    p = np.array([r[k] for r in rows], float)
    g = np.array([r[k + "_gt"] for r in rows], float)
    m = np.isfinite(p) & np.isfinite(g)
    e = (p[m] - g[m]) / g[m]
    return 100 * e.mean(), 100 * np.abs(e).mean()


HEAD = (f"{'setting':24}{'Num exact':>10}{'+-2':>5}{'CD':>6}{'CV bias':>9}{'CV':>6}{'HEX':>6}"
        f"{'aligned':>9}{'med r':>7}{'slope':>7}{'cell err':>10}")


def line(label, rows):
    n = np.array([r["num"] for r in rows])
    t = np.array([r["num_gt"] for r in rows])
    rr = np.array([r["r"] for r in rows])
    ok = np.isfinite(rr)
    med = (lambda k: np.median([r[k] for r in rows if np.isfinite(r[k])])) if ok.any() else None
    print(f"{label:24}{int((n == t).sum()):>7}/{len(rows):<3}{int((abs(n - t) <= 2).sum()):>3}"
          f"{stats(rows, 'cd')[1]:>6.1f}{stats(rows, 'cv')[0]:>+9.1f}{stats(rows, 'cv')[1]:>6.1f}"
          f"{stats(rows, 'hex')[1]:>6.1f}{int(ok.sum()):>9}"
          + (f"{med('r'):>7.2f}{med('slope'):>7.2f}{100 * med('cell_err'):>9.1f}%" if med else f"{'-':>7}{'-':>7}{'-':>10}"))


print("\n" + HEAD)
base = run(lambda d: (lambda r: None if r is None else dict(num=r["num"], cd=r["cd"], cv=r["cv"], hex=r["hex"]))(
    sd.analyse(d["pts"], UM)))
line("Voronoi, global weights", base)
base_l = run(lambda d: (lambda r: None if r is None else dict(num=r["num"], cd=r["cd"], cv=r["cv_local"], hex=r["hex"]))(
    sd.analyse(d["pts"], UM)))
line("Voronoi, local weights", base_l)

best = None
for od in ("own", "mid"):
    for p in (1, 2, 3, 4, 5):
        rows = run(lambda d, p=p, od=od: pm.analyse_patent(d["pts"], UM, p, 100.0, od, cache=d["cache"]))
        line(f"Center, p={p}, outer={od}", rows)
        err = stats(rows, "cv")[1]
        if best is None or err < best[0]:
            best = (err, p, od, rows)

print(f"\nSame walls, but the cell is built by clipping a large square with every wall (construction = clip):")
print(HEAD)
for od in ("own", "mid"):
    for p in (2, 3, 4):
        rows = run(lambda d, p=p, od=od: pm.analyse_patent(d["pts"], UM, p, 100.0, od, "clip", cache=d["cache"]))
        line(f"Clip, p={p}, outer={od}", rows)

_, bp, bod, brows = best
print(f"\nBest CV error above: p = {bp}, outer = {bod}. Sweeping the outer-cell angle at that setting:")
print(HEAD)
for a_out in (90.0, 100.0, 110.0):
    rows = run(lambda d, a=a_out: pm.analyse_patent(d["pts"], UM, bp, a, bod, cache=d["cache"]))
    line(f"angle {a_out:.0f}", rows)

with open("patent_scores.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(brows[0]))
    w.writeheader()
    w.writerows(brows)
print(f"\nwrote patent_scores.csv for p = {bp}, outer = {bod}, angle 100 ({len(brows)} fields)")
print("The best setting was picked on these same fields, so treat its error as optimistic.")
