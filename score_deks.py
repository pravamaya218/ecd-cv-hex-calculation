"""Score the marker-based CD / CV / HEX / Num model against the instrument's A*.dat values.

Run from the project folder (needs only numpy and scipy):

    python score_deks.py "SAMPLE - DEKS-EYEBANK/Reader1"
    python score_deks.py "SAMPLE - DEKS-EYEBANK/Reader1" --um 0.7716049

File formats assumed (checked against your listing, see the diagnostics it prints):
    A<k>.dat   first line = number of cells, then that many areas (um^2).
               If a further block of n small integers follows, it is read as sides.
    C<k>.dat   first line = number of markers, then x and y on alternating lines (pixels).
    P<k>.dat   four numbers (left, top, width, height) - not used for scoring.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, Delaunay, cKDTree

UM_DEFAULT = 0.768       # um per pixel, calibrated on DEKS Reader1 (the original run_all.py used 0.7716049)
CLEARANCE_PX = 2.0       # peel filter: how close to the hull counts as "on it"
MIN_PEELED = 99          # peel filter OFF (99). The original value was 3; off scored better on DEKS
MIN_EDGE_PX = 0.0        # merge polygon edges shorter than this (try 0.5, 1, 2 for HEX)
KNN = 16                 # neighbours used to cut each cell
BOX = 1e5                # half-size of the starting square for each cell
MIN_AREAS, MIN_POINTS = 20, 10
ALPHA_K = None           # None = convex hull only; e.g. 2.2 also rejects cells next to gaps


# ------------------------------------------------------------------ reading
def _numbers(path):
    return np.array(Path(path).read_text().split(), float)


def read_A(path):
    t = _numbers(path)
    n = int(t[0])
    vals = t[1:]
    areas = vals[:n]
    tail = vals[n:n + n]
    sides = None
    if (len(tail) == n and n > 0 and np.all(tail == np.round(tail))
            and tail.min() >= 3 and tail.max() <= 12):
        sides = tail.astype(int)
    return areas, sides, len(vals) - n


def read_C(path):
    t = _numbers(path)
    n = int(t[0])
    v = t[1:]
    m = min(n, len(v) // 2)
    return v[:2 * m].reshape(m, 2), n, len(v)


# ----------------------------------------------------------------- geometry
def clip(poly, nrm, off):
    """Keep the part of the polygon where nrm . x <= off."""
    d = poly @ nrm - off
    out = []
    for i in range(len(poly)):
        j = (i + 1) % len(poly)
        a, b, da, db = poly[i], poly[j], d[i], d[j]
        if da <= 0:
            out.append(a)
        if (da < 0 and db > 0) or (da > 0 and db < 0):
            out.append(a + (b - a) * (da / (da - db)))
    return np.array(out) if out else np.empty((0, 2))


def area(p):
    if len(p) < 3:
        return 0.0
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def dedupe(p, tol):
    if len(p) < 3:
        return p
    keep = [p[0]]
    for q in p[1:]:
        if np.hypot(*(q - keep[-1])) > tol:
            keep.append(q)
    if len(keep) > 1 and np.hypot(*(keep[0] - keep[-1])) <= tol:
        keep.pop()
    return np.array(keep)


def power_cells(pts, w):
    """Cell for every marker. w = 0 gives plain Voronoi; otherwise a power diagram."""
    k = min(KNN + 1, len(pts))
    _, nbr = cKDTree(pts).query(pts, k=k)
    c = pts.mean(0)
    base = np.array([[c[0] - BOX, c[1] - BOX], [c[0] + BOX, c[1] - BOX],
                     [c[0] + BOX, c[1] + BOX], [c[0] - BOX, c[1] + BOX]])
    cells = []
    for i, p in enumerate(pts):
        poly = base
        for j in nbr[i]:
            if j == i:
                continue
            nrm = pts[j] - p
            off = nrm @ (p + pts[j]) / 2 + (w[i] - w[j]) / 2
            poly = clip(poly, nrm, off)
            if len(poly) < 3:
                break
        cells.append(dedupe(poly, max(1e-6, MIN_EDGE_PX)))
    return cells


def interior(pts, cells, eps=1e-7):
    """True where every corner of the cell is inside the convex hull of all markers."""
    eq = ConvexHull(pts).equations
    ok = np.zeros(len(pts), bool)
    for i, c in enumerate(cells):
        if len(c) >= 3:
            ok[i] = bool(((c @ eq[:, :2].T + eq[:, 2]) <= eps).all())
    return ok


def alpha_mask(pts, k):
    """False for every marker touching a Delaunay triangle whose longest edge > k * pitch.

    Those long triangles span gaps (folds, dark areas, concave patch edges), so a
    marker touching one is on the real boundary of the marked patch.
    """
    tri = Delaunay(pts)
    pitch = np.median(cKDTree(pts).query(pts, k=2)[0][:, 1])
    P = pts[tri.simplices]
    longest = np.linalg.norm(P - np.roll(P, 1, axis=1), axis=2).max(axis=1)
    ok = np.ones(len(pts), bool)
    ok[np.unique(tri.simplices[longest > k * pitch])] = False
    return ok


def keep_cells(pts, cells):
    """Hull peel, then drop cells that sit on the hull with >= MIN_PEELED peeled neighbours."""
    inside = interior(pts, cells)
    if ALPHA_K is not None:
        inside &= alpha_mask(pts, ALPHA_K)
    idx = np.flatnonzero(inside)
    if len(idx) == 0:
        return idx
    eq = ConvexHull(pts).equations
    clear = np.array([-(cells[i] @ eq[:, :2].T + eq[:, 2]).max() for i in idx])
    indptr, indices = Delaunay(pts).vertex_neighbor_vertices
    kept = set(idx.tolist())
    peeled = np.array([sum(1 for j in indices[indptr[i]:indptr[i + 1]] if j not in kept)
                       for i in idx])
    return idx[~((clear < CLEARANCE_PX) & (peeled >= MIN_PEELED))]


# -------------------------------------------------------------------- model
def _base(pts):
    """Plain Voronoi cells, kept cells, areas, typical spacing, cells whose area can be trusted."""
    cl = power_cells(pts, np.zeros(len(pts)))
    keep = keep_cells(pts, cl)
    Vall = np.array([area(c) for c in cl])
    pitch = np.median(cKDTree(pts).query(pts, k=2)[0][:, 1])
    ok = interior(pts, cl) & (Vall > 0)
    return cl, keep, Vall, pitch, ok


def local_weights(pts, Vall, ok, pitch, coef):
    """w_i = coef * pitch^2 * log(V_i / mean area of its neighbours): large next to small."""
    ptr, ind = Delaunay(pts).vertex_neighbor_vertices
    Vg = np.where(ok, Vall, np.nan)
    w = np.zeros(len(pts))
    for i in np.flatnonzero(ok):
        nb = Vg[ind[ptr[i]:ptr[i + 1]]]
        nb = nb[np.isfinite(nb)]
        if len(nb):
            w[i] = coef * pitch ** 2 * np.log(Vall[i] / nb.mean())
    return w


def cv_from_weights(pts, w, um):
    """Rebuild the diagram with weights w, keep interior cells, return (CV, number kept)."""
    cw = power_cells(pts, w)
    kw = keep_cells(pts, cw)
    if len(kw) < 3:
        return np.nan, len(kw)
    A = np.array([area(cw[i]) for i in kw]) * um ** 2
    return A.std(ddof=0) / A.mean(), len(kw)


def local_sweep(pts, um, coefs):
    """CV of the local-weight diagram for each coefficient in coefs."""
    cl, keep, Vall, pitch, ok = _base(pts)
    return [cv_from_weights(pts, local_weights(pts, Vall, ok, pitch, c), um)[0] for c in coefs]


def analyse(pts, um, coef=1.0):
    cl, keep, Vall, pitch, ok = _base(pts)
    if len(keep) < 5:
        return None
    V = Vall[keep]
    sides = np.array([len(cl[i]) for i in keep])
    ave_px = V.mean()
    ave = ave_px * um ** 2

    # global weights (original rule): compare every cell with the frame average
    w = np.zeros(len(pts))
    w[ok] = pitch ** 2 * np.log(Vall[ok] / ave_px)
    cv, num_w = cv_from_weights(pts, w, um)
    # local weights: compare every cell with its own neighbours
    cv_local, _ = cv_from_weights(pts, local_weights(pts, Vall, ok, pitch, coef), um)

    hull = ConvexHull(pts)
    return dict(num=len(keep), num_w=num_w, ave=ave, cd=1e6 / round(ave),
                n_hull=len(hull.vertices), perim_pitch=hull.area / pitch, pitch_px=pitch,
                cv=cv, cv_local=cv_local, cv_plain=V.std(ddof=0) / V.mean(),
                hex=float((sides == 6).mean()))


# ------------------------------------------------------------------ scoring
def frames(root):
    dirs = sorted({p.parent for p in root.rglob("A*.dat")},
                  key=lambda d: (int(d.name) if d.name.isdigit() else 10 ** 9, str(d)))
    for d in dirs:
        afiles = sorted(d.glob("A*.dat"),
                        key=lambda p: int(p.stem[1:]) if p.stem[1:].isdigit() else 0)
        for a in afiles:
            c = d / f"C{a.stem[1:]}.dat"
            if c.exists():
                yield d.name, a.stem[1:], a, c


def stat(rows, key):
    p = np.array([r[key] for r in rows], float)
    t = np.array([r[key.split("_")[0] + "_gt"] for r in rows], float)
    m = np.isfinite(p) & np.isfinite(t)
    if not m.any():
        return None
    e = (p[m] - t[m]) / t[m]
    return 100 * e.mean(), 100 * np.abs(e).mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default="SAMPLE - DEKS-EYEBANK/Reader1")
    ap.add_argument("--um", type=float, default=UM_DEFAULT)
    ap.add_argument("--alpha", type=float, default=None,
                    help="concave-boundary threshold in pitches (try 2.0, 2.2, 2.5, 3)")
    ap.add_argument("--local-coef", type=float, default=1.0,
                    help="strength of the neighbour-based weights (see --sweep)")
    ap.add_argument("--sweep", action="store_true",
                    help="try several local-weight strengths; choose each one without the held-out image")
    args = ap.parse_args()
    global ALPHA_K
    ALPHA_K = args.alpha
    root = Path(args.root)

    rows, shown, store = [], 0, []
    for img, k, a_path, c_path in frames(root):
        areas, sides, extra = read_A(a_path)
        pts, n_decl, n_vals = read_C(c_path)
        if shown < 3:
            print(f"[check] image {img}.{k}: A declares {len(areas)} areas, "
                  f"{extra} further values (sides found: {sides is not None}); "
                  f"C declares {n_decl} markers, {n_vals} values -> {len(pts)} points")
            shown += 1
        if len(areas) < MIN_AREAS or len(pts) < MIN_POINTS:
            continue
        r = analyse(pts, args.um, args.local_coef)
        if r is None:
            continue
        rows.append(dict(
            frame=f"{img}.{k}", n_markers=len(pts),
            n_hull=r["n_hull"], perim_pitch=r["perim_pitch"], pitch_px=r["pitch_px"],
            num=r["num"], num_gt=len(areas), num_w=r["num_w"],
            cd=r["cd"], cd_gt=1e6 / areas.mean(),
            cv=r["cv"], cv_local=r["cv_local"], cv_gt=areas.std(ddof=0) / areas.mean(),
            cv_plain=r["cv_plain"],
            hex=r["hex"], hex_gt=float((sides == 6).mean()) if sides is not None else np.nan))
        store.append((img, pts))

    if not rows:
        print("No frames scored. Check the folder path and the [check] lines above.")
        return
    out = Path("scores.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"\n{len(rows)} frames scored (um/px = {args.um}, alpha = {args.alpha})")
    print(f"{'':9}{'bias %':>9}{'MAPE %':>9}")
    for key, label in (("cd", "CD"), ("cv", "CV glob"), ("cv_local", "CV local"), ("hex", "HEX")):
        s = stat(rows, key)
        print(f"{label:9}" + (f"{s[0]:>+9.1f}{s[1]:>9.1f}" if s else "   no ground truth"))
    n = np.array([r["num"] for r in rows])
    t = np.array([r["num_gt"] for r in rows])
    print(f"Num exact {int((n == t).sum())}/{len(rows)}, within 2: {int((abs(n - t) <= 2).sum())}/{len(rows)}")
    err = {}
    for a, b in zip(n - t, [0] * len(n)):
        err[int(a)] = err.get(int(a), 0) + 1
    print("Num error counts (pred - truth):", sorted(err.items()))
    d = (n - t).astype(float)
    print("corr(Num error, feature):", {
        f: round(float(np.corrcoef([r[f] for r in rows], d)[0, 1]), 2)
        for f in ("n_markers", "n_hull", "perim_pitch", "pitch_px")})
    print(f"weighted set differs from plain set on {sum(r['num_w'] != r['num'] for r in rows)} frames")
    print(f"wrote {out.resolve()}")
    print(f"(CV glob = weights from each cell vs the frame average; "
          f"CV local = vs its neighbours, strength {args.local_coef})")
    if args.sweep:
        sweep(store, rows, args.um)


def sweep(store, rows, um, coefs=(0.5, 1.0, 1.5, 2.0, 3.0, 4.0)):
    """CV error for several local-weight strengths, plus an honest leave-one-image-out score."""
    gt = np.array([r["cv_gt"] for r in rows])
    imgs = np.array([img for img, _ in store])
    M = np.array([local_sweep(pts, um, coefs) for _, pts in store])
    E = (M - gt[:, None]) / gt[:, None]
    print(f"\nLocal-weight strength sweep ({len(rows)} frames):")
    print(f"{'strength':>9}{'bias %':>9}{'MAPE %':>9}")
    for j, c in enumerate(coefs):
        print(f"{c:>9}{100 * np.nanmean(E[:, j]):>+9.1f}{100 * np.nanmean(np.abs(E[:, j])):>9.1f}")
    g = np.array([r["cv"] for r in rows])
    print(f"{'global':>9}{100 * np.mean((g - gt) / gt):>+9.1f}{100 * np.mean(np.abs(g - gt) / gt):>9.1f}   (current weights)")
    picked, errs = [], []
    for im in sorted(set(imgs), key=lambda s: (len(s), s)):
        tr, te = imgs != im, imgs == im
        j = int(np.nanargmin(np.nanmean(np.abs(E[tr]), axis=0)))
        picked.append(coefs[j])
        errs.extend(np.abs(E[te, j]))
    print(f"Strength chosen without the held-out image, then scored on it: MAPE {100 * np.nanmean(errs):.1f}% "
          f"(strengths picked: {sorted(set(picked))})")


if __name__ == "__main__":
    main()