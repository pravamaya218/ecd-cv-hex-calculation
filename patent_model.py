"""Center Method cell outlines, after Konan patent US 5,523,213 (Takahashi and Kasahara, 1996).

Input: marker positions. Output: one polygon per counted cell. For every marker:

 1. Neighbour candidates = markers within radius_k times the distance to its nearest marker.
 2. Sort them by angle around the marker.
 3. Pass 1, successive pairs: if the angle between them is below a_pair degrees, drop the
    one farther from the marker. If an angle is above a_outer degrees the marker is an
    outer cell (not counted).
 4. Pass 2, successive triples: if the first and third are closer to each other than the
    middle one is to the marker, drop the middle one, unless the angle first-marker-third
    is above a_triple degrees.
 5. d_i = mean distance from the marker to its remaining neighbours.
 6. For each neighbour j, the wall is the line perpendicular to the segment i-j through the
    point that divides it in the ratio d_i^p : d_j^p (counted from i).
 7. Cell vertices = intersections of the walls of successive neighbours.

Parameter ranges from the patent: radius_k 2.5-3.5, a_pair 20-40, a_outer 90-110,
a_triple 120-140, p between 2 and 4 ("second to fourth power"). The patent does not give
exact values, so the defaults here are the values used in its worked example (3, 30, 100,
128) and p = 3.

Assumptions made where the patent text is not explicit (see patent_check.py):
 - angles are measured as the counter-clockwise gap between successive neighbours;
 - after a point is dropped the same position is tested again;
 - outer cells still get a full neighbour set and a mean distance d_i, because their
   neighbours need it for their walls (outer_d="own"); outer_d="mid" puts those walls
   halfway instead. Both are tried by patent_check.py.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def _clip(poly, nrm, off):
    """Keep the part of the polygon where nrm . x <= off (one straight cut)."""
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


def periph(pts, radius_k=3.0, a_pair=30.0, a_triple=128.0):
    """Neighbour sets for every marker. Returns (S, max_gap, d_bar).

    max_gap[i] is the largest angle (degrees) between successive neighbours seen in
    pass 1; a marker is an outer cell when max_gap exceeds a_outer, so a_outer can be
    changed afterwards without recomputing.
    """
    pts = np.asarray(pts, float)
    tree = cKDTree(pts)
    nn = tree.query(pts, k=2)[0][:, 1]
    S_all, max_gap = [], np.zeros(len(pts))
    for i, P in enumerate(pts):
        idx = [j for j in tree.query_ball_point(P, radius_k * nn[i]) if j != i]
        if len(idx) < 3:
            S_all.append(idx)
            max_gap[i] = 360.0
            continue
        ang = {j: np.arctan2(pts[j][1] - P[1], pts[j][0] - P[0]) for j in idx}
        dist = {j: float(np.hypot(*(pts[j] - P))) for j in idx}
        S = sorted(idx, key=lambda j: ang[j])

        def gap(a, b):
            return float(np.degrees(ang[b] - ang[a]) % 360.0)

        m, worst = 0, 0.0
        while len(S) >= 2 and m < len(S):                       # pass 1: pairs
            j = (m + 1) % len(S)
            a, b = S[m], S[j]
            g = gap(a, b)
            worst = max(worst, g)
            if g < a_pair:
                drop = m if dist[a] > dist[b] else j
                S.pop(drop)
                if drop < m:
                    m -= 1
            else:
                m += 1
        m = 0
        while len(S) >= 3 and m < len(S):                       # pass 2: triples
            L = len(S)
            a, b, c = S[m], S[(m + 1) % L], S[(m + 2) % L]
            if gap(a, c) > a_triple:
                m += 1
            elif float(np.hypot(*(pts[a] - pts[c]))) < dist[b]:
                drop = (m + 1) % L
                S.pop(drop)
                if drop < m:
                    m -= 1
            else:
                m += 1
        S_all.append(S)
        max_gap[i] = worst if len(S) >= 3 else 360.0
    d_bar = np.array([np.mean([np.hypot(*(pts[j] - pts[i])) for j in s]) if len(s) else nn[i]
                      for i, s in enumerate(S_all)])
    return S_all, max_gap, d_bar


def polygons(pts, S, max_gap, d_bar, p=3.0, a_outer=100.0, outer_d="own", construction="vertices"):
    """Cell polygon for every counted marker (None for outer or degenerate cells).

    construction: "vertices" takes the corners as the intersections of the walls of successive
    neighbours (the patent's wording). "clip" cuts a large square with every wall of the
    remaining neighbours, so the cell is the area they enclose.

    outer_d: how a wall against an outer cell is placed. "own" uses that cell's mean
    distance like any other; "mid" puts the wall halfway (ratio 1:1), because the patent
    does not compute a mean distance for outer cells.
    """
    pts = np.asarray(pts, float)
    is_outer = (max_gap > a_outer) | np.array([len(s) < 3 for s in S])
    out = []
    for i, s in enumerate(S):
        if max_gap[i] > a_outer or len(s) < 3:
            out.append(None)
            continue
        P = pts[i]
        lines = []
        for j in s:
            n = pts[j] - P
            di = d_bar[i] ** p
            dj = di if (outer_d == "mid" and is_outer[j]) else d_bar[j] ** p
            t = di / (di + dj)
            lines.append((n, float(n @ P + t * (n @ n))))
        if construction == "clip":
            poly = P + 1e4 * np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
            for n, c in lines:
                poly = _clip(poly, n, c)
                if len(poly) < 3:
                    break
            out.append(poly if len(poly) >= 3 else None)
            continue
        verts = []
        for k in range(len(s)):
            n1, c1 = lines[k]
            n2, c2 = lines[(k + 1) % len(s)]
            det = n1[0] * n2[1] - n1[1] * n2[0]
            if abs(det) < 1e-9:
                continue
            verts.append(((c1 * n2[1] - c2 * n1[1]) / det, (n1[0] * c2 - n2[0] * c1) / det))
        out.append(np.array(verts) if len(verts) >= 3 else None)
    return out


def shoelace(v):
    x, y = v[:, 0], v[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def analyse_patent(pts, um, p=3.0, a_outer=100.0, outer_d="own", construction="vertices", cache=None, **kw):
    """Num, mean area, CD, CV, HEX from the Center Method outlines.

    cache: pass the tuple returned by periph() to reuse it for several p or a_outer.
    Returns a dict that also holds idx (counted markers, ascending) and areas_um2 in
    that order, for cell-by-cell comparison with an A file.
    """
    pts = np.asarray(pts, float)
    S, max_gap, d_bar = cache if cache is not None else periph(pts, **kw)
    polys = polygons(pts, S, max_gap, d_bar, p, a_outer, outer_d, construction)
    idx = np.array([i for i, v in enumerate(polys) if v is not None])
    if len(idx) < 3:
        return None
    A = np.array([shoelace(polys[i]) for i in idx]) * um ** 2
    sides = np.array([len(polys[i]) for i in idx])
    ave = A.mean()
    return dict(num=len(idx), ave=ave, cd=1e6 / round(ave), cv=A.std(ddof=0) / ave,
                hex=float((sides == 6).mean()), idx=idx, areas_um2=A, sides=sides)
