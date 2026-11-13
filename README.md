# Explanation — CD, CV, HEX and Num from marked cell centres

Input is a set of marked cell centres and a pixel scale. Output is the four
numbers the CellChek D panel prints. Scored on **37 DEKS eyebank fields**
(Reader1) against the instrument's own per-cell data. The scale, the peel
filter and the local-weight strength were chosen on these same fields, so the
scores are not fully held out (see Limits).

```python
from score_deks import analyse
r = analyse(points_xy, um_per_px)    # -> num, ave, cd, cv, cv_local, hex  (None if < 5 cells kept)
```

```
python score_deks.py "SAMPLE - DEKS-EYEBANK/Reader1"            # score every field, write scores.csv
python score_deks.py "SAMPLE - DEKS-EYEBANK/Reader1" --sweep    # also sweep the local-weight strength
```

|  | CD | CV | HEX | fields |
| --- | ---: | ---: | ---: | ---: |
| CellChek D exports  | 0.2% | 6.1% | 3.7% | 2 |
| DEKS Reader1, original settings (scale 0.7716, extra filter on) | 1.0% | 5.4% | 4.0% | 37 |
| DEKS Reader1, current settings, global weights | **0.3%** | 5.5% | **3.9%** | 37 |
| DEKS Reader1, current settings, local weights | **0.3%** | **4.9%** | **3.9%** | 37 |

MAPE. Bias with current settings: CD +0.0%, CV +1.9% (global) / +0.8% (local),
HEX −2.3%. Num is exact on 18/37 with the original settings and **21/37** with the
current ones, and within two cells on 35/37 either way. The "original settings"
row reproduces the original README's DEKS figures from a rewrite of the helper
functions, so the two code bases agree on that data.

## The method

**1 · Tessellate.** A plain Voronoi cell for every marker. Area by the shoelace
formula on the cell's corners.

**2 · Peel.** Keep a cell only if *all* of its corners fall inside the convex
hull of every marker. Cells that cross it are dropped whole, area included —
never clipped. The original method then dropped a survivor that is both within
2 px of the hull *and* has three or more wall-neighbours peeled. **On DEKS that
second rule is switched off** (`MIN_PEELED = 99`): with it off, Num is exact on
21/37 fields instead of 18/37 and the total absolute Num error falls from 33 to
25, in both halves of the images. The rule is kept as an option because the
original README introduced it to fix `Num` on the five CellChek panel frames,
which are not in this data set and have not been re-tested.

**3 · Num, Ave, CD and HEX, from the plain Voronoi.**

```
Num = number of cells kept
Ave = mean kept-cell area, µm²              (pixels × 0.768 µm/px)
CD  = 1e6 / round(Ave)
HEX = fraction of kept cells with exactly six walls
```

**4 · Re-tessellate with weights, for CV only.** Two weight rules are implemented.

```
pitch = median nearest-neighbour distance between markers

global (original):  w_i = pitch² · log(V_i / Ave)               markers inside the hull
local  (new):       w_i = pitch² · log(V_i / mean V of the neighbours of i)
                    w_i = 0                                      markers on the rim

wall between i and j:  n·x = n·(p_i+p_j)/2 + (w_i−w_j)/2,   n = p_j − p_i
```

Each wall moves off the midpoint toward the *smaller* cell. The global rule
compares every cell with the frame average; the local rule compares it with the
cells around it.

**5 · SD and CV, from the weighted cells.**

```
SD = sqrt( Σ (A_i − mean)² / n )      divide by n, not n−1
CV = SD / mean
```

The population divisor and the pentagon test (six markers, panel printed
`SD 0`) are from the original README and were not re-run on this data.

## Why the weights

**The topology is a plain Voronoi.** On the 37 DEKS fields the instrument's
side-count shares are 4:0.5%, 5:21.9%, 6:60.1%, 7:16.6%, 8:0.9%, against ours
0.4 / 21.8 / 59.0 / 17.8 / 1.0. On the fields where the A file lines up with our
cells, the side counts agree on 96% of cells (median). Merging short walls to
fix the small HEX gap made it worse at every threshold, so no wall-length rule is
used.

**The areas are not.** The A files list the instrument's per-cell areas in
marker order (median correlation with our Voronoi areas 0.90 on 19 fields).
Against those, instrument areas are about twice as spread out as plain Voronoi
areas: regression slope median 1.76, CV ratio about 2. That gap is what the
weights have to reproduce.

**What predicts the stretch.** Regressing `log(instrument area / frame mean)` on
our cell features, pooled over the 19 aligned fields (1,361 cells):

| Features | R² | Coefficients |
| --- | ---: | --- |
| own area only | 0.79 | own = 1.81 |
| + number of sides | 0.90 | own = 1.45, sides−6 = +0.19 |
| + own / neighbour area | 0.89 | own = 1.03, neighbour term = 1.33 |
| + both | 0.92 | own = 1.07, sides−6 = +0.13, neighbour term = 0.83 |

With the neighbour term in, the exponent on a cell's own size is about 1 and the
extra stretch comes from how much bigger the cell is than its neighbours. This is
why the local rule exists. A one-parameter version,
`area = V · (V / mean neighbour V)^γ`, gives γ = 1.36 (1.34 to 1.38 with each
image left out) and a CV error of 4.8%, against 5.5% for the global weights, with
the strength chosen without the image being scored. The sides term has the same
sign as the `0.045(n−6)` term in the original weight regression, which the
original code did not use.

**How much to trust that.** The local weights are closer to the instrument's CV
than the global ones on only **18 of 37** fields. Split by patch size, on the 31
fields with 70 or more markers both give about 4.2%; on the 6 fields with fewer
than 70 markers the global weights give 12.1% and the local ones 7.8%. The local
rule's advantage rests on a few small-patch fields, so it is not shown to be a
better model than the global one.

**The rim weights are muted for the same reason as before.** A marker on the edge
of the marked patch has a one-sided neighbourhood, so its cell runs off
unbounded and `log(V/Ave)` hands it a large positive weight; both rules set rim
weights to zero. The original README measured CV error going from 8% to 48%
without this; that was not re-measured here.

The following claims are from the original README (CellChek exports) and have
**not** been re-tested on this data: the walls sitting about 0.7 µm off the
midpoint, the weight regression `w = 0.026 + 0.045(n−6) + 1.071·log(V/Ave)`
with R² 0.47, and hexagons sitting at weight +0.007.

## What else was checked

- **Scale.** At 0.7716 µm/px CD was low by 0.9% on every field (MAPE 1.0%). At
  0.768 the bias is 0.0% and MAPE 0.3%. This is a calibration chosen on the data.
- **Concave boundary** (cells next to gaps rejected): no gain at any threshold.
- **Num error follows patch size.** It correlates +0.53 with the hull perimeter
  (in cell spacings) and +0.44 with marker count, not with the number of hull
  corners (−0.13). The instrument removes about 0.91 cells per cell spacing of
  perimeter (−2.5); we remove 0.84 (+0.3) with the filter off. Remaining errors
  are mostly ±1 cell.
- **The P box does not explain Num.** Each `P` file is the analysis box (left,
  top, width, height); every marker lies inside it in all 39 fields. The box
  size varies (13 fields differ in width). Fields whose patch comes within one
  cell spacing of the box edge have mean |Num error| 0.69, against 0.67 for the
  others.

## Limits

- The scale, the peel filter and γ were chosen on the fields they are scored on.
  Leave-one-image-out checks were used for the CV rules, but the sample is small
  (37 fields, 14 images, and fields from one image are not independent).
- Nothing here was run on the CellChek D exports. The "none of which was used to
  set anything" claim in the original README does not hold for the settings
  above.
- CV is still off by about 4% on larger patches and 8–12% on small ones, and
  about 10% of the variance in instrument areas is not explained by marker
  positions.
- Two fields with matching counts (2.3 and 7.3) do not line up cell by cell with
  the A file.

Reference for the cell model (Voronoi from marked centres, shoelace areas, cells
with incomplete outer rings excluded, ECD = inner cells / their total area): Li,
D., Li, Z., Duan, H. et al., *npj Digital Medicine* 8, 791 (2025),
https://doi.org/10.1038/s41746-025-02239-0. It does not describe a CV weighting.


*******************************************************************************

For Center Method:

Please run patent_model.py then run patent_check.py to verify how the center model is working

frame_metrics.ipynb notebook to test the performance in %bias and %MAPE for center method

*********************************************************************************

