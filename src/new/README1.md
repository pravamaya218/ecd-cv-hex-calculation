This project works out how a CellChek D specular microscope turns the operator's marked cell centres into the four numbers it prints for the corneal endothelium, and tests how closely a marker-only model reproduces them.

Value	Meaning
Num	number of cells counted
CD	cell density, cells/mm² (1,000,000 / mean cell area in µm²)
CV	coefficient of variation of cell area (SD / mean, SD divided by n)
HEX	share of counted cells with exactly six sides

Status: a working model and a set of diagnostics, scored on 37 fields of the DEKS eyebank sample (Reader1). Several settings were chosen on that same data, so the numbers below are not fully held out (see Limitations).

Data

Each field has three files, in SAMPLE - DEKS-EYEBANK/Reader1/DAT_*/<image>/:

File	Content	How it was established
A<k>.dat	the instrument's answer: count n, then n cell areas (µm²), then n side counts	areas are certain; the side-count block was inferred from file sizes and is supported by 96% median agreement with our Voronoi side counts
C<k>.dat	count n, then the n marker positions as alternating x and y lines (pixels)	confirmed by all markers lying inside the P box
P<k>.dat	analysis box: left, top, width, height (pixels)	confirmed: 100% of markers fall inside the box in all 39 fields

k is the field number within an image (an image has two or three fields). The A files contain only the cells the instrument kept, in marker order, so the instrument's count is smaller than the number of markers.

Method
Voronoi. Each marker gets the region closer to it than to any other marker; the area is computed with the shoelace formula.
Peel. A cell is kept only if all its corners lie inside the convex hull of all markers.
Num, CD, HEX come from the kept cells: their count, 1e6 / round(mean area), and the share with six sides.
CV comes from a second, weighted tessellation (a power diagram) in which each wall is shifted toward the smaller cell, because a plain Voronoi gives roughly half the instrument's CV. Two weight rules are implemented:
global: each cell compared with the frame average, w = pitch² · log(V / mean V);
local: each cell compared with the mean of its neighbours, w = pitch² · log(V / mean V of neighbours).

Cells at the edge of the marked patch are excluded because their true area cannot be measured. This matches the published description of the approach (Li et al., 2025; see References).

Settings in force
Setting	Value	Why
Scale	0.768 µm/px	removed a systematic 0.9% CD bias that appeared at 0.7716
Extra peel filter (MIN_PEELED)	off (99)	scored better on DEKS; kept as an option because it was introduced for the CellChek panel frames, which are not in this folder
Minimum wall length	0 px	merging short walls made HEX worse
Concave boundary (--alpha)	off	no gain at any threshold
Results (37 fields, scale 0.768, peel filter off)
Value	Bias	Typical error (MAPE)
CD	+0.0%	0.3%
HEX	−2.3%	3.9%
CV, global weights	+1.9%	5.5%
CV, local weights (strength 1.0, not tuned)	+0.8%	4.9%
Num	exact on 21 of 37, within 2 on 35 of 37	

CV by patch size: on the 31 fields with 70 or more markers both weightings give about 4.2%; on the 6 fields with fewer than 70 markers, global gives 12.1% and local 7.8%. The local rule is closer to the instrument on 18 of 37 fields, so no consistent improvement over the global rule has been shown; its gain rests on a few small-patch fields.

What the diagnostics found
Topology is plain Voronoi. Side counts agree with the instrument on 96% of cells (median), and the side-count histograms are within about 1 point.
Areas are stretched. On fields where the A file lines up with our cells (19 fields, median area correlation 0.90), instrument areas are about twice as spread out as plain Voronoi areas (median slope 1.76; CV about 2×).
Neighbour contrast matters. Regressing log instrument area on log own area gives R² 0.79; adding the cell's size relative to its neighbours or its side count raises it to about 0.90; both give 0.92. With the neighbour term in, the exponent on own area is about 1.0 and the neighbour term about 1.3 (γ = 1.36, stable between 1.34 and 1.38 when each image is left out).
Num error follows patch size, not the box. The error correlates 0.53 with the hull perimeter (in cell spacings) and 0.44 with marker count, but not with hull corners (−0.13) or the P box (−0.16 with the gap to the box edge).
The extra peel filter lowered accuracy on DEKS (Num exact 18 → 21 and total absolute error 33 → 25 when off), with the improvement present in both halves of the images.
How to run

Run everything from the project folder (the scripts read paths relative to it). Needs numpy, scipy and matplotlib.

powershell
python score_deks.py "SAMPLE - DEKS-EYEBANK/Reader1"            # writes scores.csv
python score_deks.py "SAMPLE - DEKS-EYEBANK/Reader1" --sweep    # also sweeps local-weight strength
python hex_sweep.py "SAMPLE - DEKS-EYEBANK/Reader1"
python boundary_check.py "SAMPLE - DEKS-EYEBANK/Reader1"        # needs scores.csv
python percell_check.py "SAMPLE - DEKS-EYEBANK/Reader1"
python stretch_check.py "SAMPLE - DEKS-EYEBANK/Reader1"
python local_model_check.py "SAMPLE - DEKS-EYEBANK/Reader1"
python p_box_check.py "SAMPLE - DEKS-EYEBANK/Reader1" 2         # needs scores.csv; 2 = image to draw

score_deks.py must stay in the same folder as the other scripts, which import it.

Files
File	Role
score_deks.py	the model and the scoring; writes scores.csv (the only CSV any of these scripts writes)
hex_sweep.py	does merging short walls fix HEX? (no)
boundary_check.py	does the Num error follow patch size or boundary?
percell_check.py	is the A file in marker order, and how do instrument areas differ from Voronoi?
stretch_check.py	which per-cell rule maps Voronoi area to instrument area?
local_model_check.py	one-parameter neighbour-contrast CV rule, fitted leaving one image out
p_box_check.py	draws the P boxes with the markers; tests whether the box affects Num; saves p_boxes_image<N>.png
legacy/	earlier scripts and CSVs that pool all fields of an image or compare the data with itself; do not use for results (in particular generate_scores.py copies ground truth into the prediction columns)
CD,CV,HEX-CAL - ... - JM091826/	the original code and README; left untouched
Limitations
Settings were chosen on the data they are scored on: the scale, the peel filter and γ. Leave-one-image-out checks were used for the CV rules, but the sample is small (37 fields from 14 images, and fields of one image are not independent).
Nothing is tested on the CellChek D exports. The original README's figures for those frames could not be reproduced here because that data is not in this folder.
Alignment of A-file areas with markers is an assumption (marker order). It is supported by the 0.90 correlation, but two fields with matching counts (2.3 and 7.3) do not line up.
CV is still off by about 4% on larger patches and 8–12% on small ones. The cell-level fit explains roughly 90% of the variance in instrument areas, so part of the instrument's behaviour is not explained by marker positions.
The local-weight CV rule is not shown to beat the global one (18 of 37 fields closer).