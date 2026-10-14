Code Walkthrough: Reproducing CellChek D Metrics from Marker Positions
Script-by-script explanation · CellChek D reverse-engineering project · 8 October 2026
Overview
The goal of the project is to reproduce the numbers a CellChek D specular microscope prints (Num, CD, CV and HEX) from the marked cell centres alone. This document goes through the code in the order you would present it: the core scoring script first, then the diagnostic scripts, then the original and older code.
1. score_deks.py (the core)
Purpose: run the model on every field and compare it with the instrument's answer.
Settings at the top. The scale (0.7716 µm/px by default, 0.768 when you pass --um), the peel-filter values, the minimum wall length, and the number of neighbours used to cut each cell (16).
Reading files.
• read_A takes the first number as the cell count, then the areas, then the side counts if a block of whole numbers between 3 and 12 follows.
• read_C takes the first number as the marker count, then the x,y pairs.
Geometry.
• power_cells is the heart of it. For each marker it starts with a huge square and cuts it with a straight line halfway toward each of its 16 nearest neighbours. What is left is that marker's Voronoi cell. If weights are given, each cut line is shifted toward the smaller cell.
• area is the shoelace formula, the same one used in the Li et al. (2025) paper.
• interior keeps a cell only if all its corners lie inside the convex hull of all markers, which removes edge cells.
• keep_cells applies that test, plus the optional peel filter (now switched off) and the optional --alpha concave boundary.
The model, analyse.
1. Build plain Voronoi cells and keep the interior ones.
2. Num is the count, mean area gives CD (1e6 / round(mean)), and HEX is the share of six-sided cells.
3. For CV, compute a weight for each cell from how big it is, rebuild the diagram with those weights, and take standard deviation divided by mean.
Scoring, main. For each field it builds ground truth from the A file (count, 1e6/mean, standard deviation/mean, share of 6-sided cells), runs analyse, writes scores.csv, and prints bias and MAPE.
Output: scores.csv plus the printed table.
2. hex_sweep.py
Question: does merging short walls fix HEX?
It imports score_deks as a module, tries minimum wall lengths of 0, 0.5, 1, 2 and 3 px, and recomputes HEX each time. It also prints the instrument's share of 4 to 8-sided cells beside ours.
Finding: merging makes HEX worse, and the side distributions already agree within about 1 point.
3. boundry_check.py
Question: does the Num error depend on patch size or boundary?
It reads scores.csv, reopens each C file, and computes the number of hull corners, the hull perimeter in units of the typical cell spacing, and the spacing itself. It correlates each with the Num error, then fits cells removed against perimeter for both the instrument and our model.
Finding: the error follows patch size (correlation about 0.5), not the number of hull corners.
4. percell_check.py
Question: are the A-file areas listed in marker order, and how do they differ from Voronoi?
It uses only fields where our kept count equals the instrument's. For each one it compares our area with the instrument's area cell by cell: correlation, slope, side agreement and both CVs.
Finding: the correlation is about 0.9, so the order matches. The instrument's areas are about twice as spread out as ours.
5. stretch_check.py
Question: what rule turns Voronoi areas into the instrument's areas?
On the aligned fields it builds three features per cell: its size relative to the frame, its size relative to its neighbours, and its side count. It regresses the instrument's log area on them and reports R². It also scores each rule's CV, always leaving one whole image out.
Finding: the neighbour comparison and the side count each add about 0.1 to R².
6. local_model_check.py
Question: does the neighbour-based rule predict CV?
The rule is area = V × (V / mean of neighbours)^γ. It fits γ on the aligned fields excluding one image, then scores that held-out image on all of its fields.
Finding: γ is about 1.36, and CV error is 4.8% against 5.5% for the old weights.
7. The original files
• metrics.py: the original model (the same recipe as analyse), but it needs the missing ecd package.
• alpha_roi.py: a concave-boundary builder using Delaunay triangles. Nothing calls it, and the --alpha option in score_deks.py does the same job more robustly.
• run_all.py: the original scoring script. It needs read_dat and deks_dat, which are not in the folder, and it calls an undefined rows_clear.
8. The older scripts
• calculate_metrics.py, full_analysis.py, calculate_cv_metrics.py: read A files only, pool all fields of an image, and compute mean, CD and CV. The last one also adds an invented XV proxy.
• compare_metrics.py: reads A, C and P files with one reader that suits only A files, so its C and P columns are wrong.
• generate_scores.py: copies ground truth into the prediction columns, so its output proves nothing.
Settings in force and caveats
Scale 0.768 µm/px, peel filter off (MIN_PEELED = 99), MIN_EDGE_PX = 0, no alpha boundary.
• The scale, the peel filter and γ were chosen using the same DEKS data they are scored on, so the numbers are not fully held out.
• The local-contrast CV rule lives in local_model_check.py, not in score_deks.py.
• Nothing has been tested on the CellChek D exports, which are not in the project folder.