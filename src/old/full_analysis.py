import os
import numpy as np
import csv
from pathlib import Path

def read_dat_file(file_path):
    """Reads a .dat file where the first line is the count."""
    try:
        with open(file_path, 'r') as f:
            lines = f.read().splitlines()
            if not lines: return []
            count_line = lines[0].strip()
            if not count_line: return []
            count = int(count_line)
            values = [float(x.strip()) for x in lines[1:] if x.strip()]
            return values[:count]
    except Exception:
        return []

def calculate_advanced_metrics(areas):
    if not areas:
        return np.nan, np.nan, 0, np.nan
    
    arr = np.array(areas)
    mean_area = np.mean(arr)
    std_area = np.std(arr, ddof=0)
    
    # Standard Metrics
    ecd = 1e6 / mean_area if mean_area != 0 else np.nan
    cv = std_area / mean_area if mean_area != 0 else np.nan
    
    # XV Proxy (Morphological Deviation)
    # Normalized based on the research paper's link between CV and stress
    xv_proxy = cv * (mean_area / 300.0) 
    
    return mean_area, ecd, cv, xv_proxy, len(arr)

def analyze_morphology(areas):
    """
    Simulates logic from stretch_check.py and hex_sweep.py 
    based on area distribution.
    """
    if not areas:
        return 0.0, 0.0
    
    arr = np.array(areas)
    # In a real scenario with vertices, we would count sides.
    # Here we use the CV as a proxy for morphological distortion (Stretch).
    stretch_index = np.std(arr) / np.mean(arr) if np.mean(arr) != 0 else 0
    
    # Hexagon percentage proxy: In healthy cells, CV is low.
    # Higher CV usually means fewer perfect hexagons.
    hex_percentage = max(0, 1.0 - (stretch_index * 2)) 
    
    return stretch_index, hex_percentage

def process_dataset(dataset_path, dataset_name):
    print(f"\nProcessing Dataset: {dataset_name}")
    print("-" * 70)
    
    data_root = None
    # Handle different folder naming conventions (Dat_db1 or DAT_DEKS_...)
    for p in dataset_path.iterdir():
        if p.is_dir() and (p.name.startswith("Dat") or p.name.startswith("DAT")):
            data_root = p
            break
    
    if not data_root:
        print(f"Could not find a data directory in {dataset_path}")
        return []

    results = []
    img_dirs = sorted([d for d in data_root.iterdir() if d.is_dir()], 
                     key=lambda x: int(x.name) if x.name.isdigit() else 999)

    for img_dir in img_dirs:
        img_id = img_dir.name
        all_areas = []
        
        # Find all A*.dat files (Actual areas)
        a_files = sorted([f for f in img_dir.glob("A*.dat")], 
                        key=lambda x: int(x.stem[1:]) if x.stem[1:].isdigit() else 0)
        
        for a_file in a_files:
            all_areas.extend(read_dat_file(a_file))
        
        if not all_areas:
            continue
            
        # 1. Standard and Functional Metrics
        mean, ecd, cv, xv, count = calculate_advanced_metrics(all_areas)
        
        # 2. Morphology Checks (Integrating stretch_check and hex_sweep logic)
        stretch, hex_perc = analyze_morphology(all_areas)
        
        # 3. Boundary Check Proxy
        # In the absence of vertex data, we estimate 'border' loss 
        # based on the ratio of cells found vs expected blocks.
        expected_cells = len(a_files) * 30 # Assumption: ~30 cells per block
        boundary_loss_est = max(0, (expected_cells - count) / expected_cells) if expected_cells > 0 else 0

        results.append({
            'Dataset': dataset_name,
            'Image_ID': img_id,
            'Cell_Count': count,
            'Mean_Area': mean,
            'ECD': ecd,
            'CV': cv,
            'XV_Proxy': xv,
            'Stretch_Index': stretch,
            'Hex_Percentage': hex_perc,
            'Boundary_Loss_Est': boundary_loss_est
        })
        print(f"Image {img_id}: ECD={ecd:.2f}, CV={cv:.4f}, XV={xv:.4f}, Hex={hex_perc:.2%}")
    
    return results

def main():
    # Updated paths to cover all sample folders
    datasets = [
        (Path("SAMPLE - DEKS-EYEBANK/Reader1"), "DEKS-Eyebank"),
        (Path("SAMPLE - Hexagonal-Various-Shape"), "Hexagonal-Various")
    ]
    
    all_final_results = []
    for path, name in datasets:
        if path.exists():
            all_final_results.extend(process_dataset(path, name))
        else:
            print(f"Path not found: {path}")

    if all_final_results:
        output_file = "final_combined_analysis.csv"
        keys = all_final_results[0].keys()
        with open(output_file, 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(all_final_results)
        print(f"\n\nSUCCESS: Integrated analysis complete. Results saved to {output_file}")
    else:
        print("\nNo data found to process.")

if __name__ == "__main__":
    main()
