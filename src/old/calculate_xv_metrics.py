import os
import numpy as np
import csv
from pathlib import Path

def read_dat_file(file_path):
    """Reads a .dat file where the first line is the count."""
    try:
        with open(file_path, 'r') as f:
            lines = f.read().splitlines()
            if not lines:
                return []
            count_line = lines[0].strip()
            if not count_line:
                return []
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
    
    # XV Approximation logic:
    # In the absence of the specific XV formula from the paper, 
    # we use the 'Morphological Deviation' proxy.
    # XV is linked to shape distortion. A common proxy for distortion 
    # in polygons is the ratio of the actual area to the area of a 
    # regular polygon (hexagon) with the same perimeter.
    # Since we only have Areas here, we use the CV as the primary 
    # indicator of morphological heterogeneity (which the paper links to XV).
    xv_proxy = cv * (mean_area / 300.0) # Normalized to a typical healthy mean area
    
    return mean_area, ecd, cv, xv_proxy, len(arr)

def process_dataset(dataset_path, dataset_name):
    print(f"\nProcessing Dataset: {dataset_name}")
    print("-" * 60)
    
    data_root = None
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
        
        a_files = sorted([f for f in img_dir.glob("A*.dat")], 
                        key=lambda x: int(x.stem[1:]) if x.stem[1:].isdigit() else 0)
        
        for a_file in a_files:
            all_areas.extend(read_dat_file(a_file))
        
        if not all_areas:
            continue
            
        mean, ecd, cv, xv, count = calculate_advanced_metrics(all_areas)
        
        results.append({
            'Dataset': dataset_name,
            'Image_ID': img_id,
            'Cell_Count': count,
            'Mean_Area': mean,
            'ECD': ecd,
            'CV': cv,
            'XV_Proxy': xv
        })
        print(f"Image {img_id}: Count={count}, ECD={ecd:.2f}, CV={cv:.4f}, XV_Proxy={xv:.4f}")
    
    return results

def main():
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
        output_file = "xv_functional_analysis.csv"
        keys = all_final_results[0].keys()
        with open(output_file, 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(all_final_results)
        print(f"\n\nSUCCESS: Functional analysis complete. Results saved to {output_file}")
    else:
        print("\nNo data found to process.")

if __name__ == "__main__":
    main()
