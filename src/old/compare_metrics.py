import os
import numpy as np
import csv
from pathlib import Path

def read_dat_file(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.read().splitlines()
            if not lines:
                return []
            count = int(lines[0].strip())
            values = [float(x.strip()) for x in lines[1:] if x.strip()]
            return values[:count]
    except Exception:
        return []

def calculate_stats(data):
    if not data:
        return np.nan, np.nan, 0
    arr = np.array(data)
    mean = np.mean(arr)
    std = np.std(arr, ddof=0)
    cv = std / mean if mean != 0 else np.nan
    return mean, cv, len(arr)

def main():
    root_dir = Path("SAMPLE - Hexagonal-Various-Shape/Dat_db1")
    comparison_results = []

    image_dirs = sorted([d for d in root_dir.iterdir() if d.is_dir()], key=lambda x: int(x.name) if x.name.isdigit() else 999)

    for img_dir in image_dirs:
        img_id = img_dir.name
        res = {'Image_ID': img_id}
        
        for p in ['A', 'C', 'P']:
            all_vals = []
            files = sorted([f for f in img_dir.glob(f"{p}*.dat")], key=lambda x: int(x.stem[1:]) if x.stem[1:].isdigit() else 0)
            for f in files:
                all_vals.extend(read_dat_file(f))
            
            m, cv, c = calculate_stats(all_vals)
            res[f"{p}_Mean"] = m
            res[f"{p}_CV"] = cv
            res[f"{p}_Count"] = c
        
        comparison_results.append(res)

    output_file = "comparison_metrics_summary.csv"
    if comparison_results:
        keys = comparison_results[0].keys()
        with open(output_file, 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(comparison_results)
        print(f"Successfully saved comparison to {output_file}")

if __name__ == "__main__":
    main()
