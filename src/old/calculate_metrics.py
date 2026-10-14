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
            # First line is the count, the rest are values
            count = int(lines[0].strip())
            values = [float(x.strip()) for x in lines[1:] if x.strip()]
            # Only return up to the count specified in the first line
            return values[:count]
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return []

def main():
    root_dir = Path("SAMPLE - Hexagonal-Various-Shape/Dat_db1")
    results = []

    # Sort directories numerically
    image_dirs = sorted([d for d in root_dir.iterdir() if d.is_dir()], key=lambda x: int(x.name) if x.name.isdigit() else 999)

    for img_dir in image_dirs:
        img_id = img_dir.name
        all_areas = []
        
        # Find all A*.dat files (A1, A2, etc.)
        a_files = sorted([f for f in img_dir.glob("A*.dat")], key=lambda x: int(x.stem[1:]) if x.stem[1:].isdigit() else 0)
        
        for a_file in a_files:
            areas = read_dat_file(a_file)
            all_areas.extend(areas)
        
        if not all_areas:
            print(f"Image {img_id}: No area data found.")
            continue
            
        # Calculate metrics
        all_areas = np.array(all_areas)
        mean_area = np.mean(all_areas)
        std_area = np.std(all_areas, ddof=0)
        
        # CD = 1,000,000 / average area
        cd = 1e6 / mean_area if mean_area != 0 else np.nan
        # CV = StdDev / Mean
        cv = std_area / mean_area if mean_area != 0 else np.nan
        
        results.append({
            'Image_ID': img_id,
            'Total_Cells': len(all_areas),
            'Mean_Area': mean_area,
            'CD': cd,
            'CV': cv
        })
        print(f"Image {img_id}: Cells={len(all_areas)}, Mean={mean_area:.2f}, CD={cd:.2f}, CV={cv:.4f}")

    # Write to CSV
    output_file = "hex_metrics_summary.csv"
    if results:
        keys = results[0].keys()
        with open(output_file, 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(results)
        print(f"\nSuccessfully saved results to {output_file}")
    else:
        print("\nNo data was processed.")

if __name__ == "__main__":
    main()
