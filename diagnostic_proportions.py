import sys
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from benchmark_tools import run_rnainverse, run_rnafold

def get_stats_for_structure(s):
    results = [] # List of (design_ok, match_ok)
    for _ in range(10):
        seq, ok = run_rnainverse(s, timeout_sec=5)
        if ok and seq:
            mfe, fold_ok = run_rnafold(seq)
            match = (fold_ok and mfe == s)
            results.append((True, match))
        else:
            results.append((False, False))
    return results

def main():
    log_file = "breaking_point_structures_tracked.txt"
    target_L, target_h = 50, 3
    structures = []
    
    with open(log_file, "r") as f:
        f.readline() # header
        for line in f:
            parts = line.strip().split(",")
            if len(parts) >= 4 and int(parts[0]) == target_L and int(parts[1]) == target_h:
                structures.append(parts[3])

    structures = structures[:20] # Limit to 20 structures for speed (200 trials total)
    print(f"Analyzing {len(structures)} structures for proportions...")

    total_designs = 0
    total_matches = 0

    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(get_stats_for_structure, s) for s in structures]
        for future in as_completed(futures):
            res_list = future.result()
            for design_ok, match in res_list:
                if design_ok:
                    total_designs += 1
                    if match:
                        total_matches += 1

    if total_designs > 0:
        fail_rate = (1 - (total_matches / total_designs)) * 100
        print(f"\nStats for L=50, h=3:")
        print(f"Total Designs attempted by RNAinverse: {total_designs}")
        print(f"Verified by RNAfold: {total_matches}")
        print(f"Proportion that FAILED verification: {fail_rate:.1f}%")
    else:
        print("No designs found.")

if __name__ == "__main__":
    main()
