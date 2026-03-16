import sys
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from experiment_breaking_point import evaluate_structure

def main():
    log_file = "breaking_point_structures_tracked.txt"
    if not os.path.exists(log_file):
        print(f"Error: {log_file} not found")
        return

    target_L = 50
    target_h = 3
    structures = []

    with open(log_file, "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 4: continue
            L, h, theta, s = parts
            if int(L) == target_L and int(h) == target_h:
                structures.append(s)

    if not structures:
        print(f"No structures found for L={target_L}, h={target_h}")
        return

    print(f"--- Rerunning {len(structures)} structures for L={target_L}, h={target_h} (10 Trials each) ---")
    
    success = 0
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(evaluate_structure, s) for s in structures]
        for i, future in enumerate(as_completed(futures)):
            res = future.result()
            success += res
            print(f"Progress: {i+1}/{len(structures)} | Success: {success}", end='\r')
            sys.stdout.flush()

    rate = (success / len(structures)) * 100
    print(f"\nFinal Success Rate for L=50, h=3: {success}/{len(structures)} ({rate:.1f}%)")

if __name__ == "__main__":
    main()
