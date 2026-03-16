import os
import sys
from generate_with_motif import countS, generateS
from benchmark_tools import benchmark_file

def generate_motif_file(filename, n_structs, L, wu, ws, wm, h, theta):
    cache = {}
    total_weight = countS(L, cache, wu, ws, wm, h, theta)
    if total_weight == 0:
        print(f"Warning: Zero weight for configuration L={L}, h={h}")
        return False
        
    structures = []
    attempts = 0
    while len(structures) < n_structs and attempts < n_structs * 10:
        attempts += 1
        s = generateS(L, cache, wu, ws, wm, h, theta)
        if s and '(' in s:
            structures.append(s)
            
    if structures:
        with open(filename, 'w') as f:
            for s in structures:
                f.write(s + "\n")
        print(f"Generated {len(structures)} structures -> {filename}")
        return True
    return False

def main():
    L = 100
    wu = 1.0
    ws = 1.0
    wm = 5.0 # High weight to ensure motifs
    theta = 3
    n_structs = 25 # Testing 25 structures per configuration across all tools
    
    # Configuration 1: Borderline breaking point (h=4)
    file_h4 = "structures_motif_h4.txt"
    if generate_motif_file(file_h4, n_structs, L, wu, ws, wm, 4, theta):
        benchmark_file(file_h4, max_structures=n_structs)
        
    # Configuration 2: Transition zone (h=5)
    file_h5 = "structures_motif_h5.txt"
    if generate_motif_file(file_h5, n_structs, L, wu, ws, wm, 5, theta):
        benchmark_file(file_h5, max_structures=n_structs)

if __name__ == "__main__":
    main()
