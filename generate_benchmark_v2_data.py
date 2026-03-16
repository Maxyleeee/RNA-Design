import sys
import os
import random

# Add current directory to path to import generate_with_motif
sys.path.append('.')
from generate_with_motif import generateS, countS

def generate_set(l, h, theta, wu, count):
    ws = 1.0  # weight of stack
    wm = 0.0  # weight of motif (baseline)
    
    cache = {}
    total_weight = countS(l, cache, wu, ws, wm, h, theta)
    print(f"  Total weight for h={h}, theta={theta}: {total_weight:.2e}")
    
    structures = []
    attempts = 0
    while len(structures) < count and attempts < count * 100:
        struct = generateS(l, cache, wu, ws, wm, h, theta)
        if struct and struct not in structures:
            structures.append(struct)
        attempts += 1
    return structures

if __name__ == "__main__":
    l = 50
    wu = 1
    theta = 3
    
    # h = 2
    print(f"Generating 50 structures for l={l}, h=2, theta={theta}, wu={wu}...")
    h2_structs = generate_set(l, 2, theta, wu, 50)
    
    # h = 3
    print(f"Generating 50 structures for l={l}, h=3, theta={theta}, wu={wu}...")
    h3_structs = generate_set(l, 3, theta, wu, 50)
    
    output_file = "benchmark_v2_structures_tracked.txt"
    with open(output_file, 'w') as f:
        f.write("# L=50, h=2, theta=3, wu=1 (50 structures)\n")
        f.write("# ID Format: h2_0..49\n")
        for i, s in enumerate(h2_structs):
            f.write(f"h2_{i} {s}\n")
        
        f.write("\n")
        f.write("# L=50, h=3, theta=3, wu=1 (50 structures)\n")
        f.write("# ID Format: h3_0..49\n")
        for i, s in enumerate(h3_structs):
            f.write(f"h3_{i} {s}\n")
            
    print(f"Successfully generated 100 structures and saved to {output_file}")
