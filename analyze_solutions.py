import re
import sys

def analyze_benchmark(filepath):
    print(f"Analyzing {filepath}...")
    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # A structure entry starts with "Structure:"
    # And successful designs usually have "Match: YES"
    
    structures = content.split("Target: ")
    # First split is before any structure
    if len(structures) <= 1:
        print("No structures found (delimiter mismatch).")
        return

    total_structures = len(structures) - 1
    any_success = 0
    
    for entry in structures[1:]:
        # Check if "Match: YES" exists anywhere in the tools list for this structure
        if "Match: YES" in entry:
            any_success += 1
            
    print(f"Total Structures: {total_structures}")
    print(f"Structures with at least one solution: {any_success}")
    if total_structures > 0:
        print(f"Success Rate: {(any_success/total_structures)*100:.1f}%")
    print("-" * 30)

if __name__ == "__main__":
    files = [
        "/home/maxyle/RNA/output/dataset_l50_h4_t3_wu1_2motifs_benchmark.txt",
        "/home/maxyle/RNA/output/dataset_l50_h4_t3_wu1_2motifs2_benchmark.txt"
    ]
    for f in files:
        analyze_benchmark(f)
