import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from generate_structures import generateS, countS
from benchmark_tools import benchmark_file

def generate_dataset(l, h, theta, wu, num_structures, output_filename):
    print(f"\nGenerating dataset: l={l}, h={h}, theta={theta}, wu={wu}, N={num_structures}")
    ws = 1.0
    cache = {}
    total_weight = countS(l, cache, wu, ws, h, theta)
    print(f"Total weight: {total_weight:.2e}")
    
    structures = []
    for i in range(num_structures):
        s = generateS(l, cache, wu, ws, h, theta)
        if s:
            structures.append(s)
            
    with open(output_filename, 'w') as f:
        for s in structures:
            f.write(s + "\n")
            
    print(f"Saved {len(structures)} structures to {output_filename}")
    return output_filename

def main():
    os.makedirs("output", exist_ok=True)
    num_structures = 10
    
    # Dataset 1: l=50, h=2, theta=3, wu=1
    f1 = generate_dataset(l=50, h=2, theta=3, wu=1.0, num_structures=num_structures, 
                          output_filename="output/dataset_l50_h2_t3_wu1.txt")
    
    # Dataset 2: l=50, h=3, theta=3, wu=1
    f2 = generate_dataset(l=50, h=3, theta=3, wu=1.0, num_structures=num_structures, 
                          output_filename="output/dataset_l50_h3_t3_wu1.txt")
                          
    # Benchmark both
    benchmark_file(f1, max_structures=num_structures)
    benchmark_file(f2, max_structures=num_structures)

if __name__ == "__main__":
    main()
