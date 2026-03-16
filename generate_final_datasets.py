import sys
import os

# Append current dir to path to import from generate_structures
sys.path.append('.')
from generate_structures import countS, generateS

configs = [
    {"L": 50, "h": 2, "theta": 3, "Wu": 1.0, "Ws": 1.0, "name": "l50_h2_t3_wu1"},
    {"L": 50, "h": 3, "theta": 3, "Wu": 1.0, "Ws": 1.0, "name": "l50_h3_t3_wu1"},
    {"L": 50, "h": 4, "theta": 3, "Wu": 1.0, "Ws": 1.0, "name": "l50_h4_t3_wu1"},
    {"L": 50, "h": 3, "theta": 3, "Wu": 0.2, "Ws": 1.0, "name": "l50_h3_t3_wu0.2"},
    {"L": 50, "h": 3, "theta": 3, "Wu": 5.0, "Ws": 1.0, "name": "l50_h3_t3_wu5"},
    {"L": 100, "h": 2, "theta": 3, "Wu": 1.0, "Ws": 1.0, "name": "l100_h2_t3_wu1"},
    {"L": 100, "h": 3, "theta": 3, "Wu": 1.0, "Ws": 1.0, "name": "l100_h3_t3_wu1"},
]

N = 100
out_dir = "output/final_datasets"
os.makedirs(out_dir, exist_ok=True)

for cfg in configs:
    L, h, theta, Wu, Ws, name = cfg["L"], cfg["h"], cfg["theta"], cfg["Wu"], cfg["Ws"], cfg["name"]
    cache = {}
    print(f"Generating {N} structures for {name}...")
    
    # Pre-calculate probabilities for generation
    total_weight = countS(L, cache, Wu, Ws, h, theta)
    if total_weight == 0:
        print(f"  Warning: Total weight is 0. Cannot generate structures for this configuration.")
        continue
        
    structures = []
    failed = 0
    for i in range(N):
        s = generateS(L, cache, Wu, Ws, h, theta)
        if s:
            structures.append(s)
        else:
            failed += 1
            
    out_file = f"{out_dir}/{name}.txt"
    with open(out_file, "w") as f:
        for s in structures:
            f.write(s + "\n")
    print(f"  -> Saved {len(structures)} to {out_file} (Failed: {failed})")

print("\nAll datasets generated successfully!")
