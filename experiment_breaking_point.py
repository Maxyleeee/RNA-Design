import time
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from generate_structures import countS, generateS
from benchmark_tools import run_rnainverse, run_rnafold

def evaluate_structure(s):
    # Run up to 10 independent trials to distinguish 'undesignable' from 'stochastic miss'
    for _ in range(10):
        seq, ok = run_rnainverse(s, timeout_sec=5)
        if ok and seq:
            mfe, fold_ok = run_rnafold(seq)
            if fold_ok and mfe == s:
                return 1
    return 0

def test_config(L, h, theta, n_structs=100):
    cache = {}
    Wu = 1.0
    Ws = 1.0
    
    total_weight = countS(L, cache, Wu, Ws, h, theta)
    if total_weight == 0:
        return 0, 0

    structures = []
    attempts = 0
    while len(structures) < n_structs and attempts < n_structs * 3:
        attempts += 1
        s = generateS(L, cache, Wu, Ws, h, theta)
        if s and '(' in s:
            structures.append(s)

    if not structures:
        return 0, 0

    with open("breaking_point_structures_tracked.txt", "a") as f:
        for s in structures:
            f.write(f"{L},{h},{theta},{s}\n")

    success = 0
    # Process using ThreadPoolExecutor for speed
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(evaluate_structure, s) for s in structures]
        for future in as_completed(futures):
            success += future.result()

    return success, len(structures)

def main():
    with open("breaking_point_structures_tracked.txt", "w") as f:
        f.write("L,h,theta,structure\n")
        
    lengths = [50]
    n_structs = 100
    h_values = [3]
    theta = 3
    
    print(f"--- Focused Rerun (L=50, h=3, theta=3, wu=1.0, 10 trials) ---")
    print(f"{'L':<4} | {'h':<4} | {'theta':<5} | {'Success Rate':<15} | {'Score'}")
    print("-" * 55)
    
    results = {}
    
    for L in lengths:
        for h in h_values:
            success, tried = test_config(L, h, theta, n_structs)
            
            if tried == 0:
                rate_str = "Impossible"
                score = -1.0
            else:
                rate = (success / tried) * 100
                rate_str = f"{success}/{tried} ({rate:.1f}%)"
                score = rate
            
            results[(h, theta)] = score
            print(f"{L:<4} | {h:<4} | {theta:<5} | {rate_str:<15} |")
            sys.stdout.flush()

if __name__ == "__main__":
    main()
