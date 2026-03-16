import time
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from generate_structures import countS, generateS
from benchmark_tools import run_rnainverse, run_rnafold

def evaluate_structure(s):
    seq, ok = run_rnainverse(s, timeout_sec=5)
    if ok and seq:
        mfe, fold_ok = run_rnafold(seq)
        if fold_ok and mfe == s:
            return 1
    return 0

def test_config(L, h, theta, wu, ws, n_structs=100):
    cache = {}
    total_weight = countS(L, cache, wu, ws, h, theta)
    if total_weight == 0:
        return 0, 0

    structures = []
    attempts = 0
    # Try slightly harder to find valid structures if weights make them rare
    while len(structures) < n_structs and attempts < n_structs * 5:
        attempts += 1
        s = generateS(L, cache, wu, ws, h, theta)
        if s and '(' in s:
            structures.append(s)

    if not structures:
        return 0, 0

    success = 0
    # Process using ThreadPoolExecutor for speed
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(evaluate_structure, s) for s in structures]
        for future in as_completed(futures):
            success += future.result()

    return success, len(structures)

def experiment_A():
    L = 100
    n_structs = 100
    h_values = [4, 5]
    theta_values = [3, 4, 5, 6, 7, 8]
    wu, ws = 1.0, 1.0
    
    print(f"--- Experiment A: Impact of theta (L={L}, Wu={wu}, Ws={ws}) ---")
    print(f"{'h':<4} | {'theta':<5} | {'Success Rate':<15}")
    print("-" * 35)
    
    for h in h_values:
        for theta in theta_values:
            success, tried = test_config(L, h, theta, wu, ws, n_structs)
            if tried == 0:
                rate_str = "Impossible"
            else:
                rate = (success / tried) * 100
                rate_str = f"{success}/{tried} ({rate:.1f}%)"
            
            print(f"{h:<4} | {theta:<5} | {rate_str:<15}")
            sys.stdout.flush()

def experiment_B():
    L = 100
    n_structs = 100
    h_values = [4, 5]
    theta = 3
    wu_values = [0.2, 0.5, 1.0, 2.0, 5.0]
    ws = 1.0
    
    print(f"\n--- Experiment B: Impact of Unpaired Weight (L={L}, theta={theta}, Ws={ws}) ---")
    print(f"{'h':<4} | {'Wu':<5} | {'Success Rate':<15}")
    print("-" * 35)
    
    for h in h_values:
        for wu in wu_values:
            success, tried = test_config(L, h, theta, wu, ws, n_structs)
            if tried == 0:
                rate_str = "Impossible"
            else:
                rate = (success / tried) * 100
                rate_str = f"{success}/{tried} ({rate:.1f}%)"
            
            print(f"{h:<4} | {wu:<5.1f} | {rate_str:<15}")
            sys.stdout.flush()

def main():
    experiment_A()
    experiment_B()

if __name__ == "__main__":
    main()
