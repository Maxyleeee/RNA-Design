import time
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import generation functions
from generate_with_motif import countS, generateS

# Import benchmark functions
from benchmark_tools import run_rnainverse, run_nemo, run_em2drnas, run_desirna, run_rnafold

# We exclude LEARNA and NUPACK as requested
TOOLS = ["RNAinverse", "NEMO", "eM2dRNAs", "DesiRNA"]
tool_funcs = {
    "RNAinverse": run_rnainverse,
    "NEMO": run_nemo,
    "eM2dRNAs": run_em2drnas,
    "DesiRNA": run_desirna,
}

def generate_structures_for_config(L, wu, ws, wm, h, theta, num_structures):
    cache = {}
    structures = []
    
    # Precompute weights to populate cache
    total_weight = countS(L, cache, wu, ws, wm, h, theta)
    if total_weight <= 0:
        print(f"Warning: Configuration L={L}, h={h}, theta={theta} yields 0 possible structures.")
        return structures
        
    print(f"Generating {num_structures} structures for L={L}, h={h}, theta={theta} (Total weight: {total_weight:.2e})")
    for i in range(num_structures):
        s = generateS(L, cache, wu, ws, wm, h, theta)
        if s:
            structures.append(s)
            
    return structures

def run_single_tool(tool_name, config_label, structure, timeout_sec):
    t0 = time.time()
    try:
        import inspect
        sig = inspect.signature(tool_funcs[tool_name])
        if 'timeout_sec' in sig.parameters:
            seq, ok = tool_funcs[tool_name](structure, timeout_sec=timeout_sec)
        else:
            seq, ok = tool_funcs[tool_name](structure)
            
        elapsed = time.time() - t0
        
        if ok and seq:
            mfe, fold_ok = run_rnafold(seq)
            match = fold_ok and (mfe == structure)
            status = "SUCCESS" if match else "MFE_MISMATCH"
            return tool_name, config_label, structure, {"seq": seq, "mfe": mfe, "match": match, "time": elapsed, "status": status}
        else:
            return tool_name, config_label, structure, {"seq": None, "mfe": None, "match": False, "time": elapsed, "status": "FAILED"}
    except Exception as e:
        elapsed = time.time() - t0
        return tool_name, config_label, structure, {"seq": None, "mfe": None, "match": False, "time": elapsed, "status": f"ERROR: {str(e)}"}

def main():
    num_samples = 25 # number of structures per configuration
    timeout = 120 # 2 minutes per tool
    
    # 1. Generate structures
    config_1 = {"L": 50, "wu": 1.0, "ws": 1.0, "wm": 0.0, "h": 2, "theta": 3, "label": "L=50, h=2, theta=3"}
    config_2 = {"L": 50, "wu": 1.0, "ws": 1.0, "wm": 0.0, "h": 3, "theta": 3, "label": "L=50, h=3, theta=3"}
    
    structs_1 = generate_structures_for_config(**{k:v for k,v in config_1.items() if k != "label"}, num_structures=num_samples)
    structs_2 = generate_structures_for_config(**{k:v for k,v in config_2.items() if k != "label"}, num_structures=num_samples)
    
    all_structs = [(config_1["label"], s) for s in structs_1] + [(config_2["label"], s) for s in structs_2]
    
    # 2. Save structures to a tracking file
    structs_file = "custom_htheta_structures.txt"
    with open(structs_file, "w") as f:
        f.write("# Format: label, structure\n")
        for label, s in all_structs:
            f.write(f"{label} | {s}\n")
    print(f"\nSaved {len(all_structs)} generated structures to {structs_file}\n")
    
    # 3. Benchmark tools
    print(f"Testing tools: {', '.join(TOOLS)}")
    print(f"Max Timeout: {timeout}s per structure per tool.\n")
    
    results = { (tool, label, struct): None for tool in TOOLS for label, struct in all_structs }
    
    futures_map = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        for label, struct in all_structs:
            for tool in TOOLS:
                future = executor.submit(run_single_tool, tool, label, struct, timeout)
                futures_map[future] = (tool, label, struct)
                
        completed = 0
        total = len(futures_map)
        
        for future in as_completed(futures_map):
            completed += 1
            tool, label, struct, res = future.result()
            results[(tool, label, struct)] = res
            print(f"[{completed}/{total}] {tool} on {label} -> {res['status']} ({res['time']:.1f}s)")
            sys.stdout.flush()
            
    # 4. Summarize and save results
    results_file = "custom_htheta_results.txt"
    with open(results_file, "w") as f:
        print("\n=== SUMMARY ===")
        f.write("=== SUMMARY ===\n")
        
        for label in [config_1["label"], config_2["label"]]:
            structs = [s for l, s in all_structs if l == label]
            if not structs: continue
            
            print(f"\nResults for {label} (N={len(structs)}):")
            f.write(f"\nResults for {label} (N={len(structs)}):\n")
            
            header = f"{'Tool':<15} {'Success Rate':<15} {'Avg Time (s)':<15}"
            print(header)
            f.write(header + "\n")
            print("-" * 50)
            f.write("-" * 50 + "\n")
            
            for tool in TOOLS:
                successes = 0
                total_time = 0
                for struct in structs:
                    res = results[(tool, label, struct)]
                    if res["match"]: successes += 1
                    total_time += res["time"]
                
                avg_time = total_time / len(structs) if len(structs) > 0 else 0
                row = f"{tool:<15} {successes}/{len(structs):<13} {avg_time:<15.1f}"
                print(row)
                f.write(row + "\n")
                
        f.write("\n\n=== DETAILED RESULTS ===\n")
        for label, struct in all_structs:
            f.write(f"\nStructure ({label}): {struct}\n")
            for tool in TOOLS:
                res = results[(tool, label, struct)]
                if res['match']:
                    f.write(f"  {tool:12s}: SUCCESS in {res['time']:.1f}s (seq={res['seq']})\n")
                else:
                    f.write(f"  {tool:12s}: {res['status']} in {res['time']:.1f}s\n")
                    
    print(f"\nDetailed results saved to {results_file}")

if __name__ == '__main__':
    main()
