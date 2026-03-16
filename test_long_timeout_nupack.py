import time
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from benchmark_tools import run_nupack, run_rnafold

TOOLS = ["NUPACK"]
tool_funcs = {
    "NUPACK": run_nupack,
}

def run_single_tool(tool_name, structure, timeout_sec):
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
            return tool_name, structure, {"seq": seq, "mfe": mfe, "match": match, "time": elapsed, "status": "OK" if match else "MFE_MISMATCH"}
        else:
            return tool_name, structure, {"seq": None, "mfe": None, "match": False, "time": elapsed, "status": "FAILED"}
    except Exception as e:
        elapsed = time.time() - t0
        return tool_name, structure, {"seq": None, "mfe": None, "match": False, "time": elapsed, "status": f"ERROR: {str(e)}"}

def main():
    L50_file = "breaking_point_structures_tracked.txt"
    L100_file = "structures_motif_h4.txt"
    timeout = 120 # 2 minutes per tool per structure
    num_samples = 5
    
    # Read L=50 structures
    l50_structs = []
    if os.path.exists(L50_file):
        with open(L50_file, "r") as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) == 4 and parts[0] == '50':
                    l50_structs.append(parts[3])
                    if len(l50_structs) >= num_samples: break

    # Read L=100 structures
    l100_structs = []
    if os.path.exists(L100_file):
        with open(L100_file, "r") as f:
            for line in f:
                s = line.strip()
                if set(s).issubset({'(', ')', '.'}):
                    l100_structs.append(s)
                    if len(l100_structs) >= num_samples: break
                    
    structures = [("L=50", s) for s in l50_structs] + [("L=100", s) for s in l100_structs]
    print(f"Testing {len(structures)} structures with a timeout of {timeout}s.")
    
    tools_to_test = TOOLS
    results = { (tool, struct): None for tool in tools_to_test for _, struct in structures }
    
    # We will submit all jobs to a thread pool
    futures_map = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        for length_label, struct in structures:
            for tool in tools_to_test:
                future = executor.submit(run_single_tool, tool, struct, timeout)
                futures_map[future] = (tool, length_label, struct)
                
        completed = 0
        total = len(futures_map)
        
        for future in as_completed(futures_map):
            completed += 1
            tool, struct, res = future.result()
            results[(tool, struct)] = res
            print(f"[{completed}/{total}] {tool} on {len(struct)}-nt struct -> {res['status']} ({res['time']:.1f}s)")
            sys.stdout.flush()
            
    # Print and save summary
    with open("long_timeout_results_nupack.txt", "w") as f:
        print("\n=== SUMMARY ===")
        f.write("=== SUMMARY ===\n")
        
        for length_label in ["L=50", "L=100"]:
            structs = [s for label, s in structures if label == length_label]
            if not structs: continue
            
            print(f"\nResults for {length_label} (N={len(structs)}):")
            f.write(f"\nResults for {length_label} (N={len(structs)}):\n")
            
            # Print table
            header = f"{'Tool':<15} {'Success Rate':<15} {'Avg Time (s)':<15}"
            print(header)
            f.write(header + "\n")
            print("-" * 50)
            f.write("-" * 50 + "\n")
            
            for tool in tools_to_test:
                successes = 0
                total_time = 0
                for struct in structs:
                    res = results[(tool, struct)]
                    if res["match"]: successes += 1
                    total_time += res["time"]
                
                avg_time = total_time / len(structs) if len(structs) > 0 else 0
                row = f"{tool:<15} {successes}/{len(structs):<13} {avg_time:<15.1f}"
                print(row)
                f.write(row + "\n")
                
                
        # Write details
        f.write("\n\n=== DETAILED RESULTS ===\n")
        for length_label, struct in structures:
            f.write(f"\nStructure ({length_label}): {struct}\n")
            for tool in tools_to_test:
                res = results[(tool, struct)]
                if res['match']:
                    f.write(f"  {tool:12s}: SUCCESS in {res['time']:.1f}s (seq={res['seq']})\n")
                else:
                    f.write(f"  {tool:12s}: {res['status']} in {res['time']:.1f}s\n")

if __name__ == '__main__':
    main()
