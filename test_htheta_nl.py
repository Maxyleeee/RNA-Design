import time
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from benchmark_tools import run_nupack, run_learna, run_rnafold

TOOLS = ["LEARNA", "NUPACK"]
tool_funcs = {
    "LEARNA": run_learna,
    "NUPACK": run_nupack,
}

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
    structs_file = "custom_htheta_structures.txt"
    if not os.path.exists(structs_file):
        print(f"File {structs_file} not found. Run the generation script first.")
        return

    all_structs = []
    with open(structs_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(" | ")
            if len(parts) == 2:
                all_structs.append((parts[0].strip(), parts[1].strip()))

    print(f"Loaded {len(all_structs)} structures from {structs_file}.")
    
    timeout = 120 # 2 minutes per tool
    
    print(f"Testing tools: {', '.join(TOOLS)}")
    print(f"Max Timeout: {timeout}s per structure per tool.\n")
    
    results = { (tool, label, struct): None for tool in TOOLS for label, struct in all_structs }
    
    futures_map = {}
    # Decrease workers slightly to avoid overloading with NUPACK/LEARNA
    with ThreadPoolExecutor(max_workers=4) as executor:
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
            
    # Summarize and save results
    results_file = "custom_htheta_results_nl.txt"
    labels_found = list(dict.fromkeys([l for l, _ in all_structs]))
    
    with open(results_file, "w") as f:
        print("\n=== SUMMARY ===")
        f.write("=== SUMMARY ===\n")
        
        for label in labels_found:
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
