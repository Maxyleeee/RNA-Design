import os
import sys
import time
import json
import subprocess
from datetime import datetime

# Add current directory to path
sys.path.append('.')
import benchmark_tools

def load_structures(filepath):
    structs = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                structs.append((parts[0], parts[1]))
    return structs

def verify_sequence(target_struct, sequence):
    if not sequence:
        return False, None
    try:
        # RNAfold to check MFE
        cmd = ["RNAfold"]
        res = subprocess.run(cmd, input=sequence, capture_output=True, text=True, timeout=5)
        # Output format: sequence\nstructure (energy)
        lines = res.stdout.strip().split('\n')
        if len(lines) >= 2:
            mfe_struct = lines[1].split()[0]
            return mfe_struct == target_struct, mfe_struct
        return False, None
    except Exception as e:
        print(f"Error in verify_sequence: {e}")
        return False, None

def run_benchmark():
    structures_file = "benchmark_v2_structures_tracked.txt"
    results_file = "benchmark_v2_results.txt"
    json_file = "benchmark_v2_results.json"
    
    if not os.path.exists(structures_file):
        print(f"Error: {structures_file} not found.")
        return

    structs = load_structures(structures_file)
    print(f"Loaded {len(structs)} structures from {structures_file}")

    tools = {
        "Baseline": benchmark_tools.run_baseline,
        "RNAinverse": benchmark_tools.run_rnainverse,
        "NEMO": benchmark_tools.run_nemo,
        "eM2dRNAs": benchmark_tools.run_em2drnas,
        "DesiRNA": benchmark_tools.run_desirna,
        "LEARNA": benchmark_tools.run_learna,
        "NUPACK": benchmark_tools.run_nupack
    }

    all_results = []
    summary = {tool: {"success": 0, "total": 0} for tool in tools}

    timeout_sec = 60
    
    # Header for results file
    with open(results_file, 'w') as f:
        f.write(f"Benchmark v2 Results (L=50 restricted)\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"Timeout: {timeout_sec}s\n")
        f.write("-" * 80 + "\n")

    for sid, struct in structs:
        print(f"\nProcessing {sid}: {struct}")
        case_result = {"id": sid, "structure": struct, "results": {}}
        
        with open(results_file, 'a') as f:
            f.write(f"\nStructure {sid}: {struct}\n")
        
        for tool_name, tool_func in tools.items():
            print(f"  Running {tool_name}...", end="", flush=True)
            start_time = time.time()
            try:
                seq, ok = tool_func(struct, timeout_sec=timeout_sec)
                elapsed = time.time() - start_time
                
                match = False
                mfe = None
                if ok and seq:
                    match, mfe = verify_sequence(struct, seq)
                
                res_obj = {
                    "sequence": seq,
                    "ok": ok,
                    "mfe_match": match,
                    "mfe_structure": mfe,
                    "time": round(elapsed, 2)
                }
                case_result["results"][tool_name] = res_obj
                
                if match:
                    summary[tool_name]["success"] += 1
                summary[tool_name]["total"] += 1
                
                status = "SUCCESS" if match else "FAILED"
                print(f" {status} ({elapsed:.1f}s)")
                
                with open(results_file, 'a') as f:
                    f.write(f"    {tool_name:12}: {status:8} ({elapsed:.1f}s) Seq: {seq if seq else 'N/A'}\n")
                    if mfe: f.write(f"      MFE: {mfe}\n")
                    
            except Exception as e:
                print(f" ERROR: {e}")
                case_result["results"][tool_name] = {"error": str(e)}
        
        all_results.append(case_result)
        
        # Save JSON after each structure to avoid data loss
        with open(json_file, 'w') as jf:
            json.dump({"summary": summary, "details": all_results}, jf, indent=2)

    # Final summary
    print("\n" + "="*40)
    print("FINAL SUMMARY")
    print("="*40)
    with open(results_file, 'a') as f:
        f.write("\n" + "="*40 + "\n")
        f.write("FINAL SUMMARY\n")
        f.write("="*40 + "\n")
        for tool in tools:
            succ = summary[tool]["success"]
            tot = summary[tool]["total"]
            rate = (succ / tot * 100) if tot > 0 else 0
            line = f"{tool:12}: {succ}/{tot} ({rate:.1f}%)"
            print(line)
            f.write(line + "\n")

if __name__ == "__main__":
    run_benchmark()
