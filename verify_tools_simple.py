import benchmark_tools
import time

test_structure = "((((....))))"

def verify_tool(name, func):
    print(f"--- Verifying {name} ---")
    start = time.time()
    seq, ok = func(test_structure)
    elapsed = time.time() - start
    
    if ok and seq:
        print(f"  [OK] Sequence: {seq} ({elapsed:.1f}s)")
        mfe, fold_ok = benchmark_tools.run_rnafold(seq)
        if fold_ok and mfe:
            match = (mfe == test_structure)
            print(f"  [OK] MFE match: {'YES' if match else 'NO'}")
            return match
        else:
            print(f"  [FAIL] RNAfold failed")
            return False
    else:
        print(f"  [FAIL] Tool failed to return sequence ({elapsed:.1f}s)")
        return False

if __name__ == "__main__":
    tools = {
        "Baseline": benchmark_tools.run_baseline,
        "RNAinverse": benchmark_tools.run_rnainverse,
        "NEMO": benchmark_tools.run_nemo,
        "eM2dRNAs": benchmark_tools.run_em2drnas,
        "DesiRNA": benchmark_tools.run_desirna,
        "LEARNA": benchmark_tools.run_learna,
        "NUPACK": benchmark_tools.run_nupack,
    }
    
    # We want to see raw output for failed tools
    # Let's override the runners to be more verbose in this script
    import subprocess
    import os

    def run_learna_verbose(structure, timeout_sec=20):
        learna_dir = "/home/maxyle/RNA/RNA_design_tools/LEARNA"
        learna_bin = "/home/maxyle/RNA/RNA_design_tools/LEARNA/venv/bin/learna"
        cmd = [learna_bin, "--target_structure", structure, "--timeout", str(timeout_sec)]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec + 10, cwd=learna_dir)
        print(f"DEBUG LEARNA STDOUT:\n{res.stdout}")
        print(f"DEBUG LEARNA STDERR:\n{res.stderr}")
        return benchmark_tools.run_learna(structure, timeout_sec)

    def run_desirna_verbose(structure, timeout_sec=20):
        # We'll just rely on the existing run_desirna but maybe print something if we can
        # Actually better to just run the command manually here
        print("DEBUG DesiRNA: Running manually...")
        desirna_dir = "/home/maxyle/RNA/RNA_design_tools/DesiRNA/DesiRNA-main"
        desirna_python = os.path.join(desirna_dir, "venv", "bin", "python3")
        desirna_script = os.path.join(desirna_dir, "DesiRNA.py")
        input_content = f">name\nBenchmark\n>seq_restr\nNNNNNNNNNNNN\n>sec_struct\n{structure}\n"
        with open("/tmp/test_desirna_simple.txt", "w") as f: f.write(input_content)
        cmd = [desirna_python, desirna_script, '-f', "/tmp/test_desirna_simple.txt", '-R', '1', '-e', '10', '-t', str(timeout_sec)]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec + 10, cwd=desirna_dir)
        print(f"DEBUG DesiRNA STDOUT:\n{res.stdout}")
        print(f"DEBUG DesiRNA STDERR:\n{res.stderr}")
        return benchmark_tools.run_desirna(structure, timeout_sec)

    tools["LEARNA"] = run_learna_verbose
    tools["DesiRNA"] = run_desirna_verbose
    
    results = {}
    for name, func in tools.items():
        results[name] = verify_tool(name, func)
        print()
        
    print("=" * 30)
    print("SUMMARY")
    for name, success in results.items():
        print(f"{name:12s}: {'SUCCESS' if success else 'FAILED'}")
