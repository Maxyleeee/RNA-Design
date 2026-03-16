import subprocess
import sys
import os
import argparse
import time
import tempfile

###############################################################################
# Tool runners — each returns (sequence_or_None, success_bool)
###############################################################################

def run_rnainverse(structure, timeout_sec=30):
    """
    Runs ViennaRNA RNAinverse on a dot-bracket structure.
    RNAinverse uses stochastic local search to find a sequence
    whose MFE fold matches the target structure.
    """
    try:
        proc = subprocess.Popen(
            ['RNAinverse', '-R-1'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True
        )
        try:
            stdout, _ = proc.communicate(input=structure + "\n", timeout=timeout_sec)
        except subprocess.TimeoutExpired:
            proc.kill(); proc.communicate()
            return None, False
        if proc.returncode != 0:
            return None, False
        lines = stdout.strip().split('\n')
        if not lines or not lines[0]:
            return None, False
        parts = lines[0].split()
        if parts and len(parts[0]) == len(structure) and all(c in 'ACGU' for c in parts[0]):
            return parts[0], True
        return None, False
    except Exception:
        return None, False


def run_nemo(structure, timeout_sec=30):
    """
    Runs NEMO (Monte Carlo search) on a dot-bracket structure.
    Extracts the designed sequence from its output.
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/NEMO/nemo/nemo", structure]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
        output = res.stdout.strip()
        # NEMO output format: "NMC: SEQUENCE score"
        for line in output.split('\n'):
            line = line.strip()
            if not line:
                continue
            if line.startswith("NMC:"):
                parts = line.split("NMC:")[1].strip().split()
                if parts:
                    seq = parts[0]
                    if len(seq) == len(structure) and all(c in 'ACGU' for c in seq):
                        return seq, True
        return None, False
    except subprocess.TimeoutExpired:
        return None, False
    except Exception:
        return None, False


def run_em2drnas(structure, timeout_sec=30):
    """
    Runs eM2dRNAs (evolutionary multi-objective) on a dot-bracket structure.
    eM2dRNAs uses a genetic algorithm with Turner energy model.
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/eM2dRNAs/src/e_m2dRNAs",
           "1", structure, "50", "10", "TURNER2004", "1", "1"]
    try:
        # eM2dRNAs may return exit code 1 but still produce valid output
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
        output = res.stdout.strip()
        for line in output.split('\n'):
            # Try each token on the line — eM2dRNAs can output seq alongside other data
            for token in line.strip().split():
                if len(token) == len(structure) and all(c in 'ACGU' for c in token):
                    return token, True
        return None, False
    except subprocess.TimeoutExpired:
        return None, False
    except Exception:
        return None, False


def run_desirna(structure, timeout_sec=60):
    """
    Runs DesiRNA (Replica Exchange Monte Carlo) on a dot-bracket structure.
    DesiRNA requires a formatted input file and its own virtualenv.
    """
    desirna_dir = "/home/maxyle/RNA/RNA_design_tools/DesiRNA/DesiRNA-main"
    desirna_python = os.path.join(desirna_dir, "venv", "bin", "python3")
    desirna_script = os.path.join(desirna_dir, "DesiRNA.py")

    # Create temp input file in DesiRNA's expected format
    seq_restr = 'N' * len(structure)
    input_content = f">name\nBenchmark\n>seq_restr\n{seq_restr}\n>sec_struct\n{structure}\n"

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, dir='/tmp') as tf:
            tf.write(input_content)
            tmp_path = tf.name

        cmd = [desirna_python, desirna_script, '-f', tmp_path,
               '-R', '1', '-e', '10', '-t', str(timeout_sec)]
        # Use stdin=subprocess.DEVNULL to prevent hangs
        # Combine stdout and stderr just in case
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec + 5,
                             cwd=desirna_dir, stdin=subprocess.DEVNULL)
        output = (res.stdout or "") + "\n" + (res.stderr or "")

        # Parse DesiRNA output — look for a line that is purely ACGU and matches length
        # Use replace('\r', '') to handle progress bars
        for line in output.replace('\r', '\n').split('\n'):
            line = line.strip()
            if line and len(line) == len(structure) and all(c in 'ACGU' for c in line):
                return line, True

        return None, False
    except subprocess.TimeoutExpired:
        return None, False
    except Exception:
        return None, False
    finally:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except:
                pass


def run_learna(structure, timeout_sec=60):
    """
    Runs LEARNA on a dot-bracket structure.
    Uses pre-trained weights from models/224_0_1.
    Note: LEARNA is run using a dedicated Conda environment 'learna_env'.
    """
    learna_dir = "/home/maxyle/RNA/RNA_design_tools/LEARNA"
    # Using the conda environment created for LEARNA
    learna_bin = "/home/maxyle/miniconda/envs/learna_env/bin/learna"
    
    cmd = [learna_bin, "--target_structure", structure, "--timeout", str(timeout_sec)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec + 10, cwd=learna_dir)
        output = res.stdout.strip()
        
        for line in output.split('\n'):
            line = line.strip()
            if '|' in line:
                parts = [p.strip() for p in line.split('|')]
                if len(parts) >= 7:
                    # In our custom table format: | Id | time | hamming | rel_hamming | sequence | structure |
                    # parts will be ['', 'Id', 'time', 'hamming', 'rel_hamming', 'sequence', 'structure', '']
                    # sequence is at index 5, which is parts[-3] if there's a trailing empty string
                    for p in parts:
                        if len(p) == len(structure) and all(c in 'ACGU' for c in p):
                            return p, True
                            
        return None, False
    except Exception as e:
        print(f"[DEBUG] LEARNA error: {e}")
        return None, False


def run_nupack(structure, timeout_sec=60):
    """
    Runs NUPACK 4.0 Design on a dot-bracket structure.
    Executes in a subprocess using Python 3 to ensure the nupack module is found.
    """
    script = f'''
import sys
try:
    import nupack
    config = nupack.Model(material="rna")
    domain = nupack.Domain("N" * {len(structure)}, name="d1")
    strand = nupack.TargetStrand([domain], name="s1")
    complex_target = nupack.TargetComplex([strand], "{structure}", name="c1")
    tube = nupack.TargetTube(on_targets={{complex_target: 1e-6}}, name="t1")
    design = nupack.Design(tubes=[tube], model=config)
    results = design.run(trials=1)
    
    if results:
        seq = None
        try:
            for d_obj, d_seq in results[0].domains.items():
                if len(str(d_seq)) == {len(structure)}:
                    seq = str(d_seq)
                    break
        except Exception:
            pass
        
        if not seq:
            try:
                analysis = results[0].to_analysis()
                for s_name in analysis.strands:
                    seq = str(analysis.strands[s_name])
                    break
            except Exception:
                pass
                
        if seq and len(seq) == {len(structure)} and all(c in "ACGU" for c in seq):
            print("NUPACK_SEQ:" + seq)
except Exception as e:
    print("NUPACK_ERROR:", e, file=sys.stderr)
'''
    try:
        # Run using the system python3 which has nupack installed
        res = subprocess.run(
            ["/usr/bin/python3", "-c", script], 
            capture_output=True, text=True, timeout=timeout_sec
        )
        for line in res.stdout.splitlines():
            if line.startswith("NUPACK_SEQ:"):
                return line.split("NUPACK_SEQ:")[1].strip(), True
        return None, False
    except subprocess.TimeoutExpired:
        return None, False
    except Exception as e:
        print(f"[DEBUG] Subprocess NUPACK error: {e}")
        return None, False

def parseSS(s):
    """Parses dot-bracket secondary structure and returns pairs."""
    stack = []
    pairs = []
    for i, char in enumerate(s):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    return pairs

def run_baseline(structure, timeout_sec=10):
    """
    Runs the baseline method requested by the user.
    Always places 'A' for unpaired bases, and random pairs for paired bases.
    """
    import random
    COMPATIBLE = {"G":"C", "C":"G"}
    NTS = list(COMPATIBLE.keys())
    
    res = ["A" for _ in structure]
    for i, j in parseSS(structure):
        c = random.choice(NTS)
        res[i], res[j] = c, COMPATIBLE[c]
    
    seq = "".join(res)
    return seq, True


def calculate_hamming_distance(s1, s2):
    """Calculates Hamming distance between two dot-bracket structures."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


###############################################################################
# Verification — RNAfold
###############################################################################

def run_rnafold(sequence, timeout_sec=10):
    """
    Runs ViennaRNA RNAfold on a sequence.
    Returns the MFE secondary structure.
    """
    try:
        proc = subprocess.Popen(
            ['RNAfold', '--noPS'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True
        )
        try:
            stdout, _ = proc.communicate(input=sequence + "\n", timeout=timeout_sec)
        except subprocess.TimeoutExpired:
            proc.kill(); proc.communicate()
            return None, False
        if proc.returncode != 0:
            return None, False
        lines = stdout.strip().split('\n')
        if len(lines) >= 2:
            fold_line = lines[1].split()
            if fold_line:
                return fold_line[0], True
        return None, False
    except Exception:
        return None, False


###############################################################################
# Main benchmark loop
###############################################################################

def benchmark_file(filepath, max_structures=10):
    if not os.path.exists(filepath):
        print(f"Error: File '{filepath}' not found.")
        return

    print(f"=" * 70)
    print(f"BENCHMARK: {filepath}")
    print(f"=" * 70)

    with open(filepath, 'r') as f:
        raw_lines = [l.strip() for l in f if l.strip()]
    
    lines = []
    motif_counts = {} # Map structure index -> count
    
    for i, l in enumerate(raw_lines):
        if "," in l:
            struct, count_str = l.split(",", 1)
            struct = struct.strip()
            # Basic validation
            if set(struct).issubset({'(', ')', '.'}):
                lines.append(struct)
                try:
                    motif_counts[len(lines)-1] = int(count_str.strip())
                except:
                    motif_counts[len(lines)-1] = 0
        elif set(l).issubset({'(', ')', '.'}):
            lines.append(l)
            motif_counts[len(lines)-1] = 0
            
    lines = lines[:max_structures]
    total = len(lines)
    print(f"Structures to test: {total}\n")

    TOOLS = ["Baseline", "RNAinverse", "NEMO", "eM2dRNAs", "DesiRNA", "LEARNA", "NUPACK"]
    tool_funcs = {
        "Baseline": run_baseline,
        "RNAinverse": run_rnainverse,
        "NEMO": run_nemo,
        "eM2dRNAs": run_em2drnas,
        "DesiRNA": run_desirna,
        "LEARNA": run_learna,
        "NUPACK": run_nupack,
    }

    # Stats tracking
    stats = {t: {
        "designed": 0, 
        "verified": 0, 
        "times_success": [], 
        "times_failure": [], 
        "total_time": 0.0,
        "distances_failure": []
    } for t in TOOLS}
    all_results = []

    for i, structure in enumerate(lines):
        print(f"[{i+1}/{total}] Target (len={len(structure)}): {structure}")
        struct_result = {"target": structure, "tools": {}}

        for tool_name in TOOLS:
            t0 = time.time()
            seq, ok = tool_funcs[tool_name](structure)
            elapsed = time.time() - t0
            stats[tool_name]["total_time"] += elapsed

            if ok and seq:
                stats[tool_name]["designed"] += 1
                mfe, fold_ok = run_rnafold(seq)
                if fold_ok and mfe:
                    match = (mfe == structure)
                    dist = 0 if match else calculate_hamming_distance(structure, mfe)
                    
                    if match:
                        stats[tool_name]["verified"] += 1
                        stats[tool_name]["times_success"].append(elapsed)
                    else:
                        stats[tool_name]["times_failure"].append(elapsed)
                        stats[tool_name]["distances_failure"].append(dist)
                        
                    print(f"  {tool_name:12s}: {seq[:30]}... → MFE match: {'YES' if match else 'NO'} ({elapsed:.1f}s, dist={dist})")
                    struct_result["tools"][tool_name] = {
                        "seq": seq, "mfe": mfe, "match": match, "time": elapsed, "dist": dist
                    }
                else:
                    print(f"  {tool_name:12s}: seq found but RNAfold failed ({elapsed:.1f}s)")
                    struct_result["tools"][tool_name] = {
                        "seq": seq, "mfe": None, "match": False, "time": elapsed, "dist": None
                    }
            else:
                print(f"  {tool_name:12s}: FAILED ({elapsed:.1f}s)")
                struct_result["tools"][tool_name] = {
                    "seq": None, "mfe": None, "match": False, "time": elapsed, "dist": None
                }

        all_results.append(struct_result)
        print()

    # Print summary table
    print("=" * 120)
    print("SUMMARY")
    print("=" * 120)
    header = f"{'Tool':<14} {'Designed':>10} {'Verified':>10} {'Match Rate':>12} {'Avg T(S)':>10} {'Avg T(F)':>10} {'Avg T(A)':>10} {'Avg D(F)':>10}"
    print(header)
    print("-" * len(header))
    for t in TOOLS:
        d = stats[t]["designed"]
        v = stats[t]["verified"]
        rate = f"{v}/{total} ({v/total*100:.0f}%)" if total > 0 else "N/A"
        
        avg_t_s = sum(stats[t]["times_success"]) / len(stats[t]["times_success"]) if stats[t]["times_success"] else 0
        avg_t_f = sum(stats[t]["times_failure"]) / len(stats[t]["times_failure"]) if stats[t]["times_failure"] else 0
        avg_t_a = stats[t]["total_time"] / total if total > 0 else 0
        avg_d_f = sum(stats[t]["distances_failure"]) / len(stats[t]["distances_failure"]) if stats[t]["distances_failure"] else 0
        
        print(f"{t:<14} {d:>5}/{total:<4} {v:>5}/{total:<4} {rate:>12} {avg_t_s:>9.1f}s {avg_t_f:>9.1f}s {avg_t_a:>9.1f}s {avg_d_f:>10.1f}")
    print()

    # Write detailed results
    out_path = filepath.replace('.txt', '_benchmark.txt')
    with open(out_path, 'w') as f:
        f.write(f"Benchmark Results for {filepath}\n")
        f.write(f"Total Structures: {total}\n\n")

        f.write(f"{'Tool':<14} {'Designed':>10} {'Verified':>10} {'Match Rate':>12} {'Avg T(S)':>10} {'Avg T(F)':>10} {'Avg T(A)':>10} {'Avg D(F)':>10}\n")
        f.write("-" * 100 + "\n")
        for t in TOOLS:
            d = stats[t]["designed"]
            v = stats[t]["verified"]
            rate = f"{v}/{total} ({v/total*100:.0f}%)" if total > 0 else "N/A"
            avg_t_s = sum(stats[t]["times_success"]) / len(stats[t]["times_success"]) if stats[t]["times_success"] else 0
            avg_t_f = sum(stats[t]["times_failure"]) / len(stats[t]["times_failure"]) if stats[t]["times_failure"] else 0
            avg_t_a = stats[t]["total_time"] / total if total > 0 else 0
            avg_d_f = sum(stats[t]["distances_failure"]) / len(stats[t]["distances_failure"]) if stats[t]["distances_failure"] else 0
            
            f.write(f"{t:<14} {d:>5}/{total:<4} {v:>5}/{total:<4} {rate:>12} {avg_t_s:>9.1f}s {avg_t_f:>9.1f}s {avg_t_a:>9.1f}s {avg_d_f:>10.1f}\n")
        f.write("\n" + "=" * 70 + "\n\n")

        for r in all_results:
            f.write(f"Target: {r['target']}\n")
            for t in TOOLS:
                tr = r["tools"][t]
                if tr["seq"]:
                    f.write(f"  {t}:\n")
                    f.write(f"    Seq: {tr['seq']}\n")
                    if tr["mfe"]:
                        f.write(f"    MFE: {tr['mfe']}\n")
                        f.write(f"    Match: {'YES' if tr['match'] else 'NO'}\n")
                    else:
                        f.write(f"    RNAfold: FAILED\n")
                else:
                    f.write(f"  {t}: FAILED\n")
                f.write(f"    Time: {tr['time']:.1f}s\n")
            f.write("-" * 60 + "\n")

    print(f"Detailed results saved to: {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Benchmark RNA design tools (RNAinverse, NEMO, eM2dRNAs, DesiRNA, LEARNA) with RNAfold verification."
    )
    parser.add_argument("filepath", help="Path to file with dot-bracket structures")
    parser.add_argument("--max", type=int, default=None, help="Max structures to test (default: All)")
    args = parser.parse_args()
    benchmark_file(args.filepath, args.max)
