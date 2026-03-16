import subprocess
import argparse
import sys
import os

def run_rnainverse(structure, timeout_sec=10):
    """
    Runs RNAinverse on a dot-bracket structure.
    Returns (sequence, success).
    """
    try:
        process = subprocess.Popen(
            ['RNAinverse'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        try:
            stdout, stderr = process.communicate(input=structure + "\n", timeout=timeout_sec)
        except subprocess.TimeoutExpired:
            process.kill()
            process.communicate()
            return None, False
        if process.returncode != 0:
            return None, False
        lines = stdout.strip().split('\n')
        if not lines or not lines[0]:
            return None, False
        parts = lines[0].split()
        if len(parts) >= 1:
            sequence = parts[0]
            if len(sequence) == len(structure) and all(c in 'ACGU' for c in sequence):
                return sequence, True
        return None, False
    except Exception as e:
        print(f"  [Error in RNAinverse]: {e}")
        return None, False

def run_rnafold(sequence, timeout_sec=5):
    """
    Runs RNAfold on a sequence.
    Returns (mfe_structure, success).
    """
    try:
        process = subprocess.Popen(
            ['RNAfold', '--noPS'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        try:
            stdout, stderr = process.communicate(input=sequence + "\n", timeout=timeout_sec)
        except subprocess.TimeoutExpired:
            process.kill()
            process.communicate()
            return None, False
        if process.returncode != 0:
            return None, False
        lines = stdout.strip().split('\n')
        if len(lines) >= 2:
            # RNAfold output: line 1 = sequence, line 2 = "structure (energy)"
            fold_line = lines[1].split()
            if len(fold_line) >= 1:
                return fold_line[0], True
        return None, False
    except Exception as e:
        print(f"  [Error in RNAfold]: {e}")
        return None, False

def verify_file(filepath, max_structures=10):
    if not os.path.exists(filepath):
        print(f"Error: File '{filepath}' not found.")
        return

    print(f"Verifying up to {max_structures} structures in: {filepath}")
    
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and set(line.strip()).issubset({'(', ')', '.'})]
        
    lines = lines[:max_structures]
    total = len(lines)
    
    print(f"Found {total} valid structures to verify.\n")
    
    results = []
    
    for i, structure in enumerate(lines):
        print(f"[{i+1}/{total}] Target: {structure}")
        seq, success_inv = run_rnainverse(structure)
        if success_inv and seq:
            print(f"  RNAinverse => {seq}")
            mfe_struct, success_fold = run_rnafold(seq)
            if success_fold and mfe_struct:
                match = (mfe_struct == structure)
                print(f"  RNAfold    => {mfe_struct}")
                print(f"  Match: {'YES' if match else 'NO'}")
                results.append((structure, seq, mfe_struct, match))
            else:
                print("  RNAfold: FAILED")
                results.append((structure, seq, None, False))
        else:
            print("  RNAinverse: FAILED (no valid sequence found)")
            results.append((structure, None, None, False))
            
    # Summary
    matches = sum(1 for r in results if r[3])
    inverses = sum(1 for r in results if r[1] is not None)
    
    print(f"\n{'='*50}")
    print(f"SUMMARY")
    print(f"{'='*50}")
    print(f"Total structures tested:    {total}")
    print(f"RNAinverse found sequence:  {inverses}/{total}")
    print(f"RNAfold exact match:        {matches}/{total}")
    if inverses > 0:
        print(f"Match rate (of designed):   {matches}/{inverses} ({matches/inverses*100:.1f}%)")
    
    # Save results
    out_filepath = filepath.replace('.txt', '_verification.txt')
    with open(out_filepath, 'w') as f:
        f.write(f"Verification Results for {filepath}\n")
        f.write(f"Total Structures Tested: {total}\n\n")
        f.write(f"RNAinverse success: {inverses}/{total}\n")
        f.write(f"RNAfold exact match: {matches}/{total}\n")
        if inverses > 0:
            f.write(f"Match rate (of designed): {matches}/{inverses} ({matches/inverses*100:.1f}%)\n")
        f.write("\n")
        
        for target, seq, mfe, match in results:
            f.write(f"Target: {target}\n")
            if seq:
                f.write(f"Seq:    {seq}\n")
                if mfe:
                    f.write(f"MFE:    {mfe}\n")
                    f.write(f"Match:  {'YES' if match else 'NO'}\n")
                else:
                    f.write("RNAfold: FAILED\n")
            else:
                f.write("RNAinverse: FAILED\n")
            f.write("-" * 50 + "\n")
            
    print(f"\nDetailed results saved to: {out_filepath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify RNA structure designability: RNAinverse -> RNAfold round-trip.")
    parser.add_argument("filepath", help="Path to file with dot-bracket structures (one per line)")
    parser.add_argument("--max", type=int, default=10, help="Max structures to test (default: 10)")
    args = parser.parse_args()
    verify_file(args.filepath, args.max)
