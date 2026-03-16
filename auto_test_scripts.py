import subprocess
import os
import sys

def test_linearbpdesign(structure, nb_samples=10):
    """
    Tests an RNA structure using LinearBPDesign.
    Inputs:
        structure: str - RNA sequence in dot-bracket notation
        nb_samples: int - Number of sequence samples to generate
    Outputs:
        tuple (int, str) - Count of valid sequences found, and raw stdout log
    """
    cmd = ["python3", "/home/maxyle/RNA/RNA_design_tools/RNAInverse/LinearBPDesign.py", "--nb_samples", str(nb_samples), "--structure", structure]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        output = res.stdout.strip()
        # Count non-empty lines that look like RNA sequences (mostly ACGU)
        valid_seqs = [line for line in output.split('\n') if all(c in 'ACGU' for c in line.strip()) and len(line.strip()) == len(structure)]
        return len(valid_seqs), output
    except Exception as e:
        return 0, str(e)

def test_nemo(structure):
    """
    Tests an RNA structure using NEMO.
    Inputs:
        structure: str - RNA sequence in dot-bracket notation
    Outputs:
        tuple (bool, str) - True if a valid solution was found, and raw stdout log
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/NEMO/nemo", structure]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        output = res.stdout.strip()
        if "Solution found" in output or "SUCCESS" in output or any(c in 'ACGU' for c in output if len(c)>10):
             return True, output
        valid_seqs = [line for line in output.split('\n') if len(line.strip().split()[0]) == len(structure) and all(c in 'ACGU' for c in line.strip().split()[0])]
        if valid_seqs:
            return True, output
        return False, output
    except Exception as e:
        return False, str(e)

def test_em2drnas(structure):
    """
    Tests an RNA structure using eM2dRNAs.
    Inputs:
        structure: str - RNA sequence in dot-bracket notation
    Outputs:
        tuple (bool, str) - True if a valid solution was found, and raw stdout log
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/eM2dRNAs/src/e_m2dRNAs", "1", structure, "10", "10", "TURNER2004", "1", "1"]
    try:
        # Give it 10 seconds timeout via its internal parameter, plus 15s subprocess timeout
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        output = res.stdout.strip()
        valid_seqs = [line for line in output.split('\n') if all(c in 'ACGU' for c in line.strip()) and len(line.strip()) == len(structure)]
        return len(valid_seqs) > 0, output
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)

if __name__ == "__main__":
    filepath = sys.argv[1]
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    structures = lines[:3]  # Take top 3 to test
    
    print(f"Testing {len(structures)} structures from {filepath}\n")
    
    results = []
    
    for idx, s in enumerate(structures):
        print(f"--- Structure {idx+1}: {s}")
        print(f"Length: {len(s)}")
        
        counts_lbp, out_lbp = test_linearbpdesign(s, nb_samples=50)
        succ_nemo, out_nemo = test_nemo(s)
        succ_em2, out_em2 = test_em2drnas(s)
        
        print(f"LinearBPDesign: {counts_lbp} samples found.")
        print(f"NEMO success: {succ_nemo}")
        print(f"eM2dRNAs success: {succ_em2}")
        print()
