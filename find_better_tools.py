import subprocess
import glob

def parse_lbp(s):
    """
    Parses an RNA structure against LinearBPDesign algorithm.
    Inputs:
        s: str - RNA dot-bracket notation string
    Outputs:
        int - Count of valid derived sequences based on structural boundaries
    """
    cmd = ["python3", "/home/maxyle/RNA/RNA_design_tools/RNAInverse/LinearBPDesign.py", "--nb_samples", "100", "--structure", s]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        valid = [line for line in res.stdout.split('\n') if all(c in 'ACGU' for c in line.strip()) and len(line.strip()) == len(s)]
        return len(set(valid))
    except: return 0

def parse_nemo(s):
    """
    Parses an RNA structure using the NEMO design software.
    Inputs:
        s: str - Target RNA dot-bracket sequence
    Outputs:
        bool - Successful mapping resolution
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/NEMO/nemo/nemo", s]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        return "NMC:" in res.stdout
    except: return False

def parse_em2(s):
    """
    Executes an MFE constraint loop against eM2dRNAs engine.
    Inputs:
        s: str - Desired nested structural layout notation
    Outputs:
        bool - Determines genetic algorithm matching success
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/eM2dRNAs/src/e_m2dRNAs", "1", s, "10", "5", "TURNER2004", "1", "1"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=8)
        return "@FOUND:" in res.stdout
    except: return False

def evaluate():
    """
    Orchestrates continuous comparative evaluations resolving tool disparities.
    Inputs:
        None
    Outputs:
        None - Generates side-by-side performance file comparisons
    """
    files = glob.glob("output/structures_with_*.txt")
    structures = []
    for f in files:
        with open(f, 'r') as file:
            structures.extend([line.strip() for line in file if line.strip() and not line.startswith('#')])
            
    structures = list(set(structures))
    
    everyone_succeeds = []
    
    for i, s in enumerate(structures):
        if len(s) > 60: continue # Focus on small valid structures to guarantee eM2dRNAs doesn't timeout
        
        lbp = parse_lbp(s)
        nemo = parse_nemo(s)
        em2 = parse_em2(s)
        
        print(f"Len: {len(s)} | LBP: {lbp}, NEMO: {nemo}, eM2: {em2}")
        
        if lbp == 0 and (nemo or em2):
            everyone_succeeds.append((s, lbp, nemo, em2))
            with open("case6_results.txt", "a") as f:
                f.write(f"Structure: {s}\nLBP: {lbp} | NEMO: {nemo} | eM2: {em2}\n\n")
            if len(everyone_succeeds) >= 10:
                break


if __name__ == "__main__":
    evaluate()
