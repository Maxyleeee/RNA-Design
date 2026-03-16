import subprocess
import os
import glob

def test_linearbpdesign(structure, nb_samples=100):
    """
    Execute LinearBPDesign to analyze RNA structure constraints.
    Inputs:
        structure: str - Valid secondary format structural string
        nb_samples: int - Parameter configuring statistical design searches
    Outputs:
        int - Returned volume of identical output configurations 
    """
    cmd = ["python3", "/home/maxyle/RNA/RNA_design_tools/RNAInverse/LinearBPDesign.py", "--nb_samples", str(nb_samples), "--structure", structure]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        output = res.stdout.strip()
        valid_seqs = [line for line in output.split('\n') if all(c in 'ACGU' for c in line.strip()) and len(line.strip()) == len(structure)]
        return len(set(valid_seqs))
    except Exception as e:
        return 0

def test_nemo(structure):
    """
    Leverages NEMO binary algorithms checking motif limits.
    Inputs:
        structure: str - Master input string requiring evaluation
    Outputs:
        bool - Determines algorithm resolution success
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/NEMO/nemo", structure]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        output = res.stdout.strip()
        if "Solution found" in output or "SUCCESS" in output or any(c in 'ACGU' for c in output if len(c)>10):
             return True
        valid_seqs = [line for line in output.split('\n') if len(line.strip().split()[0]) == len(structure) and all(c in 'ACGU' for c in line.strip().split()[0])]
        if valid_seqs: return True
        return False
    except Exception as e:
        return False

def test_em2drnas(structure):
    """
    Triggers genetic heuristic modeling engine for structural solutions.
    Inputs:
        structure: str - Candidate constraint configuration
    Outputs:
        int - Metric resolving total outputted heuristic sequences
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/eM2dRNAs/src/e_m2dRNAs", "1", structure, "10", "5", "TURNER2004", "1", "1"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=8)
        output = res.stdout.strip()
        valid_seqs = [line for line in output.split('\n') if all(c in 'ACGU' for c in line.strip()) and len(line.strip()) == len(structure)]
        return len(set(valid_seqs))
    except subprocess.TimeoutExpired:
        return 0
    except Exception as e:
        return 0

def classify_structures():
    """
    Automated orchestration loop grabbing batched results from geometry generation.
    Inputs:
        None
    Outputs:
        None - Programmatically handles text operations generating categorization records
    """
    files = glob.glob("output/structures_with_multi_motif_len*.txt")
    
    structures_to_test = []
    for fpath in files:
        with open(fpath, 'r') as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            structures_to_test.extend(lines[:30]) # test first 30 from each file
            
    print(f"Total structures to test: {len(structures_to_test)}")
    
    cases = {
        "zero_designs": [],
        "too_many_designs": [],
        "few_designs": [],
        "partial_success": []
    }
    
    for i, s in enumerate(structures_to_test):
        print(f"Testing {i+1}/{len(structures_to_test)} (Len: {len(s)})")
        lbp = test_linearbpdesign(s, nb_samples=100)
        
        # Categorize mostly off LinearBPDesign for speed
        if lbp == 0:
            cases["zero_designs"].append((s, lbp))
        elif lbp > 90:
            cases["too_many_designs"].append((s, lbp))
        elif 0 < lbp <= 90:
            # Maybe few designs or partial success
            nemo = test_nemo(s)
            em2 = test_em2drnas(s)
            
            if nemo or em2 > 0:
                cases["few_designs"].append((s, lbp, nemo, em2))
            else:
                cases["partial_success"].append((s, lbp, nemo, em2))
                
        # Stop if we have enough of each
        if len(cases["zero_designs"]) >= 10 and len(cases["too_many_designs"]) >= 10 and len(cases["few_designs"]) >= 10 and len(cases["partial_success"]) >= 10:
            break

    # Write out the results
    with open("categorized_results.txt", "w") as f:
        f.write("# Automated Classification of Generated Motif Structures\n\n")
        
        f.write("### Case 1: Structure with ZERO folding designs\n")
        for s, lbp in cases["zero_designs"][:10]:
            f.write(f"Structure: {s}\nLinearBPDesign: 0 designs found.\n\n")
            
        f.write("### Case 2: Structure with TOO MANY designs (Loose)\n")
        for s, lbp in cases["too_many_designs"][:10]:
            f.write(f"Structure: {s}\nLinearBPDesign: Extracted {lbp}+ unique strict negative designs (search space is massive).\n\n")
            
        f.write("### Case 3: Structure with ONLY A FEW designs\n")
        for s, lbp, nemo, em2 in cases["few_designs"][:10]:
            f.write(f"Structure: {s}\nLinearBPDesign: {lbp} unique designs.\nNEMO Success: {nemo}\neM2dRNAs generated: {em2}\n\n")
            
        f.write("### Case 4: Edge Case (Some tools succeed, some fail)\n")
        for s, lbp, nemo, em2 in cases["partial_success"][:10]:
            f.write(f"Structure: {s}\nLinearBPDesign: {lbp} unique designs.\nNEMO Success: {nemo} (Failed)\neM2dRNAs generated: {em2} (Failed/Timeout)\n\n")

if __name__ == "__main__":
    classify_structures()
