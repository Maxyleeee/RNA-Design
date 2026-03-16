import subprocess

# Let's test varying degrees of constraints
test_structures = {
    "zero_designs": "((((.((((((((....))))))))((((((((....))))))))((((((((....))))))))))))", # Extremely dense, 71 chars
    "too_many_designs": ".........................(((....))).........................", # Very loose, 60 chars
    "few_designs": "((((...(((((....)))))...))))", # Moderate, 28 chars
    "edge_case_some_tools_fail": "(((((..(((....))).(((((....)))))..)))))...." # 43 chars, partially constrained
}

def test_linearbpdesign(structure, nb_samples=10000):
    """
    Tests LinearBPDesign execution explicitly for targeted edge cases.
    Inputs:
        structure: str - Secondary structural bounds constraint
        nb_samples: int - Max limit threshold for exploration
    Outputs:
        int - Total quantity of valid unique negative design responses
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
    Evaluates NEMO binary on specific test instances.
    Inputs:
        structure: str - Secondary structural bounds constraint
    Outputs:
        bool - Boolean resolving successful generation of an identical structural match sequence
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/NEMO/nemo", structure]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        output = res.stdout.strip()
        if "Solution found" in output or "SUCCESS" in output or any(c in 'ACGU' for c in output if len(c)>10):
             return True
        valid_seqs = [line for line in output.split('\n') if len(line.strip().split()[0]) == len(structure) and all(c in 'ACGU' for c in line.strip().split()[0])]
        if valid_seqs:
            return True
        return False
    except Exception as e:
        return False

def test_em2drnas(structure):
    """
    Evaluates genetic performance stringency of eM2dRNAs.
    Inputs:
        structure: str - Secondary structural requirement layout
    Outputs:
        int - Count of successful generations retrieved by testing parameters
    """
    cmd = ["/home/maxyle/RNA/RNA_design_tools/eM2dRNAs/src/e_m2dRNAs", "1", structure, "50", "10", "TURNER2004", "1", "1"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        output = res.stdout.strip()
        valid_seqs = [line for line in output.split('\n') if all(c in 'ACGU' for c in line.strip()) and len(line.strip()) == len(structure)]
        return len(set(valid_seqs))
    except subprocess.TimeoutExpired:
        return 0
    except Exception as e:
        return 0

for name, s in test_structures.items():
    print(f"Testing {name}: {s} (length: {len(s)})")
    lbp = test_linearbpdesign(s)
    nemo = test_nemo(s)
    em_count = test_em2drnas(s)
    print(f"LinearBPDesign unique generated: {lbp}")
    print(f"NEMO success: {nemo}")
    print(f"eM2dRNAs generated: {em_count}")
    print("-" * 40)
