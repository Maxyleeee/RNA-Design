import RNA
import random
import os

def parseSS(s):
    stack = []
    pairs = []
    for i, char in enumerate(s):
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            pairs.append((j, i))
    return pairs

COMPATIBLE = {"G": "C", "C": "G", "A": "U", "U": "A"}
GC_COMPATIBLE = {"G": "C", "C": "G"}

def MFE(sequence):
    md = RNA.md()
    md.uniq_ML = 1
    fc = RNA.fold_compound(sequence, md)
    (ss, mfe) = fc.mfe()
    return ss

def baseline(s, gc_only=True, a_loops=True):
    if a_loops:
        res = ["A" for _ in s]
    else:
        res = [random.choice(["A", "U", "C", "G"]) for _ in s]
    
    target_compatible = GC_COMPATIBLE if gc_only else COMPATIBLE
    nts = list(target_compatible.keys())
    for i, j in parseSS(s):
        c = random.choice(nts)
        res[i], res[j] = c, target_compatible[c]
    return "".join(res)

def main():
    log_file = "structures_motif_h4.txt"
    structures = []
    with open(log_file, "r") as f:
        for line in f:
            l = line.strip()
            if l and not l.startswith("#") and set(l).issubset({"(", ")", "."}):
                structures.append(l)

    print(f"Comparison of Loop Alphabets (L=50, h=3, GC-Helices):")
    for mode_name, a_loops in [("A-only Loops (Supervisor style)", True), ("Random Loops (Natural complexity)", False)]:
        designs = 0
        total = len(structures)
        for s in structures:
            seq = baseline(s, gc_only=True, a_loops=a_loops)
            if MFE(seq) == s:
                designs += 1
        print(f"  {mode_name}: {designs}/{total} ({100.*designs/total:.1f}%) OK")

if __name__ == "__main__":
    main()
