import RNA
import random
import os
import glob

def parseSS(s):
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

GC_COMPATIBLE = {"G": "C", "C": "G"}

def MFE(sequence):
    md = RNA.md()
    md.uniq_ML = 1 # Using supervisor's setting
    fc = RNA.fold_compound(sequence, md)
    (ss, mfe) = fc.mfe()
    return ss

def baseline(s):
    res = ["A" for _ in s]
    nts = list(GC_COMPATIBLE.keys())
    for i, j in parseSS(s):
        c = random.choice(nts)
        res[i], res[j] = c, GC_COMPATIBLE[c]
    return "".join(res)

def main():
    files = glob.glob("/home/maxyle/RNA/*.txt") + glob.glob("/home/maxyle/RNA/output/*.txt")
    print(f"{'Filename':<40} | {'Success Rate'}")
    print("-" * 60)
    
    for fpath in files:
        structures = []
        try:
            with open(fpath, "r") as f:
                for line in f:
                    l = line.strip()
                    if l and not l.startswith("#") and set(l).issubset({"(", ")", "."}) and len(l) > 10:
                        structures.append(l)
        except:
            continue
            
        if not structures:
            continue
            
        designs = 0
        total = len(structures)
        for s in structures:
            seq = baseline(s)
            if MFE(seq) == s:
                designs += 1
        
        rate = (designs / total) * 100
        print(f"{os.path.basename(fpath):<40} | {designs}/{total} ({rate:.1f}%)")

if __name__ == "__main__":
    main()
