import sys

filepath = 'output/structures_with_multi_motif_len100_loose.txt'
with open(filepath, 'r') as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith('#'): continue
        valid = True
        count = 0
        for c in s:
            if c == '(': count += 1
            elif c == ')':
                count -= 1
                if count < 0:
                    valid = False
                    break
        print(f"Len: {len(s)} | Valid: {valid and count == 0} | String: {s}")
