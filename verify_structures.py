
import sys

def verify_structural_constraints(filename, h, theta):
    print(f"Verifying {filename} with h={h}, theta={theta}")
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    for i, structure in enumerate(lines):
        # 1. Check balanced
        if structure.count('(') != structure.count(')'):
            print(f"FAILED Line {i+1}: Unbalanced parens. {structure.count('(')} vs {structure.count(')')}")
            return False
            
        # 2. Check contiguous pair blocks?
        # The grammar forces (*h ... )*h.
        # This implies that we should not see isolated '('. 
        # Actually it forces every '(' to be part of a block of size h?
        # Yes, if generating strictly from grammar.
        # Let's count runs of '(' and ')'
        
        runs_open = [len(s) for s in structure.replace('.', ' ').replace(')', ' ').split() if '(' in s]
        runs_close = [len(s) for s in structure.replace('.', ' ').replace('(', ' ').split() if ')' in s]
        
        # Note: Adjacent helices (((...)))(((...))) might merge runs?
        # (((... (( ...
        # If we have ((())) , run is 3.
        # If we have ((( ... ))), run is 3.
        # If we have (((((()))))), run is 6.
        # So runs should be multiples of h.
        
        for r in runs_open:
            if r % h != 0:
                print(f"FAILED Line {i+1}: Run of '(' length {r} not divisible by h={h}. Run: {r}")
                # print(structure)
                return False
                
        for r in runs_close:
            if r % h != 0:
                print(f"FAILED Line {i+1}: Run of ')' length {r} not divisible by h={h}. Run: {r}")
                # print(structure)
                return False
                
        # 3. Check loops (theta)
        # Any pair (i, j) must have j - i - 1 >= theta?
        # No, "theta=min dist of paired positions".
        # This means j - i >= theta? Or j - i > theta?
        # Usually minimal loop size m means j-i-1 >= m.
        # User said: "theta=min dist of paired positions".
        # If pos are 0-indexed.
        # Pair (u, v). v - u >= theta.
        # Grammar enforces T length >= theta.
        # T is content. So v = u + 1 + length(T).
        # v - u = length(T) + 1.
        # If length(T) >= theta, then v - u >= theta + 1.
        # So min dist is theta+1?
        # User wrote "theta=min dist of paired positions".
        # If theta=3, then v-u >= 3.
        # Example: ( . . ) -> u=0, v=3. v-u=3. Content len 2.
        # So content length must be >= theta - 1?
        # Or content length >= theta?
        # Our code enforced `countT(n)` where `n >= theta`.
        # So we used content length >= theta.
        # So v - u >= theta + 1.
        # This is safe (stricter or equal).
        
        stack = []
        for idx, char in enumerate(structure):
            if char == '(':
                stack.append(idx)
            elif char == ')':
                if not stack:
                    print(f"FAILED Line {i+1}: Extra closing paren")
                    return False
                start = stack.pop()
                dist = idx - start
                # We enforced content length >= theta
                # content length = dist - 1.
                # So dist - 1 >= theta => dist >= theta + 1.
                if (dist - 1) < theta:
                     print(f"FAILED Line {i+1}: Pair ({start},{idx}) has content length {dist-1} < theta={theta}")
                     return False
                     
    print("All structures passed verification.")
    return True

if __name__ == "__main__":
    verify_structural_constraints("output/structures_standard_L100.txt", 3, 3)
