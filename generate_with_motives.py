import random
import sys
import os

# Increase recursion depth for deep structures
sys.setrecursionlimit(20000)

def countS(n, cache, wu, ws, wm, h, theta, motifs):
    """
    Counts total weight of structures of length n starting with S.
    Inputs:
        n: int - Length of the sequence
        cache: dict - Memoization dict to store previous calculations
        wu: float - Weight for unpaired bases
        ws: float - Weight for base pairs (stacks)
        wm: float - Weight for motifs from Motif array
        h: int - Minimum helix length
        theta: int - Minimum loop length
        motifs: list - Array of textual motif strings with any number of stars
    Outputs:
        float - Total weight of structures
    S -> . S
    S -> (*h T )*h S  (Weight ws**h)
    S -> eps
    """
    if n == 0: return 1.0
    if ("S", n) in cache: return cache[("S", n)]
    
    val = 0.0
    
    # 1. . S
    if n >= 1:
        val += wu * countS(n-1, cache, wu, ws, wm, h, theta, motifs)
    
    # 2. (*h T )*h S
    if n >= 2*h + theta:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h + 1):
            val += w_stack_start * countT(k, cache, wu, ws, wm, h, theta, motifs) * countS(n - 2*h - k, cache, wu, ws, wm, h, theta, motifs)
            
    cache[("S", n)] = val
    return val

def countT_stars(s, n, cache, wu, ws, wm, h, theta, motifs):
    """
    Counts total weight for s copies of T distributed over length n.
    This effectively resolves the problem of variable holes natively inside sequences.
    T_s -> T T_{s-1}
    T_0 -> eps
    """
    if s == 0:
        return 1.0 if n == 0 else 0.0
        
    key = ("T_stars", s, n)
    if key in cache: return cache[key]
    
    val = 0.0
    if n >= s * theta:
        for k in range(theta, n - (s-1)*theta + 1):
            val += countT(k, cache, wu, ws, wm, h, theta, motifs) * countT_stars(s-1, n-k, cache, wu, ws, wm, h, theta, motifs)
        
    cache[key] = val
    return val

def countT(n, cache, wu, ws, wm, h, theta, motifs):
    """
    Counts total weight of structures inside a helix (Grammar T).
    Inputs:
        n: int - Length of internal sequence block
        cache: dict - Memoization cache
        wu: float - Weight for unpaired bases
        ws: float - Weight for base pairs
        wm: float - Weight for motifs
        h: int - Minimum helix length
        theta: int - Minimum loop length
        motifs: list - Array of textual motif strings
    Outputs:
        float - Total weight of internal loop substructures
    T -> . S
    T -> (*h T )*h . S  (Weight ws**h)
    T -> (*h T )*h (*h T )*h S (Weight ws**2h)
    T -> ( T )  (Weight ws)
    T -> .* theta (Weight wu**theta)
    T -> Motif S (Weight wm)
    """
    if n < theta: return 0.0
    if ("T", n) in cache: return cache[("T", n)]
    
    val = 0.0
    
    # 1. T -> . S
    if n >= 1:
        val += wu * countS(n-1, cache, wu, ws, wm, h, theta, motifs)
        
    # 2. T -> (*h T )*h . S
    if n >= 2*h + theta + 1:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h - 1 + 1):
            rem = n - (2*h + k + 1)
            val += w_stack_start * countT(k, cache, wu, ws, wm, h, theta, motifs) * wu * countS(rem, cache, wu, ws, wm, h, theta, motifs)

    # 3. T -> (*h T )*h (*h T )*h S
    if n >= 4*h + 2*theta:
        w_stack_start_double = ws**(2 * max(0, h-1))
        for k1 in range(theta, n - 4*h - theta + 1):
            cT_k1 = countT(k1, cache, wu, ws, wm, h, theta, motifs)
            for k2 in range(theta, n - 4*h - k1 + 1):
                rem = n - 4*h - k1 - k2
                val += w_stack_start_double * cT_k1 * countT(k2, cache, wu, ws, wm, h, theta, motifs) * countS(rem, cache, wu, ws, wm, h, theta, motifs)

    # 4. T -> ( T )
    if n >= 2 + theta:
        val += ws * countT(n-2, cache, wu, ws, wm, h, theta, motifs)
        
    # 5. T -> .* theta
    if n == theta:
        val += wu**theta

    # 6. Dynamic generic motifs via generalized recursion
    for motif in motifs:
        stars = motif.count('*')
        L = len(motif) - stars
        if n >= L + stars * theta:
            for k in range(L + stars * theta, n + 1):
                val += wm * countT_stars(stars, k - L, cache, wu, ws, wm, h, theta, motifs) * countS(n - k, cache, wu, ws, wm, h, theta, motifs)

    cache[("T", n)] = val
    return val

def generateS(n, cache, wu, ws, wm, h, theta, motifs):
    """
    Generates a top-level random RNA structure string of length n.
    Outputs:
        (str, int) - (Generated dot-bracket string, motif count)
    """
    if n == 0: return "", 0
    total = countS(n, cache, wu, ws, wm, h, theta, motifs)
    if total <= 0: return "." * n, 0
    r = random.random() * total
    
    # 1. . S
    if n >= 1:
        term = wu * countS(n-1, cache, wu, ws, wm, h, theta, motifs)
        if r < term:
            rest_s, rest_c = generateS(n-1, cache, wu, ws, wm, h, theta, motifs)
            return "." + rest_s, rest_c
        r -= term
        
    # 2. (*h T )*h S
    if n >= 2*h + theta:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h + 1):
            term = w_stack_start * countT(k, cache, wu, ws, wm, h, theta, motifs) * countS(n - 2*h - k, cache, wu, ws, wm, h, theta, motifs)
            if r < term:
                inner_s, inner_c = generateT(k, cache, wu, ws, wm, h, theta, motifs)
                rest_s, rest_c = generateS(n - 2*h - k, cache, wu, ws, wm, h, theta, motifs)
                return "(" * h + inner_s + ")" * h + rest_s, inner_c + rest_c
            r -= term
            
    return "." * n, 0

def generateT_stars(s, n, cache, wu, ws, wm, h, theta, motifs):
    """
    Generates an array of inner sequence parts filling s holes.
    Outputs:
        (list of str, int) - (List of structures, total motif count)
    """
    if s == 0:
        return ([] if n == 0 else None), 0
        
    total = countT_stars(s, n, cache, wu, ws, wm, h, theta, motifs)
    if total <= 0: return None, 0
    
    r = random.random() * total
    if n >= s * theta:
        for k in range(theta, n - (s-1)*theta + 1):
            term = countT(k, cache, wu, ws, wm, h, theta, motifs) * countT_stars(s-1, n-k, cache, wu, ws, wm, h, theta, motifs)
            if r < term:
                inner_s, inner_c = generateT(k, cache, wu, ws, wm, h, theta, motifs)
                rest_inners, rest_c = generateT_stars(s-1, n-k, cache, wu, ws, wm, h, theta, motifs)
                if rest_inners is not None:
                    return [inner_s] + rest_inners, inner_c + rest_c
                return None, 0
            r -= term
    return None, 0

def generateT(n, cache, wu, ws, wm, h, theta, motifs):
    """
    Generates a random internal structure dot-bracket string of length n.
    Outputs:
        (str, int) or (None, 0)
    """
    if n < theta: return None, 0
    total = countT(n, cache, wu, ws, wm, h, theta, motifs)
    if total <= 0: return None, 0
    r = random.random() * total
    
    # 1. . S
    if n >= 1:
        term = wu * countS(n-1, cache, wu, ws, wm, h, theta, motifs)
        if r < term:
            rest_s, rest_c = generateS(n-1, cache, wu, ws, wm, h, theta, motifs)
            return "." + rest_s, rest_c
        r -= term
        
    # 2. (*h T )*h . S
    if n >= 2*h + theta + 1:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h - 1 + 1):
            rem = n - (2*h + k + 1)
            term = w_stack_start * countT(k, cache, wu, ws, wm, h, theta, motifs) * wu * countS(rem, cache, wu, ws, wm, h, theta, motifs)
            if r < term:
                inner_s, inner_c = generateT(k, cache, wu, ws, wm, h, theta, motifs)
                rest_s, rest_c = generateS(rem, cache, wu, ws, wm, h, theta, motifs)
                return "(" * h + inner_s + ")" * h + "." + rest_s, inner_c + rest_c
            r -= term

    # 3. (*h T )*h (*h T )*h S
    if n >= 4*h + 2*theta:
        w_stack_start_double = ws**(2 * max(0, h-1))
        for k1 in range(theta, n - 4*h - theta + 1):
            cT_k1 = countT(k1, cache, wu, ws, wm, h, theta, motifs)
            for k2 in range(theta, n - 4*h - k1 + 1):
                rem = n - 4*h - k1 - k2
                term = w_stack_start_double * cT_k1 * countT(k2, cache, wu, ws, wm, h, theta, motifs) * countS(rem, cache, wu, ws, wm, h, theta, motifs)
                if r < term:
                    inner1_s, inner1_c = generateT(k1, cache, wu, ws, wm, h, theta, motifs)
                    inner2_s, inner2_c = generateT(k2, cache, wu, ws, wm, h, theta, motifs)
                    rest_s, rest_c = generateS(rem, cache, wu, ws, wm, h, theta, motifs)
                    return "(" * h + inner1_s + ")" * h + "(" * h + inner2_s + ")" * h + rest_s, inner1_c + inner2_c + rest_c
                r -= term

    # 4. ( T )
    if n >= 2 + theta:
        term = ws * countT(n-2, cache, wu, ws, wm, h, theta, motifs)
        if r < term:
            inner_s, inner_c = generateT(n-2, cache, wu, ws, wm, h, theta, motifs)
            return "(" + inner_s + ")", inner_c
        r -= term
        
    # 5. .* theta
    if n == theta:
        term = wu**theta
        if r < term:
            return "." * theta, 0
        r -= term

    # 6. Dynamic generic motifs
    for motif in motifs:
        stars = motif.count('*')
        L_rem = len(motif) - stars
        if n >= L_rem + stars * theta:
            for k in range(L_rem + stars * theta, n + 1):
                term = wm * countT_stars(stars, k - L_rem, cache, wu, ws, wm, h, theta, motifs) * countS(n - k, cache, wu, ws, wm, h, theta, motifs)
                if r < term:
                    inners, inners_c = generateT_stars(stars, k - L_rem, cache, wu, ws, wm, h, theta, motifs)
                    rest_s, rest_c = generateS(n - k, cache, wu, ws, wm, h, theta, motifs)
                    
                    motif_parts = motif.split('*')
                    motif_str = motif_parts[0]
                    if inners:
                        for idx, inner in enumerate(inners):
                            motif_str += inner + motif_parts[idx+1]
                    
                    return motif_str + rest_s, 1 + inners_c + rest_c
                r -= term
                
    return None, 0

def decompose_helices(ss):
    """
    Decomposes an RNA sequence string and logs its constituent helices computationally.
    """
    p = []
    res,H,C = {},{},{}
    for i,c in enumerate(ss):
        if c=="(":
            p.append(i)
        elif c== ")":
            j = p.pop()
            a,b = j,i
            ii,jj = a+1,b-1
            if (ii,jj) in res:
                hid = res[(ii,jj)]
                res[(a,b)] = hid
                C[hid] += 1
            else:
                hid = 1+len(H)
                res[(a,b)] = hid
                H[hid] = (a,b)
                C[hid] = 1
    return H,C

if __name__ == "__main__":
    import random
    
    # We will generate 50 structures for each motif and combine them into one file
    SELECTED_MOTIFS = [
        "(.((.(*))(*)))",    # 1 star, 2 possibilities
        "(((*).)(*).)"    # 2 stars, 2 possibilities
    ]

    L = 50
    Wu = 1.0
    Ws = 1.0
    Wm = 10.0
    h = 4
    theta = 3
    N_PER_MOTIF = 50
    
    filename = "output/dataset_l50_h4_t3_wu1_motifs40.txt"
    os.makedirs("output", exist_ok=True)
    
    print(f"Starting batch generation of {N_PER_MOTIF} structures per motif into {filename}...")
    
    with open(filename, "w") as f:
        for motif in SELECTED_MOTIFS:
            cache = {}
            # structures will now store (string, count) tuples
            results = []
            print(f"Generating 50 structures for motif: {motif}")
            total_weight = countS(L, cache, Wu, Ws, Wm, h, theta, [motif])
            
            for i in range(N_PER_MOTIF):
                s, c = generateS(L, cache, Wu, Ws, Wm, h, theta, [motif])
                if s:
                    results.append((s, c))
                else:
                    print(f"  Warning: Failed to generate structure {i}")
            
            for s, c in results:
                if c > 0:
                    f.write(f"{s},{c}\n")
                
    print(f"Saved total 100 structures to {filename}")
    print("\nBatch generation complete.")
