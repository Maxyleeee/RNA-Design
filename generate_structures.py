import random
import sys

# Increase recursion depth for deep structures
sys.setrecursionlimit(20000)

def countS(n, cache, wu, ws, h, theta):
    """
    Counts total weight of structures of length n starting with S.
    S -> . S(n-1)
    S -> (*h T(k) )*h S(n-2h-k)
    S -> eps
    """
    if n == 0: return 1.0
    if ("S", n) in cache: return cache[("S", n)]
    
    val = 0.0
    
    # 1. Unpaired: . S
    if n >= 1:
        val += wu * countS(n-1, cache, wu, ws, h, theta)
    
    # 2. Paired: (*h T )*h S
    # Requires 2*h + k + m = n
    # Helix uses 2*h. T uses k. S uses m.
    # T requires min length theta. So k >= theta.
    # m >= 0.
    # max k = n - 2*h.
    if n >= 2*h + theta:
        for k in range(theta, n - 2*h + 1):
            term = ws * countT(k, cache, wu, ws, h, theta) * countS(n - 2*h - k, cache, wu, ws, h, theta)
            val += term
            
    cache[("S", n)] = val
    return val

def countT(n, cache, wu, ws, h, theta):
    """
    Counts total weight of structures of length n inside a helix (Grammar T).
    Constraint: n >= theta (min pair distance).
    
    T -> . S
    T -> (*h T )*h . S
    T -> (*h T )*h (*h T )*h S
    """
    if n < theta: return 0.0
    if ("T", n) in cache: return cache[("T", n)]
    
    val = 0.0
    
    # 1. T -> . S
    # Requires 1 + m = n => m = n - 1
    # n >= theta is guaranteed by start check.
    # S(n-1) valid if n-1 >= 0. Since n>=theta>=0, n-1 could be -1 if theta=0? 
    # Usually theta >= 1. If theta=3, n>=3, n-1>=2.
    if n >= 1:
        val += wu * countS(n-1, cache, wu, ws, h, theta)
        
    # 2. T -> (*h T )*h . S
    # Structure: Helix(size 2h+k) . S(size m)
    # Total n = (2h + k) + 1 + m = 2h + k + 1 + m
    # k >= theta. m >= 0.
    # Min length required: 2h + theta + 1.
    if n >= 2*h + theta + 1:
        # Loop over possible k (size of inner T)
        # max k: n - 2h - 1 (since m=0).
        for k in range(theta, n - 2*h - 1 + 1):
            rem_after_helix_dot = n - (2*h + k + 1)
            term = ws * countT(k, cache, wu, ws, h, theta) * wu * countS(rem_after_helix_dot, cache, wu, ws, h, theta)
            val += term

    # 3. T -> (*h T )*h (*h T )*h S
    # Structure: H1(size 2h+k1) H2(size 2h+k2) S(size m)
    # Total n = (2h+k1) + (2h+k2) + m = 4h + k1 + k2 + m
    # k1, k2 >= theta. m >= 0.
    # Min length: 4h + 2*theta.
    if n >= 4*h + 2*theta:
        # We can optimize this sum if needed, but for length ~100-200, nested loops are acceptable.
        # Loop k1
        # max k1: n - (2h + theta) - (rem_min=0) - 2h = n - 4h - theta.
        for k1 in range(theta, n - 4*h - theta + 1):
            min_k2 = theta
            max_k2 = n - 4*h - k1 # since m >= 0
            
            # Precompute term1
            term1 = ws * countT(k1, cache, wu, ws, h, theta)
            
            for k2 in range(min_k2, max_k2 + 1):
                rem_after_H1_H2 = n - 4*h - k1 - k2
                term = term1 * ws * countT(k2, cache, wu, ws, h, theta) * countS(rem_after_H1_H2, cache, wu, ws, h, theta)
                val += term

    cache[("T", n)] = val
    return val

def generateS(n, cache, wu, ws, h, theta):
    if n == 0: return ""
    
    total = countS(n, cache, wu, ws, h, theta)
    if total == 0: return None # Should not happen for S unless constraints impossible
    r = random.random() * total
    
    # 1. . S
    if n >= 1:
        term = wu * countS(n-1, cache, wu, ws, h, theta)
        if r < term:
            rest = generateS(n-1, cache, wu, ws, h, theta)
            return "." + rest if rest is not None else "."
        r -= term
        
    # 2. (*h T )*h S
    if n >= 2*h + theta:
        for k in range(theta, n - 2*h + 1):
            term = ws * countT(k, cache, wu, ws, h, theta) * countS(n - 2*h - k, cache, wu, ws, h, theta)
            if r < term:
                inner = generateT(k, cache, wu, ws, h, theta)
                rest = generateS(n - 2*h - k, cache, wu, ws, h, theta)
                # Construct (*h ... )*h
                return "(" * h + inner + ")" * h + rest
            r -= term
            
    # Fallback/Rounding
    # If we are here, it means floating point issues or 0 total.
    # Return unpaired as default if possible
    if n >= 1: return "." + generateS(n-1, cache, wu, ws, h, theta)
    return ""

def generateT(n, cache, wu, ws, h, theta):
    # n >= theta check done by caller usually, but logic holds.
    # Note: T is strictly length n.
    if n < theta: return None 
    
    total = countT(n, cache, wu, ws, h, theta)
    if total == 0: return None
    r = random.random() * total
    
    # 1. . S
    if n >= 1:
        term = wu * countS(n-1, cache, wu, ws, h, theta)
        if r < term:
            rest = generateS(n-1, cache, wu, ws, h, theta)
            return "." + rest
        r -= term
        
    # 2. (*h T )*h . S
    if n >= 2*h + theta + 1:
        for k in range(theta, n - 2*h - 1 + 1):
            rem = n - (2*h + k + 1)
            term = ws * countT(k, cache, wu, ws, h, theta) * wu * countS(rem, cache, wu, ws, h, theta)
            if r < term:
                inner = generateT(k, cache, wu, ws, h, theta)
                rest = generateS(rem, cache, wu, ws, h, theta)
                return "(" * h + inner + ")" * h + "." + rest
            r -= term

    # 3. (*h T )*h (*h T )*h S
    if n >= 4*h + 2*theta:
        for k1 in range(theta, n - 4*h - theta + 1):
            term1_val = ws * countT(k1, cache, wu, ws, h, theta)
            
            for k2 in range(theta, n - 4*h - k1 + 1):
                rem = n - 4*h - k1 - k2
                term = term1_val * ws * countT(k2, cache, wu, ws, h, theta) * countS(rem, cache, wu, ws, h, theta)
                if r < term:
                    inner1 = generateT(k1, cache, wu, ws, h, theta)
                    inner2 = generateT(k2, cache, wu, ws, h, theta)
                    rest = generateS(rem, cache, wu, ws, h, theta)
                    return "(" * h + inner1 + ")" * h + "(" * h + inner2 + ")" * h + rest
                r -= term
                
    return None # Should not happen if total > 0

def decompose_helices(ss):
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
    # Define configurations: (Length, Weight Unpaired, Weight Stack, Name, Min Helix h, Min Loop theta)
    # Defaulting h=3, theta=3 for existing configs as a reasonable starting point
    configs = [
        {"length": 100, "w_unpaired": 1.0, "w_stack": 1.0, "name": "standard", "h": 3, "theta": 3},
        {"length": 50,  "w_unpaired": 0.5, "w_stack": 5.0, "name": "short_structured", "h": 2, "theta": 3},
        {"length": 200, "w_unpaired": 2.0, "w_stack": 0.5, "name": "long_unstructured", "h": 3, "theta": 3},
        {"length": 100, "w_unpaired": 5.0, "w_stack": 1.0, "name": "unpaired_bias", "h": 3, "theta": 3},
        {"length": 100, "w_unpaired": 1.0, "w_stack": 5.0, "name": "stack_bias", "h": 4, "theta": 3},
    ]
    
    N_STRUCTURES = 20

    print(f"Starting batch generation of {N_STRUCTURES} structures per config...")

    for cfg in configs:
        L = cfg["length"]
        Wu = cfg["w_unpaired"]
        Ws = cfg["w_stack"]
        h = cfg.get("h", 3)
        theta = cfg.get("theta", 3)
        name = cfg["name"]
        
        # Clear cache for each config because counts depend on weights and params
        cache = {} 
        structures = []
        
        print(f"Generating set '{name}': Length={L}, Wu={Wu}, Ws={Ws}, h={h}, theta={theta}")
        
        # Precompute counts
        total_weight = countS(L, cache, Wu, Ws, h, theta)
        print(f"  Total weight for L={L}: {total_weight:.2e}")
        
        for i in range(N_STRUCTURES):
            s = generateS(L, cache, Wu, Ws, h, theta)
            if s:
                structures.append(s)
            else:
                print(f"  Warning: Failed to generate structure {i}")
        
        # Save to output/ folder
        filename = f"output/structures_{name}_L{L}.txt"
        # Ensure directory exists
        import os
        os.makedirs("output", exist_ok=True)
        
        with open(filename, "w") as f:
            for s in structures:
                f.write(s + "\n")
        print(f"Saved {len(structures)} structures to {filename}")

        # Decompose and print helices for each structure
        for ss in structures:
            H,C = decompose_helices(ss)
            print(ss)
            for hid in H.keys():
                i,j = H[hid]
                nbbps = C[hid]
                k,l = i+nbbps-1,j-(nbbps-1)
                print("  Helix %s:"%hid,(i,j),"->",(k,l),"%s BPs"%(nbbps))


    print("\nBatch generation complete.")


