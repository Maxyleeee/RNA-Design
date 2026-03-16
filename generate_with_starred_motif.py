import random
import sys
import os

# Increase recursion depth for deep structures
sys.setrecursionlimit(20000)

def countS(n, cache, wu, ws, wm, h, theta):
    """
    Counts total computational weight for structures starting with grammar S.
    Inputs:
        n: int - Target sequence length
        cache: dict - Memoization cache
        wu: float - Unpaired base weight
        ws: float - Stacked base pair weight
        wm: float - Fixed starred motif weight
        h: int - Min helix threshold
        theta: int - Min loop length
    Outputs:
        float - Total path weight count
    S -> . S
    S -> (*h T )*h S  (Weight ws**h)
    S -> eps
    """
    if n == 0: return 1.0
    if ("S", n) in cache: return cache[("S", n)]
    
    val = 0.0
    
    # 1. . S
    if n >= 1:
        val += wu * countS(n-1, cache, wu, ws, wm, h, theta)
    
    # 2. (*h T )*h S
    if n >= 2*h + theta:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h + 1):
            val += w_stack_start * countT(k, cache, wu, ws, wm, h, theta) * countS(n - 2*h - k, cache, wu, ws, wm, h, theta)
            
    cache[("S", n)] = val
    return val

def countT(n, cache, wu, ws, wm, h, theta):
    """
    Counts internal loops (Grammar T) computation weight.
    Inputs:
        n: int - Length of internal sub sequence
        cache: dict - Memoization cache
        wu: float - Unpaired base weight
        ws: float - Stacked base pair weight
        wm: float - Fixed starred motif weight
        h: int - Min helix length
        theta: int - Min loop distance
    Outputs:
        float - Derived weight for T nodes
    T -> . S
    T -> (*h T )*h . S  (Weight ws**h)
    T -> (*h T )*h (*h T )*h S (Weight ws**2h)
    T -> ( T )  (Weight ws)
    T -> .* theta (Weight wu**theta)
    T -> Motif S      where Motif is (( T ).(....)) (Weight wm)
    """
    if n < theta: return 0.0
    if ("T", n) in cache: return cache[("T", n)]
    
    val = 0.0
    
    # 1. T -> . S
    if n >= 1:
        val += wu * countS(n-1, cache, wu, ws, wm, h, theta)
        
    # 2. T -> (*h T )*h . S
    if n >= 2*h + theta + 1:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h - 1 + 1):
            rem = n - (2*h + k + 1)
            val += w_stack_start * countT(k, cache, wu, ws, wm, h, theta) * wu * countS(rem, cache, wu, ws, wm, h, theta)

    # 3. T -> (*h T )*h (*h T )*h S
    if n >= 4*h + 2*theta:
        w_stack_start_double = ws**(2 * max(0, h-1))
        for k1 in range(theta, n - 4*h - theta + 1):
            cT_k1 = countT(k1, cache, wu, ws, wm, h, theta)
            for k2 in range(theta, n - 4*h - k1 + 1):
                rem = n - 4*h - k1 - k2
                val += w_stack_start_double * cT_k1 * countT(k2, cache, wu, ws, wm, h, theta) * countS(rem, cache, wu, ws, wm, h, theta)

    # 4. T -> ( T )
    if n >= 2 + theta:
        val += ws * countT(n-2, cache, wu, ws, wm, h, theta)
        
    # 5. T -> .* theta
    if n == theta:
        val += wu**theta

    # 6. T -> Motif S
    if n >= 11 + theta:
        for k in range(theta, n - 11 + 1):
            val += wm * countT(k, cache, wu, ws, wm, h, theta) * countS(n - 11 - k, cache, wu, ws, wm, h, theta)

    cache[("T", n)] = val
    return val

def generateS(n, cache, wu, ws, wm, h, theta):
    """
    Stochastically generates a master RNA sequence string using grammar parameters.
    Inputs:
        n: int - String length
        cache: dict - Precalculated weights dictionary
        wu: float - Unpaired weight configuration
        ws: float - Base pair stacking configuration
        wm: float - Starred motif fixed constraint weight
        h: int - Helix setting limits
        theta: int - Hairpin gap configuration
    Outputs:
        str - Master dot bracket format sequence layout
    """
    if n == 0: return ""
    total = countS(n, cache, wu, ws, wm, h, theta)
    if total <= 0: return "." * n 
    r = random.random() * total
    
    # 1. . S
    if n >= 1:
        term = wu * countS(n-1, cache, wu, ws, wm, h, theta)
        if r < term:
            return "." + generateS(n-1, cache, wu, ws, wm, h, theta)
        r -= term
        
    # 2. (*h T )*h S
    if n >= 2*h + theta:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h + 1):
            term = w_stack_start * countT(k, cache, wu, ws, wm, h, theta) * countS(n - 2*h - k, cache, wu, ws, wm, h, theta)
            if r < term:
                inner = generateT(k, cache, wu, ws, wm, h, theta)
                rest = generateS(n - 2*h - k, cache, wu, ws, wm, h, theta)
                return "(" * h + inner + ")" * h + rest
            r -= term
            
    return "." + generateS(n-1, cache, wu, ws, wm, h, theta)

def generateT(n, cache, wu, ws, wm, h, theta):
    """
    Stochastically plots an internal structural variation constrained within parent loops.
    Inputs:
        n: int - Sequence snippet length
        cache: dict - System weight registry
        wu: float - Unpaired configuration scale
        ws: float - Secondary pairs scale measure
        wm: float - Fixed starred motif weight
        h: int - Strict min helix bound
        theta: int - Strict min loop stretch 
    Outputs:
        str or None - Processed dot-bracket sub component
    """
    if n < theta: return None 
    total = countT(n, cache, wu, ws, wm, h, theta)
    if total <= 0: return None
    r = random.random() * total
    
    # 1. . S
    if n >= 1:
        term = wu * countS(n-1, cache, wu, ws, wm, h, theta)
        if r < term:
            return "." + generateS(n-1, cache, wu, ws, wm, h, theta)
        r -= term
        
    # 2. (*h T )*h . S
    if n >= 2*h + theta + 1:
        w_stack_start = ws**max(0, h-1)
        for k in range(theta, n - 2*h - 1 + 1):
            rem = n - (2*h + k + 1)
            term = w_stack_start * countT(k, cache, wu, ws, wm, h, theta) * wu * countS(rem, cache, wu, ws, wm, h, theta)
            if r < term:
                inner = generateT(k, cache, wu, ws, wm, h, theta)
                rest = generateS(rem, cache, wu, ws, wm, h, theta)
                return "(" * h + inner + ")" * h + "." + rest
            r -= term

    # 3. (*h T )*h (*h T )*h S
    if n >= 4*h + 2*theta:
        w_stack_start_double = ws**(2 * max(0, h-1))
        for k1 in range(theta, n - 4*h - theta + 1):
            cT_k1 = countT(k1, cache, wu, ws, wm, h, theta)
            for k2 in range(theta, n - 4*h - k1 + 1):
                rem = n - 4*h - k1 - k2
                term = w_stack_start_double * cT_k1 * countT(k2, cache, wu, ws, wm, h, theta) * countS(rem, cache, wu, ws, wm, h, theta)
                if r < term:
                    inner1 = generateT(k1, cache, wu, ws, wm, h, theta)
                    inner2 = generateT(k2, cache, wu, ws, wm, h, theta)
                    rest = generateS(rem, cache, wu, ws, wm, h, theta)
                    return "(" * h + inner1 + ")" * h + "(" * h + inner2 + ")" * h + rest
                r -= term

    # 4. ( T )
    if n >= 2 + theta:
        term = ws * countT(n-2, cache, wu, ws, wm, h, theta)
        if r < term:
            return "(" + generateT(n-2, cache, wu, ws, wm, h, theta) + ")"
        r -= term
        
    # 5. .* theta
    if n == theta:
        return "." * theta
        
    # 6. T -> Motif S
    if n >= 11 + theta:
        for k in range(theta, n - 11 + 1):
            term = wm * countT(k, cache, wu, ws, wm, h, theta) * countS(n - 11 - k, cache, wu, ws, wm, h, theta)
            if r < term:
                inner_T = generateT(k, cache, wu, ws, wm, h, theta)
                rest_S = generateS(n - 11 - k, cache, wu, ws, wm, h, theta)
                motif_str = "((" + inner_T + ").(....))"
                return motif_str + rest_S
            r -= term
                
    return None

def decompose_helices(ss):
    """
    Programmatically logs coordinates representing topological RNA shapes.
    Inputs:
        ss: str - Secondary sequence geometry notation
    Outputs:
        tuple (dict, dict) - Identified helix indices mapping and total length map
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
    # Define configurations: (Length, Weight Unpaired, Weight Stack, Weight Motif, Name, Min Helix h, Min Loop theta)
    configs = [
        {"length": 100, "w_unpaired": 1.0, "w_stack": 2.0, "w_motif": 5.0, "name": "starred_motif_bias_strong_stacks", "h": 3, "theta": 3},
        {"length": 100, "w_unpaired": 1.0, "w_stack": 0.5, "w_motif": 5.0, "name": "starred_motif_bias_weak_stacks", "h": 3, "theta": 3},
        {"length": 100, "w_unpaired": 1.0, "w_stack": 1.2, "w_motif": 1.0, "name": "starred_balanced_motif", "h": 3, "theta": 3},
    ]
    
    N_STRUCTURES = 20

    print(f"Starting batch generation of {N_STRUCTURES} structures per config...")

    for cfg in configs:
        L = cfg["length"]
        Wu = cfg["w_unpaired"]
        Ws = cfg["w_stack"]
        Wm = cfg["w_motif"]
        h = cfg.get("h", 3)
        theta = cfg.get("theta", 3)
        name = cfg["name"]
        
        # Clear cache for each config because counts depend on weights and params
        cache = {} 
        structures = []
        
        print(f"Generating set '{name}': Length={L}, Wu={Wu}, Ws={Ws}, Wm={Wm}, h={h}, theta={theta}")
        
        # Precompute counts
        total_weight = countS(L, cache, Wu, Ws, Wm, h, theta)
        print(f"  Total weight for L={L}: {total_weight:.2e}")
        
        for i in range(N_STRUCTURES):
            s = generateS(L, cache, Wu, Ws, Wm, h, theta)
            if s:
                structures.append(s)
            else:
                print(f"  Warning: Failed to generate structure {i}")
        
        # Save to output/ folder
        filename = f"output/structures_with_{name}.txt"
        os.makedirs("output", exist_ok=True)
        
        # Write config metadata string + the structures
        with open(filename, "w") as f:
            f.write(f"# Generation settings: Length={L}, Wu(unpaired)={Wu}, Ws(stacks)={Ws}, Wm(motif)={Wm}, h(min helix)={h}, theta(min loop)={theta}\n")
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
