import random
import sys

# Increase recursion depth for deep structures
sys.setrecursionlimit(20000)

def countS(n, cache, wu, ws):
    """
    Counts total weight of structures of length n.
    S -> . S(n-1)
    S -> ( U(k) ) S(n-2-k)
    """
    if n == 0: return 1.0
    if ("S", n) in cache: return cache[("S", n)]
    
    # 1. Unpaired
    val = wu * countS(n-1, cache, wu, ws)
    
    # 2. Paired: (U) S
    # U must be non-empty (k >= 1)
    if n >= 2:
        for k in range(1, n - 1):
            val += ws * countU(k, cache, wu, ws) * countS(n - 2 - k, cache, wu, ws)
            
    cache[("S", n)] = val
    return val

def countU(n, cache, wu, ws):
    """
    Counts total weight of structures "inside pairs" (Grammar U).
    U is identical to S but cannot be empty (n >= 1).
    U -> . S(n-1)
    U -> ( U(k) ) S(n-2-k)
    """
    if n == 0: return 0.0 # U cannot be empty
    if ("U", n) in cache: return cache[("U", n)]
    
    # Same recurrence as S for n >= 1
    val = wu * countS(n-1, cache, wu, ws)
    if n >= 2:
        for k in range(1, n - 1):
            val += ws * countU(k, cache, wu, ws) * countS(n - 2 - k, cache, wu, ws)
            
    cache[("U", n)] = val
    return val

def generateS(n, cache, wu, ws):
    if n == 0: return ""
    
    total = countS(n, cache, wu, ws)
    r = random.random() * total
    
    # 1. Unpaired
    term = wu * countS(n-1, cache, wu, ws)
    if r < term:
        return "." + generateS(n-1, cache, wu, ws)
    r -= term
    
    # 2. Paired
    for k in range(1, n - 1):
        term = ws * countU(k, cache, wu, ws) * countS(n - 2 - k, cache, wu, ws)
        if r < term:
            return "(" + generateU(k, cache, wu, ws) + ")" + generateS(n - 2 - k, cache, wu, ws)
        r -= term
        
    # Fallback (should be rare/impossible with consistent weights)
    return "." + generateS(n-1, cache, wu, ws)

def generateU(n, cache, wu, ws):
    if n == 0: return ""
    
    total = countU(n, cache, wu, ws)
    r = random.random() * total
    
    # 1. Unpaired
    term = wu * countS(n-1, cache, wu, ws)
    if r < term:
        return "." + generateS(n-1, cache, wu, ws)
    r -= term
    
    # 2. Paired
    for k in range(1, n - 1):
        term = ws * countU(k, cache, wu, ws) * countS(n - 2 - k, cache, wu, ws)
        if r < term:
            return "(" + generateU(k, cache, wu, ws) + ")" + generateS(n - 2 - k, cache, wu, ws)
        r -= term
        
    # Fallback
    return "." + generateS(n-1, cache, wu, ws)

if __name__ == "__main__":
    # Define configurations: (Length, Weight Unpaired, Weight Stack, Name)
    configs = [
        {"length": 100, "w_unpaired": 1.0, "w_stack": 1.0, "name": "standard"},
        {"length": 50,  "w_unpaired": 0.5, "w_stack": 5.0, "name": "short_structured"},
        {"length": 200, "w_unpaired": 2.0, "w_stack": 0.5, "name": "long_unstructured"},
        {"length": 100, "w_unpaired": 5.0, "w_stack": 1.0, "name": "unpaired_bias"},
        {"length": 100, "w_unpaired": 1.0, "w_stack": 5.0, "name": "stack_bias"},
    ]
    
    N_STRUCTURES = 20

    print(f"Starting batch generation of {N_STRUCTURES} structures per config...")

    for cfg in configs:
        L = cfg["length"]
        Wu = cfg["w_unpaired"]
        Ws = cfg["w_stack"]
        name = cfg["name"]
        
        # Clear cache for each config because counts depend on weights
        cache = {} 
        structures = []
        
        print(f"Generating set '{name}': Length={L}, Wu={Wu}, Ws={Ws}")
        for i in range(N_STRUCTURES):
            s = generateS(L, cache, Wu, Ws)
            if s:
                structures.append(s)
        
        # Save to output/ folder
        filename = f"output/structures_{name}_L{L}.txt"
        with open(filename, "w") as f:
            for s in structures:
                f.write(s + "\n")
        print(f"Saved {len(structures)} structures to {filename}")

    print("\nBatch generation complete.")
