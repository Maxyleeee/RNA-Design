import random
import sys

# Increase recursion depth for deep structures
sys.setrecursionlimit(10000)

def countS(n,cache,weight_unpaired,weight_stack):
    if ("S",n) not in cache:
        if n == 0:
            val = 1
        else:
            val = weight_unpaired*countS(n-1,cache,weight_unpaired,weight_stack)
            if n>1:
                val += countS(n-2,cache,weight_unpaired,weight_stack)
            if n>2:
                for i in range(1,n-1):
                    val += countT(i,cache,weight_unpaired,weight_stack)*countS(n-2-i,cache,weight_unpaired,weight_stack)
        cache[("S",n)] = val
    return cache[("S",n)]
    
def generateS(n,cache,weight_unpaired,weight_stack):
    if n == 0:
        return ""
    else:
        r = random.random() * countS(n,cache,weight_unpaired,weight_stack)
        r -= weight_unpaired * countS(n-1,cache,weight_unpaired,weight_stack)
        if r<0:
            return "." + generateS(n-1,cache,weight_unpaired,weight_stack)
        if n>1:
            r -= countS(n-2,cache,weight_unpaired,weight_stack)
            if r<0:
                return "()" + generateS(n-2,cache,weight_unpaired,weight_stack)
        if n>2:
            for i in range(1,n-1):
                r -= countT(i,cache,weight_unpaired,weight_stack)*countS(n-2-i,cache,weight_unpaired,weight_stack)
                if r<0:
                    return "("+generateT(i,cache,weight_unpaired,weight_stack)+")"+generateS(n-2-i,cache,weight_unpaired,weight_stack)
    return None

def countT(n,cache,weight_unpaired,weight_stack):
    if ("T",n) not in cache:
        if n == 0:
            val = 0
        else:
            val = weight_unpaired*countS(n-1,cache,weight_unpaired,weight_stack)
            if n==2:
                val += weight_stack
            if n>2:
                val += countT(n-2,cache,weight_unpaired,weight_stack)
            if n>2:
                val += weight_stack*countT(n-2,cache,weight_unpaired,weight_stack)
            if n>3:
                for i in range(1,n-2):
                    val += countT(i,cache,weight_unpaired,weight_stack)*countT(n-2-i,cache,weight_unpaired,weight_stack)
        cache[("T",n)] = val
    return cache[("T",n)]
    
def generateT(n,cache,weight_unpaired,weight_stack):
    if n == 0:
        return "#"
    else:
        r = random.random()*countT(n,cache,weight_unpaired,weight_stack)
        r -= weight_unpaired*countS(n-1,cache,weight_unpaired,weight_stack)
        if r<0:
            return "."+generateS(n-1,cache,weight_unpaired,weight_stack)
        if n==2:
            r -= weight_stack
            if r<0:
                return "()"
        if n>2:
            r -= countT(n-2,cache,weight_unpaired,weight_stack)
            if r<0:
                return "()"+generateT(n-2,cache,weight_unpaired,weight_stack)
        if n>2:
            r -= weight_stack*countT(n-2,cache,weight_unpaired,weight_stack)
            if r<0:
                return "("+generateT(n-2,cache,weight_unpaired,weight_stack)+")"
        if n>3:
            for i in range(1,n-2):
                r -= countT(i,cache,weight_unpaired,weight_stack)*countT(n-2-i,cache,weight_unpaired,weight_stack)
                if r<0:
                    return "("+generateT(i,cache,weight_unpaired,weight_stack)+")"+generateT(n-2-i,cache,weight_unpaired,weight_stack)
    return None

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
        
        filename = f"structures_{name}_L{L}.txt"
        with open(filename, "w") as f:
            for s in structures:
                f.write(s + "\n")
        print(f"Saved {len(structures)} structures to {filename}")

    print("\\nBatch generation complete.")
