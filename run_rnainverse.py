import subprocess
import argparse
import sys
import os

def run_rnainverse(structure, timeout_sec=2):
    """
    Runs RNAinverse on a single dot-bracket structure via WSL.
    Returns (sequence, success)
    """
    try:
        process = subprocess.Popen(
            ['wsl', 'RNAinverse'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            stdout, stderr = process.communicate(input=structure + "\n", timeout=timeout_sec)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            # print(f"Timeout ({timeout_sec}s) for structure: {structure}")
            return None, False
        
        if process.returncode != 0:
            print(f"Error running RNAinverse for {structure}: {stderr}")
            return None, False
            
        lines = stdout.strip().split('\n')
        if not lines or not lines[0]:
            return None, False
            
        parts = lines[0].split()
        if len(parts) >= 1:
            sequence = parts[0]
            if len(sequence) == len(structure):
                 return sequence, True
                 
        return None, False
        
    except FileNotFoundError:
        print("Error: 'wsl' or 'RNAinverse' not found. Ensure ViennaRNA is installed in WSL.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None, False

def analyze_file(filepath, max_structures=10, timeout_sec=2):
    """
    Reads a file containing RNA structures (one per line) and runs RNAinverse on each.
    """
    if not os.path.exists(filepath):
        print(f"Error: File '{filepath}' not found.")
        return

    print(f"Analyzing up to {max_structures} structures in: {filepath}")
    
    with open(filepath, 'r') as f:
        # Ignore empty lines and lines that aren't valid dot-bracket
        lines = [line.strip() for line in f if line.strip() and set(line.strip()).issubset({'(', ')', '.'})]
        
    if not lines:
        print("No valid structures found in the file.")
        return
        
    # Limit to max_structures to save time
    lines = lines[:max_structures]
    total = len(lines)
    successes = 0
    results = []
    
    print(f"Found {total} valid structures to test. Running RNAinverse with {timeout_sec}s timeout...")
    
    for i, structure in enumerate(lines):
        print(f"\rProcessing {i+1}/{total} ", end="")
        sys.stdout.flush()
        
        seq, success = run_rnainverse(structure, timeout_sec)
        if success:
            successes += 1
            results.append((structure, seq, True))
        else:
            results.append((structure, None, False))
            
    print("\n\n--- Analysis Results ---")
    print(f"Total Structures: {total}")
    print(f"Designable Structures: {successes} ({successes/total*100:.2f}%)")
    
    # Save results to a file
    out_filepath = filepath.replace('.txt', '_rnainverse_results.txt')
    with open(out_filepath, 'w') as f:
        f.write(f"Total Structures: {total}\n")
        f.write(f"Designable Structures: {successes} ({successes/total*100:.2f}%)\n\n")
        for struct, seq, success in results:
            if success:
                f.write(f"{struct}  =>  {seq}\n")
            else:
                f.write(f"{struct}  =>  FAILED\n")
                
    print(f"Detailed results saved to: {out_filepath}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate RNA structure designability using RNAinverse.")
    parser.add_argument("filepath", help="Path to the file containing RNA structures.")
    args = parser.parse_args()
    
    analyze_file(args.filepath)
