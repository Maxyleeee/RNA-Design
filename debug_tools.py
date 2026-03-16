"""
Diagnostic script to check why NUPACK and LEARNA are failing.
Runs one simple structure and prints ALL output/errors.
"""
import subprocess
import os
import sys
import traceback

# A simple, known-designable structure
test_structure = "((((....))))"

print("=" * 60)
print("DIAGNOSING NUPACK")
print("=" * 60)

# 1) Check if nupack is importable
try:
    import nupack
    print(f"[OK] nupack imported successfully, version: {getattr(nupack, '__version__', 'unknown')}")
except ImportError as e:
    print(f"[FAIL] Cannot import nupack: {e}")
except Exception as e:
    print(f"[FAIL] Error importing nupack: {e}")

# 2) Try to actually run NUPACK design
try:
    import nupack
    config = nupack.Model(material='rna')
    print(f"[OK] nupack.Model created")
    
    # Create a domain with all Ns
    domain = nupack.Domain('N' * len(test_structure), name='d1')
    # Create a strand from the domain
    strand = nupack.TargetStrand([domain], name='s1')
    print(f"[OK] TargetStrand created")
    
    complex_target = nupack.TargetComplex([strand], test_structure, name='c1')
    print(f"[OK] TargetComplex created")
    
    # Create a tube containing the complex
    tube = nupack.TargetTube(on_targets={complex_target: 1e-6}, name='t1')
    design = nupack.Design(tubes=[tube], model=config)
    
    print("[OK] Design object created")
    results = design.run(trials=1)
    if results:
        analysis = results[0].to_analysis()
        # In NUPACK 4, analysis.strands might be a mapping. 
        # Try to get the first sequence found.
        seq = None
        try:
            # If it's a mapping of name -> sequence
            for s_name in analysis.strands:
                seq = analysis.strands[s_name]
                if seq: break
        except:
            try:
                # If it's indexable
                seq = analysis.strands[0]
            except:
                print(f"[DEBUG] contents of analysis.strands: {analysis.strands}")
                raise
        
        if seq:
            print(f"[OK] Designed sequence: {seq}")
        else:
            print("[FAIL] Could not extract sequence from strands")
    else:
        print("[FAIL] results is empty")
except Exception as e:
    print(f"[FAIL] NUPACK design failed: {str(e)}")
    traceback.print_exc()

print()
print("=" * 60)
print("DIAGNOSING LEARNA")
print("=" * 60)

learna_dir = "/home/maxyle/RNA/RNA_design_tools/LEARNA"
learna_bin = os.path.join(learna_dir, "venv", "bin", "learna")
weights_path = os.path.join(learna_dir, "models", "224_0_1")

print(f"learna_bin exists: {os.path.exists(learna_bin)}")
print(f"weights_path exists: {os.path.exists(weights_path)}")

# List the venv/bin to see what's there
venv_bin = os.path.join(learna_dir, "venv", "bin")
if os.path.exists(venv_bin):
    print(f"Contents of {venv_bin}:")
    for f in sorted(os.listdir(venv_bin)):
        print(f"  {f}")
else:
    print(f"[FAIL] venv/bin does not exist at {venv_bin}")
    # Try to find learna
    if os.path.exists(learna_dir):
        print(f"Contents of {learna_dir}:")
        for f in sorted(os.listdir(learna_dir)):
            print(f"  {f}")

# Try running LEARNA
cmd = [learna_bin, "--target_structure", test_structure, "--timeout", "10"]
print(f"\nRunning: {' '.join(cmd)}")
try:
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
    print(f"Return code: {res.returncode}")
    print(f"STDOUT:\n{res.stdout}")
    print(f"STDERR:\n{res.stderr}")
except FileNotFoundError:
    print(f"[FAIL] learna binary not found at {learna_bin}")
except subprocess.TimeoutExpired:
    print(f"[FAIL] LEARNA timed out after 20s")
except Exception as e:
    print(f"[FAIL] Error running LEARNA:")
    traceback.print_exc()
