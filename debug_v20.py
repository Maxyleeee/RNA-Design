import benchmark_tools
import traceback
import sys

test_structure = "((((....))))"

print("="*60)
print("DEBUGGING NUPACK")
print("="*60)
try:
    seq, ok = benchmark_tools.run_nupack(test_structure)
    print(f"Result: {seq}, {ok}")
except Exception:
    traceback.print_exc()

print("\n" + "="*60)
print("DEBUGGING LEARNA")
print("="*60)
try:
    seq, ok = benchmark_tools.run_learna(test_structure)
    print(f"Result: {seq}, {ok}")
except Exception:
    traceback.print_exc()

print("\n" + "="*60)
print("DEBUGGING DESIRNA")
print("="*60)
try:
    seq, ok = benchmark_tools.run_desirna(test_structure, timeout_sec=20)
    print(f"Result: {seq}, {ok}")
except Exception:
    traceback.print_exc()
