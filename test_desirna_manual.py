import benchmark_tools
import sys

struct = "((.(((((((...)))(((...)))..((.((((....))))))))))))"
print(f"Testing DesiRNA on L=50 struct: {struct}")
try:
    res = benchmark_tools.run_desirna(struct, timeout_sec=60)
    print(f"DesiRNA result: {res}")
except Exception as e:
    import traceback
    traceback.print_exc()
