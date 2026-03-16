import benchmark_tools
import sys

struct = "((.(((((((...)))(((...)))..((.((((....))))))))))))"
print(f"Testing NUPACK debug on: {struct}")
try:
    res = benchmark_tools.run_nupack(struct, timeout_sec=20)
    print(f"NUPACK result: {res}")
except Exception as e:
    import traceback
    traceback.print_exc()
