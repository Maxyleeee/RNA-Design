"""
Generate per-tool success/failure classification files from benchmark .txt results.
Each file lists which structures each tool succeeded and failed on.
"""
import os
import re

BENCHMARK_FILES = [
    "structures_motif_h4_benchmark.txt",
    "structures_motif_h5_benchmark.txt",
]

TOOLS = ["RNAinverse", "NEMO", "eM2dRNAs", "DesiRNA", "LEARNA", "NUPACK"]

def parse_benchmark(filepath):
    """Parse a benchmark output file and return per-structure, per-tool results."""
    results = []
    current = None
    current_tool = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            # Detect target structure line
            m = re.match(r"^Target:\s*(.+)$", line)
            if m:
                if current:
                    results.append(current)
                current = {"target": m.group(1), "tools": {}}
                current_tool = None
                continue
            
            if current is None:
                continue
                
            # Detect tool name lines (indented with 2 spaces)
            for tool in TOOLS:
                if re.match(rf"^  {re.escape(tool)}:", line):
                    current_tool = tool
                    # Check if line itself says FAILED (single-line format)
                    if "FAILED" in line:
                        current["tools"][tool] = "FAILED"
                        current_tool = None
                    break
            
            # Check Match: YES/NO on indented sub-lines
            if current_tool and re.match(r"^    Match:\s", line):
                if "YES" in line:
                    current["tools"][current_tool] = "SUCCESS"
                else:
                    current["tools"][current_tool] = "NO_MATCH"
                current_tool = None

    if current:
        results.append(current)
    return results


def write_tool_classification(results, config_name):
    """Write per-tool classification files: successes and failures."""
    tool_success = {t: [] for t in TOOLS}
    tool_fail = {t: [] for t in TOOLS}

    for r in results:
        for tool in TOOLS:
            status = r["tools"].get(tool, "FAILED")
            if status == "SUCCESS":
                tool_success[tool].append(r["target"])
            else:
                tool_fail[tool].append(r["target"])

    for tool in TOOLS:
        tool_safe = tool.replace("/", "_")
        fname = f"classification_{config_name}_{tool_safe}.txt"
        total = len(results)
        with open(fname, 'w') as f:
            f.write(f"=== TOOL: {tool} | Config: {config_name} ===\n")
            f.write(f"Successes: {len(tool_success[tool])}/{total}\n")
            f.write(f"Failures:  {len(tool_fail[tool])}/{total}\n\n")
            
            f.write(f"--- SUCCEEDED ({len(tool_success[tool])}) ---\n")
            for s in tool_success[tool]:
                f.write(s + "\n")
            
            f.write(f"\n--- FAILED ({len(tool_fail[tool])}) ---\n")
            for s in tool_fail[tool]:
                f.write(s + "\n")
        
        print(f"  Written: {fname}")


def write_summary_file(all_results, config_name):
    """Write a combined summary file showing which tools solved each structure."""
    fname = f"classification_{config_name}_summary.txt"
    with open(fname, 'w') as f:
        f.write(f"=== SUMMARY: {config_name} ({len(all_results)} structures) ===\n\n")
        f.write(f"{'Structure':100s}  " + "  ".join(f"{t:8s}" for t in TOOLS) + "\n")
        f.write("-" * (100 + 3 * len(TOOLS) + 10 * len(TOOLS)) + "\n")
        
        for r in all_results:
            struct_abbrev = r["target"][:97] + "..." if len(r["target"]) > 100 else r["target"]
            row = f"{struct_abbrev:100s}"
            for tool in TOOLS:
                status = r["tools"].get(tool, "FAILED")
                icon = "YES" if status == "SUCCESS" else ("---" if status == "FAILED" else "NO ")
                row += f"  {icon:8s}"
            f.write(row + "\n")
    print(f"  Written: {fname}")


def write_per_structure_file(all_results, config_name):
    """Write a tracking file with every structure and all tool results."""
    fname = f"classification_{config_name}_per_structure.txt"
    with open(fname, 'w') as f:
        f.write(f"Per-Structure Tool Results: {config_name}\n")
        f.write("=" * 80 + "\n\n")
        for i, r in enumerate(all_results, 1):
            f.write(f"[{i:03d}] {r['target']}\n")
            for tool in TOOLS:
                status = r["tools"].get(tool, "FAILED")
                icon = "✓ SOLVED" if status == "SUCCESS" else ("✗ FAILED" if status == "FAILED" else "~ NO_MFE_MATCH")
                f.write(f"  {tool:12s}: {icon}\n")
            f.write("\n")
    print(f"  Written: {fname}")


def main():
    for bf in BENCHMARK_FILES:
        if not os.path.exists(bf):
            print(f"Benchmark file not found: {bf}, skipping.")
            continue
        
        config = bf.replace("structures_", "").replace("_benchmark.txt", "")
        print(f"\nProcessing {bf} (config: {config})")
        
        results = parse_benchmark(bf)
        print(f"  Parsed {len(results)} structures")
        
        write_tool_classification(results, config)
        write_summary_file(results, config)
        write_per_structure_file(results, config)
    
    print("\nDone! All classification files written.")


if __name__ == "__main__":
    main()
