import os
import re

def parse_benchmark_blocks(filepath):
    """Parses a benchmark file into a list of (target_structure, block_text) tuples."""
    if not os.path.exists(filepath):
        print(f"INFO: File not found: {filepath}")
        return []
        
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Blocks are separated by 60 hyphens
    raw_blocks = content.split("-" * 60)
    parsed = []
    
    for rb in raw_blocks:
        rb = rb.strip()
        if "Target:" in rb:
            # Strip any leading header text from the first target block
            idx = rb.find("Target:")
            clean_block = rb[idx:].strip()
            
            match = re.search(r"Target: ([.()]+)", clean_block)
            if match:
                struct = match.group(1)
                parsed.append((struct, clean_block))
    
    return parsed

def is_motif_less(block_text):
    """Checks if a block corresponds to a motif-less structure (Baseline Match: YES)."""
    if "Baseline:" not in block_text:
        return False
    
    baseline_section = block_text.split("Baseline:")[1]
    if "RNAinverse:" in baseline_section:
        baseline_section = baseline_section.split("RNAinverse:")[0]
    
    return "Match: YES" in baseline_section

def generate_summary(blocks):
    """Parses stats from details blocks and builds a summary table."""
    TOOLS = ["Baseline", "RNAinverse", "NEMO", "eM2dRNAs", "DesiRNA", "LEARNA", "NUPACK"]
    stats = {t: {"designed": 0, "verified": 0, "times_success": [], "times_failure": [], "total_time": 0.0} for t in TOOLS}
    
    for _, block in blocks:
        for tool in TOOLS:
            # Extract tool section
            # Pattern looks for 'ToolName:' followed by indented lines until the next 'ToolName:' or separator
            pattern = rf"  {tool}:(.*?)(?=\n  [A-Z]|\n-|\Z)"
            match = re.search(pattern, block, re.DOTALL)
            
            if match:
                section = match.group(1)
                # If there's a Seq, it was designed
                if "Seq:" in section:
                    stats[tool]["designed"] += 1
                    
                    # Check for match and time
                    is_match = "Match: YES" in section
                    time_match = re.search(r"Time: ([\d.]+)s", section)
                    t_val = float(time_match.group(1)) if time_match else 0.0
                    
                    stats[tool]["total_time"] += t_val
                    if is_match:
                        stats[tool]["verified"] += 1
                        stats[tool]["times_success"].append(t_val)
                    else:
                        stats[tool]["times_failure"].append(t_val)
                elif "FAILED" in section:
                    # Tool execution failed entirely
                    time_match = re.search(r"Time: ([\d.]+)s", section)
                    t_val = float(time_match.group(1)) if time_match else 0.0
                    stats[tool]["total_time"] += t_val

    total = len(blocks)
    summary = "SUMMARY (Generated from 100 structures)\n" + "=" * 105 + "\n"
    header = f"{'Tool':<14} {'Designed':>10} {'Verified':>10} {'Match Rate':>12} {'Avg T(S)':>10} {'Avg T(F)':>10} {'Avg T(A)':>10}"
    summary += header + "\n" + "-" * len(header) + "\n"
    
    for t in TOOLS:
        d = stats[t]["designed"]
        v = stats[t]["verified"]
        rate = f"{v}/{total} ({v/total*100:.0f}%)" if total > 0 else "N/A"
        avg_t_s = sum(stats[t]["times_success"]) / len(stats[t]["times_success"]) if stats[t]["times_success"] else 0
        avg_t_f = sum(stats[t]["times_failure"]) / len(stats[t]["times_failure"]) if stats[t]["times_failure"] else 0
        avg_t_a = stats[t]["total_time"] / total if total > 0 else 0
        summary += f"{t:<14} {d:>5}/{total:<4} {v:>5}/{total:<4} {rate:>12} {avg_t_s:>9.1f}s {avg_t_f:>9.1f}s {avg_t_a:>9.1f}s\n"
    
    return summary

def merge():
    base_path = "output"
    if not os.path.exists(base_path):
        base_path = "."
    
    print(f"Merging results from directory: {os.path.abspath(base_path)}")
    
    # === Dataset 1 ===
    f1_old = os.path.join(base_path, "dataset_l50_h4_t3_wu1_2motifs_benchmark.txt")
    f1_new = os.path.join(base_path, "dataset_l50_h4_t3_wu1_2motifs_part2_benchmark.txt")
    
    blocks1_old = parse_benchmark_blocks(f1_old)
    valid_old1 = blocks1_old[:70] 
    blocks1_new = parse_benchmark_blocks(f1_new)
    
    final_blocks1 = valid_old1 + blocks1_new
    print(f"Dataset 1: Found {len(valid_old1)} old + {len(blocks1_new)} new = {len(final_blocks1)} total structures.")
    
    if len(final_blocks1) > 0:
        summary1 = generate_summary(final_blocks1)
        out1 = os.path.join(base_path, "dataset_l50_h4_t3_wu1_2motifs_final_100_details.txt")
        with open(out1, 'w', encoding='utf-8') as f:
            f.write(summary1 + "\n\n")
            f.write("DETAILS\n" + "=" * 80 + "\n\n")
            for _, b in final_blocks1:
                f.write(b + "\n" + "-" * 60 + "\n")
        print(f"Successfully generated: {out1}")
    
    # === Dataset 2 ===
    f2_old = os.path.join(base_path, "dataset_l50_h4_t3_wu1_2motifs2_benchmark.txt")
    f2_new = os.path.join(base_path, "dataset_l50_h4_t3_wu1_2motifs2_part2_benchmark.txt")
    
    blocks2_old = parse_benchmark_blocks(f2_old)
    # Correctly filter motif-embedded ones
    valid_old2 = [ (s,b) for s,b in blocks2_old if not is_motif_less(b) ]
    valid_old2 = valid_old2[:86]
    
    blocks2_new = parse_benchmark_blocks(f2_new)
    final_blocks2 = valid_old2 + blocks2_new
    print(f"Dataset 2: Found {len(valid_old2)} old (motif-embedded) + {len(blocks2_new)} new = {len(final_blocks2)} total structures.")

    if len(final_blocks2) > 0:
        summary2 = generate_summary(final_blocks2)
        out2 = os.path.join(base_path, "dataset_l50_h4_t3_wu1_2motifs2_final_100_details.txt")
        with open(out2, 'w', encoding='utf-8') as f:
            f.write(summary2 + "\n\n")
            f.write("DETAILS\n" + "=" * 80 + "\n\n")
            for _, b in final_blocks2:
                f.write(b + "\n" + "-" * 60 + "\n")
        print(f"Successfully generated: {out2}")

if __name__ == "__main__":
    merge()
    print("\nMerging process finished.")
