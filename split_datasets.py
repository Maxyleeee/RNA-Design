import os

def split_dataset(filepath):
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return

    basename = os.path.splitext(filepath)[0]
    file_count1 = f"{basename}_count1.txt"
    file_count2plus = f"{basename}_count2plus.txt"

    c1_list = []
    c2_list = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if "," in line:
                struct, count_str = line.split(",", 1)
                try:
                    count = int(count_str.strip())
                except:
                    count = 0
                
                if count == 1:
                    c1_list.append(line)
                elif count >= 2:
                    c2_list.append(line)
            else:
                # Fallback for lines without count (shouldn't happen with new generation)
                pass

    if c1_list:
        with open(file_count1, 'w') as f:
            for item in c1_list:
                f.write(item + "\n")
        print(f"Saved {len(c1_list)} structures to {file_count1}")
    
    if c2_list:
        with open(file_count2plus, 'w') as f:
            for item in c2_list:
                f.write(item + "\n")
        print(f"Saved {len(c2_list)} structures to {file_count2plus}")

if __name__ == "__main__":
    datasets = [
        "/home/maxyle/RNA/output/dataset_l50_h4_t3_wu1_motifs30.txt",
        "/home/maxyle/RNA/output/dataset_l50_h4_t3_wu1_motifs40.txt",
        "/home/maxyle/RNA/output/dataset_l50_h4_t3_wu1_motifs50.txt"
    ]
    for d in datasets:
        split_dataset(d)
