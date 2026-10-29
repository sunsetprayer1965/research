import pandas as pd
import numpy as np
df = pd.read_csv("restriction_fragments_mappable_HindIII.txt",sep=",",skiprows=2,names=['chr', 'left', 'mid', 'right','mappable'])
print(df)
max_right_per_chr = df.groupby("chr")["right"].max()
# Create a bed file with 10,000 resolution for chr 1 to 5
bed_rows = []
id_counter = -1

for chr_num, max_right in max_right_per_chr.items():
    if chr_num >8:
        continue
    for start in range(1, max_right + 1, 10000):
        end = start + 9999
        id_counter += 1
        bed_rows.append([f"chr{chr_num:02d}", start, end, id_counter])

# Create DataFrame for the BED file
bed_df = pd.DataFrame(bed_rows, columns=["chr", "start", "end", "id"])

# Save to a file
bed_file_name = "duan.SC.10000.raw_sub.bed"
bed_df.to_csv(bed_file_name, sep="\t", index=False, header=False)

print(f"BED file '{bed_file_name}' has been created.")
