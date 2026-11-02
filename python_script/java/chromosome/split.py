import os

def split_chromosomes(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    with open(input_file, 'r') as f:
        lines = f.readlines()

    chrom_data = []
    chrom_index = 1

    for line in lines:
        if line.strip() == "":
            if chrom_data:
                output_path = os.path.join(output_dir, f"{chrom_index}.txt")
                with open(output_path, 'w') as out_file:
                    out_file.writelines(chrom_data)
                print(f"✅ Written: {output_path}")
                chrom_index += 1
                chrom_data = []
        else:
            chrom_data.append(line)

    # Write the final chromosome, if any.
    if chrom_data:
        output_path = os.path.join(output_dir, f"{chrom_index}.txt")
        with open(output_path, 'w') as out_file:
            out_file.writelines(chrom_data)
        print(f"✅ Written: {output_path}")

    print(f"🎉 Done! Total chromosomes split: {chrom_index - 1}\n")


# ==========================
# Process multiple files in one run.
# ==========================

base_path = "."

files = [
    "MDS.coordinate.txt",
    "NMDS.coordinate.txt",
    "PM1.coordinate.txt",
    "PM2.coordinate.txt"
]

for file in files:
    input_path = os.path.join(base_path, file)
    output_folder = file.replace(".coordinate.txt", "_split")
    split_chromosomes(input_path, output_folder)
