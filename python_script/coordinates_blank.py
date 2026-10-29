import pandas as pd
import numpy as np

import pandas as pd

# ADDED: 需要处理的4个文件（注意大小写要与磁盘一致）
files = ['MDS', 'NMDS', 'PM1', 'PM2']

# Define the chromosome grouping criteria
chromosome_counts = [24,82,32,154,58,28,110,57,44,75,67,108,93,79,110,95]  # Number of rows for each chromosome

for base in files:
    output_data = []  # 每个文件单独的输出缓冲
    # Load your data (replace 'your_file.csv' with your actual file path)
    # Assuming your file is tab-delimited and has no header
    df = pd.read_csv(f'{base}.structure', header=None, delim_whitespace=True)

    # Initialize the final DataFrame to store the grouped data
    final_data = pd.DataFrame()

    # Start row index
    start_index = 0

    # Loop through each chromosome count and slice rows
    for count in chromosome_counts:
        # Extract rows for the current chromosome
        chromosome_data = df.iloc[start_index:start_index + count]
        chromosome_data = chromosome_data.dropna()
        # Append the chromosome rows to the output data
        output_data.extend(chromosome_data.values.tolist())
        
        # Append a blank row (as an empty list)
        output_data.append([])
        
        # Update the start index for the next chromosome
        start_index += count

    # Save the final output
    # Write the final data to a CSV file
    out_name = f'{base}.coordinate.txt'  # ADDED: 输出文件按输入前缀命名
    with open(out_name, 'w', newline='') as f:
        for row in output_data:
            if row:  # Write regular rows
                f.write(','.join(map(str, row)) + '\n')
            else:  # Write a truly blank row
                f.write('\n')

    print(f"{base}: Data grouped by chromosomes with proper blank rows saved successfully -> {out_name}!")

# nan_row_indices = df[df.isnull().any(axis=1)].index

# # Print the indices
# print("Indices of rows containing NaN values:", nan_row_indices.tolist())
