import numpy as np
from scipy.stats import chi2
import itertools
import os

files = [
    "MDS_split_linking_frequency.txt",
    "NMDS_split_linking_frequency.txt",
    "PM1_split_linking_frequency.txt",
    "PM2_split_linking_frequency.txt"
]

folder = "/data"
data_dict = {}

for file in files:
    path = os.path.join(folder, file)   # /data/MDS_linking_number.txt
    name = os.path.splitext(file)[0]
    data = np.loadtxt(path)
    data = data * 10
    data_dict[name] = data

N = 1000
df = 120

print("\n====== Pairwise Pearson Chi-square Test ======\n")

for (name1, data1), (name2, data2) in itertools.combinations(data_dict.items(), 2):
    numerator = (data1 - data2) ** 2
    denominator = (data1 + data2) * (2 * N - (data1 + data2))
    mask = denominator != 0
    X2 = 2 * N * np.sum(numerator[mask] / denominator[mask])
    p_value = 1 - chi2.cdf(X2, df)

    print(f"Comparing: {name1} vs {name2}")
    print(f"  Pearson Chi-square X² = {X2:.4f}")
    print(f"  p-value = {p_value:.6f}\n")

# docker run -it --rm `
# -v "C:\Users\weilu\Desktop\research william luo\python_script:/app" `
# -v "C:\Users\weilu\Desktop\research william luo\dataset_linking_number\HindIII+Msel(no FDR)\linking:/data" `
# pastis python /app/chi-square.py