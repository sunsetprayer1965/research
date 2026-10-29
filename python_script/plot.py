import numpy as np
import matplotlib.pyplot as plt

# ===== Step 1. Load data =====
file_A = "output_chroms_msplPM1_linking_number.txt"   # x-values (e.g. 1vs2, 1vs3, ...)
file_B = "output_chroms_msplPM2_linking_number.txt"   # y-values (e.g. 2vs1, 3vs1, ...)

x = np.loadtxt(file_A)
y = np.loadtxt(file_B)

if len(x) != 120 or len(y) != 120:
    raise ValueError(f"Expected 120 numbers in each file, got {len(x)} and {len(y)}")

# ===== Step 2. Plot scatter =====
plt.figure(figsize=(6,6))
plt.scatter(x, y, c='steelblue', edgecolors='k')

plt.xlabel("PM1_linking_number")
plt.ylabel("PM2_linking_number")
plt.grid(True, linestyle="--", alpha=0.7)
plt.xlim(0, 120)
plt.ylim(0, 120)




plt.tight_layout()
plt.savefig("PM1 vs PM2.png", dpi=300)
