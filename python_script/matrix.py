import csv
import matplotlib.pyplot as plt
import numpy as np

N = 16
roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII",
         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]

def read_numbers(path):
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read().replace(',', ' ')
    nums = [float(x) for x in text.split()]
    if len(nums) != N*(N-1)//2:
        raise ValueError(f"{path} 应该有 {N*(N-1)//2} 个数字，实际 {len(nums)}")
    return nums

def order_upper():
    for i in range(N-1):
        for j in range(i+1, N):
            yield i, j

def order_lower():
    for j in range(N-1):
        for i in range(j+1, N):
            yield i, j

# === 1. 读入两个 txt ===
a = read_numbers("Panel_1_linking_frequency.txt")  # 上三角
b = read_numbers("panel1-thesis.txt")              # 下三角

# === 2. 填充矩阵 ===
M = [[None for _ in range(N)] for _ in range(N)]
k = 0
for i, j in order_upper():
    M[i][j] = a[k]; k += 1
k = 0
for i, j in order_lower():
    M[i][j] = b[k]; k += 1

# === 3. 保存为 CSV，并加罗马数字 ===
with open("matrix_with_roman.csv", "w", newline='', encoding="utf-8") as f:
    writer = csv.writer(f)
    # 表头
    writer.writerow([""] + roman)
    # 每行加行标签
    for i, row in enumerate(M):
        writer.writerow([roman[i]] + ["" if v is None else f"{v:.1f}" for v in row])

# === 4. 作图 ===
grid = np.array([[float('nan') if v is None else v for v in row] for row in M])
plt.imshow(grid, interpolation='nearest', cmap='coolwarm')
plt.colorbar(label="Linking Frequency")
plt.xticks(ticks=range(N), labels=roman, rotation=45)
plt.yticks(ticks=range(N), labels=roman)
plt.title("Matrix from A.txt and B.txt")
plt.tight_layout()
plt.savefig("matrix_with_roman.png", dpi=300)
