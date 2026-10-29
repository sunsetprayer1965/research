import numpy as np
import pandas as pd
import os

# ========= 参数 =========
folder = "/data"
files = [f"{i}.txt" for i in range(1, 13)]
labels = [str(i) for i in range(1, 13)]

N = 1000

# ========= 更稳的读取函数：按行读，每行取一个数字 =========
def read_single_column_txt(path):
    values = []
    with open(path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()

            # 跳过空行
            if not line:
                continue

            parts = line.split()

            # 如果一行里不止一个数字，直接报错，方便定位
            if len(parts) != 1:
                raise ValueError(
                    f"{os.path.basename(path)} 第 {line_num} 行不是单个数字，而是: {parts}"
                )

            values.append(float(parts[0]))

    return np.array(values, dtype=float)

# ========= 读取全部数据 =========
data_dict = {}
length_dict = {}

for file in files:
    path = os.path.join(folder, file)
    name = os.path.splitext(file)[0]

    data = read_single_column_txt(path)
    data = data * 10

    data_dict[name] = data
    length_dict[name] = len(data)

# ========= 检查长度 =========
print("Each file length:")
for name in labels:
    print(f"{name}.txt -> {length_dict[name]}")

unique_lengths = set(length_dict.values())
if len(unique_lengths) != 1:
    raise ValueError("不是所有 txt 的长度都一致，请检查文件。")

# 如果你就是要求必须 120 行，也可以强检查
expected_len = 120
for name in labels:
    if length_dict[name] != expected_len:
        raise ValueError(f"{name}.txt 实际读取到 {length_dict[name]} 行，不是 120 行。")

# ========= 创建结果矩阵（先全空） =========
n = len(labels)
result_matrix = np.full((n, n), "", dtype=object)

# ========= 只计算下三角 =========
for i in range(n):
    for j in range(i):   # 只保留左下角：i > j
        data1 = data_dict[labels[i]]
        data2 = data_dict[labels[j]]

        numerator = (data1 - data2) ** 2
        denominator = (data1 + data2) * (2 * N - (data1 + data2))

        valid = denominator != 0
        X2 = 2 * N * np.sum(numerator[valid] / denominator[valid])

        result_matrix[i, j] = round(X2, 4)

# ========= 转成 DataFrame =========
df_result = pd.DataFrame(result_matrix, index=labels, columns=labels)

# ========= 导出 Excel =========
output_path = "/app/chi_square_lower_triangle.xlsx"
df_result.to_excel(output_path)

print(f"\nExcel saved to: {output_path}")