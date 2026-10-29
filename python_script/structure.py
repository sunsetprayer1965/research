import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
from scipy import interpolate
import os

def interpolate_chromosomes(X, lengths, eps=1e-1, smooth=0):
    print(f"Debug: Starting interpolation with lengths: {lengths}")
    print(f"Debug: Input X shape: {X.shape}")
    smoothed_X = [[] for _ in range(len(lengths))]  # Reconstruct all chromosomes
    mask = np.invert(np.isnan(X[:, 0]))
    begin, end = 0, 0
    
    for i, length in enumerate(lengths):
        end += length
        print(f"Debug: Processing chromosome {i + 1}")
        print(f"Debug: Current indices: begin={begin}, end={end}")
        
        if end > X.shape[0]:
            print(f"Error: 'end' index ({end}) exceeds the size of X ({X.shape[0]}).")
            return None
        
        x, y, z = X[begin:end, 0], X[begin:end, 1], X[begin:end, 2]
        
        if len(x) == 0:
            break
        
        indx = np.invert(np.isnan(x))
        m = np.arange(len(x))[indx]
        
        f_x = interpolate.Rbf(m, x[indx], smooth=smooth)
        f_y = interpolate.Rbf(m, y[indx], smooth=smooth)
        f_z = interpolate.Rbf(m, z[indx], smooth=smooth)
        
        for j in range(length - 1):
            m = np.arange(j, j + 1, 0.1)
            try:
                smoothed_X[i].append([
                    f_x(np.arange(m.min(), m.max(), 0.1)),
                    f_y(np.arange(m.min(), m.max(), 0.1)),
                    f_z(np.arange(m.min(), m.max(), 0.1))
                ])
            except Exception as e:
                print(f"Error while interpolating chromosome {i + 1} at segment {j}: {e}")
                return None
        
        begin = end
    
    return smoothed_X

# 16条染色体长度
theLengths = np.array([24, 82, 32, 154, 58, 28, 110, 57, 44, 75, 67, 108, 93, 79, 110, 95])

# === 只新增：循环处理四个文件 ===
files = ["PM1.structure", "PM2.structure", "mds.structure", "nmds.structure"]

for fname in files:
    if not os.path.exists(fname):
        print(f"[WARN] {fname} 不存在（注意大小写）——跳过")
        continue

    # Load coordinate data（默认按任意空白分隔更稳）
    coordinates = np.loadtxt(fname)

    # 逐条切分出16条染色体坐标
    chromosome_data = []
    start, end = 0, 0
    for length in theLengths:
        end += length
        chromosome_data.append((
            coordinates[start:end, 0],
            coordinates[start:end, 1],
            coordinates[start:end, 2]
        ))
        start = end

    # RBF 插值
    chroms = interpolate_chromosomes(coordinates, theLengths)
    if chroms is None:
        print(f"[ERROR] 插值失败：{fname}，已跳过绘图。")
        continue

    # 绘图
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # 配色与标签
    colors = [
        'r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange', 
        'purple', 'pink', 'lime', 'teal', 'brown', 'navy', 'gold', 'grey'
    ]
    labels = [f'Chromosome {i + 1}' for i in range(16)]

    for i, (x, y, z) in enumerate(chromosome_data):
        ax.scatter(x, y, z, color=colors[i], label=labels[i], s=30)
        for segment in chroms[i]:
            ax.plot(segment[0], segment[1], segment[2], color=colors[i])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()

    stem = fname.rsplit(".", 1)[0]
    title = f'inferred_3d_structures_{stem}_HindIII+fdr1%_william'
    plt.title(title)

    # 先保存，再show（防止无显示环境导致保存失败）
    out_png = f'{title}.png'
    plt.savefig(out_png, dpi=300)
    print(f"[OK] 保存图片 -> {out_png}")

    # 若在本地有GUI，依然可以显示；没有也不会影响保存
    plt.show()
    plt.close()
