import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
matplotlib.use('Agg')

from iced.datasets import load_sample_yeast
from iced import filter, normalization
from pastis.datasets import generate_dataset_from_distances
from sklearn.metrics.pairwise import euclidean_distances
from pastis.optimization.mds import estimate_X

# counts, lengths = load_sample_yeast()
# counts = filter.filter_low_counts(counts, remove_all_zeros_loci=True, sparsity=False)
# print(counts[0])
lengths = np.array([31])
df = pd.read_csv('chromosome-first_5.matrix', sep='\t', header=None)

# Filter rows where both the first and second column values are between 106 and 137 (inclusive).
df = df[(df[0] >= 106) & (df[0] <= 137) & (df[1] >= 106) & (df[1] <= 137)]
matrix_df = df.pivot(index=0, columns=1, values=2)

# Define full list of loci (106 to 137 inclusive)
loci = list(range(106, 138))

matrix_df = matrix_df.reindex(index=loci, columns=loci)

# Create a copy to build the symmetric matrix.
symmetric_matrix_df = matrix_df.copy()

# Fill in the lower triangle using the upper triangle values.
# We assume that our data is provided in the upper triangle.
for i in range(len(loci)):
    for j in range(len(loci)):
        if i < j:
            # Get the upper triangle value
            value = symmetric_matrix_df.iloc[i, j]
            # Copy it to the corresponding lower triangle position.
            symmetric_matrix_df.iloc[j, i] = value

# Convert the symmetric DataFrame to a list of lists.
counts = symmetric_matrix_df.values.tolist()
counts = np.array(counts)
counts = normalization.ICE_normalization(counts)

# structure = estimate_X(counts)
# distances = euclidean_distances(structure, structure)

# ##############################################################################
# # # Now, generate the contact count matrix from the distance matrix
# # negative_binomial_count = generate_dataset_from_distances(
# #     distances, alpha=-3, distribution="NegativeBinomial")

# poisson_count = generate_dataset_from_distances(
#     distances, alpha=-3, distribution="Poisson", lengths=lengths,
#     alpha_inter=-2)

# ##############################################################################
# Plot the resulting contact count matrices
# fig, axes = plt.subplots(figsize=(7, 7),
#                          ncols=2, nrows=2, tight_layout=True)

def plot_matrix(ax, matrix, cmap="RdBu_r", norm=None, lengths=None):
    if lengths is None:
        lengths = np.array([matrix.shape[0]])
    ax.matshow(matrix, cmap=cmap, norm=norm,
               extent=(0, lengths.sum(), 0, lengths.sum()))
    [ax.axhline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
    [ax.axvline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
    ax.set_xticks([])
    ax.set_yticks([])

# # # First, plot the original contact count matrix
# plot_matrix(counts, norm=colors.SymLogNorm(1), lengths=lengths)
# axes[0, 0].set_title("Contact counts", fontweight="bold")

# # Plot the distance matrix used as input to the generated data function
# plot_matrix(axes[0, 1], distances, cmap="viridis", lengths=lengths)
# axes[0, 1].set_title("Distance matrix", fontweight="bold")

# # The negative binomial counts
# plot_matrix(axes[1, 0], negative_binomial_count, norm=colors.SymLogNorm(1), lengths=lengths)
# axes[1, 0].set_title("Simulated dataset (NB)", fontweight="bold")

# The Poisson dataset
# plot_matrix(axes[1, 1], poisson_count, norm=colors.SymLogNorm(1), lengths=lengths)
# axes[1, 1].set_title("Simulated dataset (Poisson)", fontweight="bold")

# Display the plots

fig, ax = plt.subplots(figsize=(7, 7), tight_layout=True)

def plot_matrix(ax, matrix, cmap="Blues", norm=None, lengths=None):
    if lengths is None:
        lengths = np.array([matrix.shape[0]])
    ax.matshow(matrix, cmap=cmap, norm=norm,
               extent=(0, lengths.sum(), 0, lengths.sum()))
    for i in lengths.cumsum():
        ax.axhline(i, linewidth=1, color="#000000")
        ax.axvline(i, linewidth=1, color="#000000")
    ax.set_xticks([])
    ax.set_yticks([])

# Plot on the single axes object 'ax'
plot_matrix(ax, counts, norm=colors.SymLogNorm(1), lengths=lengths)
ax.set_title("Contact counts", fontweight="bold")
plt.show()
plt.savefig('/app/matrix_plot.png')   # Use 'r' to avoid escape character issues or use double backslashes
plt.show()
