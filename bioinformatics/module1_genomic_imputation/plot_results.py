"""
plot_results.py
---------------
Produces reproducible summary figures from results/benchmark_results.csv.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_DIR / "results"
RESULTS_FILE = RESULTS_DIR / "benchmark_results.csv"

df = pd.read_csv(RESULTS_FILE)

# ---- Figure 1: heatmap of dosage r2 ----
pivot = df.pivot(
    index="panel_size",
    columns="depth",
    values="mean_dosage_r2",
)

fig, ax = plt.subplots(figsize=(7, 5))
im = ax.imshow(
    pivot.values,
    aspect="auto",
    cmap="viridis",
    vmin=0.5,
    vmax=1.0,
)

ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels([f"{d}x" for d in pivot.columns])
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index)

ax.set_xlabel("Target sequencing depth")
ax.set_ylabel("Reference panel size (haplotypes)")
ax.set_title(
    "Imputation accuracy (dosage $r^2$)\n"
    "Li-Stephens HMM, low-recombination simulated panel"
)

for i in range(pivot.shape[0]):
    for j in range(pivot.shape[1]):
        value = pivot.values[i, j]
        ax.text(
            j,
            i,
            f"{value:.2f}",
            ha="center",
            va="center",
            color="white" if value < 0.85 else "black",
            fontsize=8,
        )

fig.colorbar(im, ax=ax, label="Dosage $r^2$")
fig.tight_layout()

heatmap_file = RESULTS_DIR / "heatmap_dosage_r2.png"
fig.savefig(heatmap_file, dpi=150)
plt.close(fig)

# ---- Figure 2: concordance vs depth ----
fig, ax = plt.subplots(figsize=(7, 5))

for panel_size, sub in df.groupby("panel_size"):
    sub = sub.sort_values("depth")
    ax.plot(
        sub["depth"],
        sub["mean_concordance"],
        marker="o",
        label=f"n={panel_size} ref haps",
    )

ax.set_xscale("log")
ax.set_xlabel("Target sequencing depth (x, log scale)")
ax.set_ylabel("Genotype concordance")
ax.set_title(
    "Imputation concordance vs. depth\n"
    "by reference panel size"
)
ax.legend(title="Panel size")
ax.grid(alpha=0.3)
fig.tight_layout()

concordance_file = RESULTS_DIR / "concordance_vs_depth.png"
fig.savefig(concordance_file, dpi=150)
plt.close(fig)

print(f"Saved {heatmap_file.relative_to(PROJECT_DIR)}")
print(f"Saved {concordance_file.relative_to(PROJECT_DIR)}")
