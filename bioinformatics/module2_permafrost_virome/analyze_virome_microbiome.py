"""
analyze_virome_microbiome.py
------------------------------
Analysis pipeline for the simulated permafrost chronosequence data
(microbiome_abundance.csv, virome_abundance.csv, sample_metadata.csv from
simulate_metagenome.py). Implements the three core analyses the project
brief calls for:

  1. Microbial community diversity/composition changes across the
     chronosequence (alpha diversity trend + beta-diversity ordination).
  2. Same for the virome.
  3. In silico virus-host linkage: predicting which bacterial taxa each
     phage infects, using co-abundance correlation across samples (a
     standard, widely-used host-prediction signal in real viral
     metagenomics, complementary to sequence-composition/CRISPR-spacer
     methods used e.g. by iPHoP/vConTACT2 on real data).

Outputs (written to results/):
  - alpha_diversity.csv, alpha_diversity_trend.png
  - beta_diversity_ordination.png (bacteria + phage, PCoA via classical MDS
    on Bray-Curtis dissimilarity)
  - predicted_host_links.csv (+ precision/recall against the simulator's
    ground truth, since this is a validated demo)
  - virus_host_correlation_heatmap.png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr


def shannon(row):
    p = row[row > 0]
    return -np.sum(p * np.log(p))


def simpson(row):
    p = row[row > 0]
    return 1 - np.sum(p ** 2)


def alpha_diversity_table(abund_df, label):
    df = pd.DataFrame({
        f"{label}_shannon": abund_df.apply(shannon, axis=1),
        f"{label}_simpson": abund_df.apply(simpson, axis=1),
        f"{label}_richness": (abund_df > 0).sum(axis=1),
    })
    return df


def bray_curtis(abund_df):
    dist = pdist(abund_df.values, metric="braycurtis")
    return squareform(dist)


def ordinate(abund_df, seed=0):
    dmat = bray_curtis(abund_df)
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=seed,
              normalized_stress="auto")
    coords = mds.fit_transform(dmat)
    return coords


def predict_virus_host_links(bact_df, phage_df, top_n=1, corr_threshold=0.6):
    """For each phage, correlate its abundance profile (across samples)
    against every bacterial taxon's profile and call the strongest positive
    correlate(s) above threshold as predicted host(s). This mirrors
    co-abundance-based host prediction used as one line of evidence in real
    virus-host linkage workflows."""
    records = []
    corr_matrix = pd.DataFrame(index=phage_df.columns, columns=bact_df.columns, dtype=float)
    for phage in phage_df.columns:
        corrs = {}
        for host in bact_df.columns:
            rho, _ = spearmanr(phage_df[phage], bact_df[host])
            corrs[host] = rho
            corr_matrix.loc[phage, host] = rho
        ranked = sorted(corrs.items(), key=lambda kv: kv[1], reverse=True)
        for host, rho in ranked[:top_n]:
            if rho >= corr_threshold:
                records.append({"phage": phage, "predicted_host": host, "spearman_rho": rho})
    return pd.DataFrame.from_records(records), corr_matrix


def main():
    bact_df = pd.read_csv("results/microbiome_abundance.csv", index_col=0)
    phage_df = pd.read_csv("results/virome_abundance.csv", index_col=0)
    metadata = pd.read_csv("results/sample_metadata.csv", index_col=0)
    true_links = pd.read_csv("results/true_host_links_for_validation.csv")

    # --- 1. Alpha diversity across the chronosequence ---
    alpha_bact = alpha_diversity_table(bact_df, "microbiome")
    alpha_phage = alpha_diversity_table(phage_df, "virome")
    alpha = alpha_bact.join(alpha_phage).join(metadata)
    alpha.to_csv("results/alpha_diversity.csv")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(alpha["relative_age"], alpha["microbiome_shannon"], c="steelblue")
    axes[0].set_xlabel("Relative sample age (0=thaw front, 1=deep Pleistocene)")
    axes[0].set_ylabel("Shannon diversity")
    axes[0].set_title("Microbiome diversity across chronosequence")
    z = np.polyfit(alpha["relative_age"], alpha["microbiome_shannon"], 1)
    xs = np.linspace(0, 1, 50)
    axes[0].plot(xs, np.polyval(z, xs), "k--", alpha=0.6)

    axes[1].scatter(alpha["relative_age"], alpha["virome_shannon"], c="indianred")
    axes[1].set_xlabel("Relative sample age (0=thaw front, 1=deep Pleistocene)")
    axes[1].set_ylabel("Shannon diversity")
    axes[1].set_title("Virome diversity across chronosequence")
    z2 = np.polyfit(alpha["relative_age"], alpha["virome_shannon"], 1)
    axes[1].plot(xs, np.polyval(z2, xs), "k--", alpha=0.6)
    fig.tight_layout()
    fig.savefig("results/alpha_diversity_trend.png", dpi=150)
    plt.close(fig)

    # --- 2. Beta diversity ordination ---
    bact_coords = ordinate(bact_df)
    phage_coords = ordinate(phage_df)
    age = metadata["relative_age"].values

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sc0 = axes[0].scatter(bact_coords[:, 0], bact_coords[:, 1], c=age, cmap="viridis", s=70)
    axes[0].set_title("Microbiome PCoA (Bray-Curtis)")
    axes[0].set_xlabel("MDS1"); axes[0].set_ylabel("MDS2")
    fig.colorbar(sc0, ax=axes[0], label="Relative age")

    sc1 = axes[1].scatter(phage_coords[:, 0], phage_coords[:, 1], c=age, cmap="viridis", s=70)
    axes[1].set_title("Virome PCoA (Bray-Curtis)")
    axes[1].set_xlabel("MDS1"); axes[1].set_ylabel("MDS2")
    fig.colorbar(sc1, ax=axes[1], label="Relative age")
    fig.tight_layout()
    fig.savefig("results/beta_diversity_ordination.png", dpi=150)
    plt.close(fig)

    # --- 3. Virus-host linkage prediction ---
    predicted, corr_matrix = predict_virus_host_links(bact_df, phage_df,
                                                        top_n=1, corr_threshold=0.5)
    predicted.to_csv("results/predicted_host_links.csv", index=False)

    merged = predicted.merge(true_links, on="phage", how="left")
    merged["correct"] = merged["predicted_host"] == merged["true_host"]
    n_predicted = len(merged)
    n_correct = merged["correct"].sum()
    n_total_phage = phage_df.shape[1]
    precision = n_correct / n_predicted if n_predicted else float("nan")
    recall = n_correct / n_total_phage
    print(f"Virus-host linkage: predicted links for {n_predicted}/{n_total_phage} phage "
          f"(threshold rho>=0.5)")
    print(f"  Precision (of predictions made): {precision:.2f}")
    print(f"  Recall (of all phage, incl. no-call): {recall:.2f}")
    merged.to_csv("results/predicted_host_links_validated.csv", index=False)

    # heatmap of correlation matrix (subset for readability if large)
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(corr_matrix.values.astype(float), aspect="auto", cmap="RdBu_r",
                    vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr_matrix.columns)))
    ax.set_xticklabels(corr_matrix.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(len(corr_matrix.index)))
    ax.set_yticklabels(corr_matrix.index, fontsize=6)
    ax.set_title("Phage-host co-abundance correlation (Spearman ρ)\nacross permafrost chronosequence samples")
    fig.colorbar(im, ax=ax, label="Spearman ρ")
    fig.tight_layout()
    fig.savefig("results/virus_host_correlation_heatmap.png", dpi=150)
    plt.close(fig)

    print("\nSaved: alpha_diversity.csv, alpha_diversity_trend.png, "
          "beta_diversity_ordination.png, predicted_host_links.csv, "
          "predicted_host_links_validated.csv, virus_host_correlation_heatmap.png")


if __name__ == "__main__":
    main()
