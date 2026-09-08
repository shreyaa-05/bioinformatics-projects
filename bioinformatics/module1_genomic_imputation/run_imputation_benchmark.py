"""
run_imputation_benchmark.py
----------------------------
Benchmarks a Li-Stephens-style haplotype-copying HMM across sequencing depths
and reference-panel sizes using the reproducible simulated dataset in data/.
"""

from pathlib import Path
import numpy as np
import pandas as pd

from simulate_data import (
    REFERENCE_FILE,
    TARGET_FILE,
    generate_and_save_dataset,
    simulate_low_coverage,
)

PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_DIR / "results"
RESULTS_FILE = RESULTS_DIR / "benchmark_results.csv"


def li_stephens_forward_backward(
    observed,
    ref_panel,
    recomb_rate=0.003,
    theta=0.001,
):
    """Run a simplified Li-Stephens forward-backward HMM for one haplotype."""
    n_ref, n_snps = ref_panel.shape

    stay = 1.0 - recomb_rate
    jump_each = recomb_rate / n_ref

    def emission(site_idx, allele_call):
        if allele_call == -1:
            return np.ones(n_ref)
        ref_alleles = ref_panel[:, site_idx]
        match = ref_alleles == allele_call
        return np.where(match, 1.0 - theta, theta)

    alpha = np.empty((n_snps, n_ref))
    alpha[0] = emission(0, observed[0]) / n_ref
    alpha[0] /= alpha[0].sum()

    for j in range(1, n_snps):
        prev = alpha[j - 1]
        predicted = stay * prev + jump_each * prev.sum()
        alpha[j] = predicted * emission(j, observed[j])
        alpha[j] /= alpha[j].sum()

    beta = np.empty((n_snps, n_ref))
    beta[-1] = 1.0

    for j in range(n_snps - 2, -1, -1):
        emission_next = emission(j + 1, observed[j + 1])
        tmp = beta[j + 1] * emission_next
        beta[j] = stay * tmp + jump_each * tmp.sum()
        beta[j] /= beta[j].sum()

    post = alpha * beta
    post /= post.sum(axis=1, keepdims=True)

    return np.einsum("jk,kj->j", post, ref_panel)


def impute_target(observed, ref_panel, recomb_rate, theta=0.001):
    dosage = li_stephens_forward_backward(
        observed, ref_panel, recomb_rate, theta
    )
    hard_call = (dosage >= 0.5).astype(np.int8)
    return dosage, hard_call


def evaluate(true_hap, hard_call, dosage, observed_mask):
    """Evaluate only sites that were originally missing."""
    missing = ~observed_mask

    if missing.sum() == 0:
        return np.nan, np.nan

    concordance = (
        hard_call[missing] == true_hap[missing]
    ).mean()

    true_values = true_hap[missing].astype(float)
    dosage_values = dosage[missing]

    if true_values.std() == 0 or dosage_values.std() == 0:
        r2 = np.nan
    else:
        r2 = np.corrcoef(true_values, dosage_values)[0, 1] ** 2

    return concordance, r2


def load_or_create_dataset():
    """Load saved data, or create it if the data files do not exist."""
    if not REFERENCE_FILE.exists() or not TARGET_FILE.exists():
        print("Saved dataset not found. Generating reproducible dataset...")
        generate_and_save_dataset()

    ref_panel = np.load(REFERENCE_FILE)
    true_targets = np.load(TARGET_FILE)

    print(
        f"Loaded reference panel: {ref_panel.shape[0]} haplotypes x "
        f"{ref_panel.shape[1]} SNPs"
    )
    print(f"Loaded target haplotypes: {true_targets.shape[0]}")

    return ref_panel, true_targets


def run_benchmark(
    depths=(0.1, 0.25, 0.5, 1.0, 2.0, 4.0),
    panel_sizes=(20, 50, 100, 200),
    recomb_rate=0.003,
    seed=2026,
):
    ref_panel_full, true_targets = load_or_create_dataset()
    rng = np.random.default_rng(seed)

    if max(panel_sizes) > ref_panel_full.shape[0]:
        raise ValueError(
            "Requested panel size exceeds the number of saved reference "
            "haplotypes."
        )

    records = []

    for panel_size in panel_sizes:
        ref_panel = ref_panel_full[:panel_size]

        for depth in depths:
            concordances = []
            r2s = []

            for true_hap in true_targets:
                observed, mask = simulate_low_coverage(
                    true_hap[None, :],
                    depth=depth,
                    rng=rng,
                )
                dosage, hard_call = impute_target(
                    observed[0],
                    ref_panel,
                    recomb_rate,
                )
                conc, r2 = evaluate(
                    true_hap,
                    hard_call,
                    dosage,
                    mask[0],
                )
                concordances.append(conc)
                r2s.append(r2)

            record = {
                "panel_size": panel_size,
                "depth": depth,
                "mean_concordance": np.nanmean(concordances),
                "mean_dosage_r2": np.nanmean(r2s),
                "n_targets": len(true_targets),
            }
            records.append(record)

            print(
                f"panel_size={panel_size:>4}  "
                f"depth={depth:>4}x  "
                f"concordance={record['mean_concordance']:.3f}  "
                f"dosage_r2={record['mean_dosage_r2']:.3f}"
            )

    return pd.DataFrame.from_records(records)


if __name__ == "__main__":
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    df = run_benchmark()
    df.to_csv(RESULTS_FILE, index=False)

    print(f"\nSaved {RESULTS_FILE.relative_to(PROJECT_DIR)}")
