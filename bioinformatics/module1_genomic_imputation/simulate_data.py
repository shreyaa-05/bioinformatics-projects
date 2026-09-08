"""
simulate_data.py
-----------------
Creates a reproducible simulated haplotype dataset for benchmarking genotype
imputation in a non-model system.

The simulation uses a founder-mosaic model:
1. Founder haplotypes are generated across biallelic SNPs.
2. Reference and held-out target haplotypes are recombinant mosaics of founders.
3. The generated ground truth is saved to disk so downstream benchmarking can
   be reproduced without relying on hidden in-memory data.
"""

from pathlib import Path
import json
import numpy as np

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"

REFERENCE_FILE = DATA_DIR / "reference_panel.npy"
TARGET_FILE = DATA_DIR / "target_haplotypes.npy"
METADATA_FILE = DATA_DIR / "dataset_metadata.json"


def simulate_founders(
    n_founders: int,
    n_snps: int,
    rng: np.random.Generator,
    maf_alpha: float = 0.4,
    maf_beta: float = 0.4,
) -> np.ndarray:
    """Simulate founder haplotypes at biallelic SNPs."""
    freqs = rng.beta(maf_alpha, maf_beta, size=n_snps)
    return (rng.random((n_founders, n_snps)) < freqs).astype(np.int8)


def mosaic_haplotypes(
    founders: np.ndarray,
    n_out: int,
    recomb_rate: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate haplotypes as recombinant mosaics of founder haplotypes."""
    n_founders, n_snps = founders.shape
    out = np.empty((n_out, n_snps), dtype=np.int8)

    for i in range(n_out):
        current_founder = rng.integers(0, n_founders)
        hap = np.empty(n_snps, dtype=np.int8)
        hap[0] = founders[current_founder, 0]

        switches = rng.random(n_snps - 1) < recomb_rate
        for j in range(1, n_snps):
            if switches[j - 1]:
                current_founder = rng.integers(0, n_founders)
            hap[j] = founders[current_founder, j]

        out[i] = hap

    return out


def simulate_low_coverage(
    target_haps: np.ndarray,
    depth: float,
    rng: np.random.Generator,
):
    """Mask sites according to a Poisson-coverage model.

    P(missing) = exp(-depth). Missing sites are encoded as -1.
    """
    p_missing = np.exp(-depth)
    mask = rng.random(target_haps.shape) >= p_missing
    observed = np.where(mask, target_haps, -1).astype(np.int8)
    return observed, mask


def build_dataset(
    n_founders: int = 25,
    n_ref_haps: int = 200,
    n_target_haps: int = 15,
    n_snps: int = 800,
    recomb_rate: float = 0.003,
    seed: int = 1,
):
    """Build one reference panel and held-out ground-truth target set."""
    rng = np.random.default_rng(seed)
    founders = simulate_founders(n_founders, n_snps, rng)
    ref_panel = mosaic_haplotypes(founders, n_ref_haps, recomb_rate, rng)
    true_targets = mosaic_haplotypes(
        founders, n_target_haps, recomb_rate, rng
    )
    return ref_panel, true_targets


def save_dataset(ref_panel: np.ndarray, true_targets: np.ndarray, metadata: dict):
    """Save simulated ground-truth data and simulation metadata."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    np.save(REFERENCE_FILE, ref_panel)
    np.save(TARGET_FILE, true_targets)

    with open(METADATA_FILE, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print(f"Saved {REFERENCE_FILE.relative_to(PROJECT_DIR)}")
    print(f"Saved {TARGET_FILE.relative_to(PROJECT_DIR)}")
    print(f"Saved {METADATA_FILE.relative_to(PROJECT_DIR)}")


def generate_and_save_dataset(
    n_founders: int = 25,
    n_ref_haps: int = 200,
    n_target_haps: int = 15,
    n_snps: int = 800,
    recomb_rate: float = 0.003,
    seed: int = 1,
):
    """Generate and save the exact dataset used by the benchmark."""
    ref_panel, true_targets = build_dataset(
        n_founders=n_founders,
        n_ref_haps=n_ref_haps,
        n_target_haps=n_target_haps,
        n_snps=n_snps,
        recomb_rate=recomb_rate,
        seed=seed,
    )

    metadata = {
        "simulation_type": "founder_mosaic_haplotype_simulation",
        "n_founders": n_founders,
        "n_reference_haplotypes": n_ref_haps,
        "n_target_haplotypes": n_target_haps,
        "n_snps": n_snps,
        "recomb_rate_per_snp_interval": recomb_rate,
        "seed": seed,
        "description": (
            "Synthetic benchmark dataset inspired by non-model population "
            "genetic settings. It is not a real white clover dataset."
        ),
    }

    save_dataset(ref_panel, true_targets, metadata)
    return ref_panel, true_targets


if __name__ == "__main__":
    ref_panel, true_targets = generate_and_save_dataset()

    print(
        f"Reference panel: {ref_panel.shape[0]} haplotypes x "
        f"{ref_panel.shape[1]} SNPs"
    )
    print(
        f"Target haplotypes (held out, ground truth): "
        f"{true_targets.shape[0]}"
    )

    demo_rng = np.random.default_rng(123)
    _, mask = simulate_low_coverage(
        true_targets, depth=1.0, rng=demo_rng
    )
    print(f"At 1x depth, fraction of sites observed: {mask.mean():.3f}")
