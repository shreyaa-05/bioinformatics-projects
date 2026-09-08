"""
simulate_metagenome.py
------------------------
Simulates paired microbiome + virome relative-abundance tables across a
permafrost age chronosequence (older, deeper Pleistocene-age samples through
to younger, shallower Holocene-age / actively-thawing samples), for
prototyping the analysis pipeline described in the project brief:
"identify changes in the microbiomes and viromes of ancient permafrost
samples spanning a chronosequence... and link viruses to their host
microorganisms."

No real metagenomic sequencing data was available at prototyping time, so
this generates synthetic but ecologically-structured community data with
three deliberately-designed signals for the pipeline to recover:

  1. Directional community turnover with age (older samples dominated by a
     "dormant/cryophile" taxon set; younger/thawing samples show blooms of a
     "fast-responder" taxon set), consistent with reactivation of dormant
     cells on thaw.
  2. A declining total-activity trend with age (older samples = lower
     effective community abundance/evenness), since ancient permafrost
     communities are less metabolically active than samples closer to the
     thaw front.
  3. Kill-the-winner viral dynamics: each simulated phage OTU is coupled to
     one specific bacterial host OTU, with phage log-abundance tracking a
     noisy linear function of host log-abundance in the same sample. This
     gives the virus-host linkage step in the analysis script a real,
     recoverable signal to detect via co-abundance correlation.
"""

import numpy as np
import pandas as pd


def simulate_chronosequence(n_samples=24, n_bacterial_taxa=40, n_phage_taxa=25,
                             seed=7):
    rng = np.random.default_rng(seed)

    # Sample ages: 0 = present-day/actively thawing, 1 = oldest (deep
    # Pleistocene). Evenly spaced along the chronosequence.
    age = np.linspace(0, 1, n_samples)

    # --- Bacterial community ---
    # Split taxa into "dormant/cryophile" (favoured at high age) and
    # "fast-responder" (favoured at low age / post-thaw) guilds.
    n_dormant = n_bacterial_taxa // 2
    n_fast = n_bacterial_taxa - n_dormant

    bact_log_abund = np.zeros((n_samples, n_bacterial_taxa))
    baseline = rng.normal(0, 0.5, n_bacterial_taxa)  # taxon-level baseline log-abundance

    for i, a in enumerate(age):
        effect = np.zeros(n_bacterial_taxa)
        effect[:n_dormant] = 2.0 * a          # dormant guild rises with age
        effect[n_dormant:] = 2.0 * (1 - a)    # fast responders rise as age -> 0 (thaw)
        noise = rng.normal(0, 0.4, n_bacterial_taxa)
        bact_log_abund[i] = baseline + effect + noise

    # Overall community activity/evenness declines with age: shrink the
    # dynamic range and add a floor of near-dormant reads for old samples.
    activity_scale = 1.0 - 0.5 * age  # 1.0 at present -> 0.5 at oldest
    bact_log_abund *= activity_scale[:, None]

    bact_abund = np.exp(bact_log_abund)
    bact_rel = bact_abund / bact_abund.sum(axis=1, keepdims=True)

    bact_taxa_names = (
        [f"Bacterium_dormant_{i+1:02d}" for i in range(n_dormant)] +
        [f"Bacterium_fastresp_{i+1:02d}" for i in range(n_fast)]
    )

    # --- Viral (phage) community, coupled to specific bacterial hosts ---
    # Each phage is assigned one host taxon (biased toward fast-responder
    # hosts, since those guilds bloom and are more likely to sustain lytic
    # infection dynamics near the thaw front).
    host_idx = rng.choice(n_bacterial_taxa, size=n_phage_taxa, replace=True,
                           p=np.concatenate([
                               np.full(n_dormant, 0.3 / n_dormant),
                               np.full(n_fast, 0.7 / n_fast),
                           ]))

    phage_log_abund = np.zeros((n_samples, n_phage_taxa))
    for p in range(n_phage_taxa):
        host = host_idx[p]
        host_series = bact_log_abund[:, host]
        # Phage tracks host abundance in the *same* sample (no lag): in a
        # cross-sectional chronosequence design, kill-the-winner dynamics
        # play out on timescales far shorter than the spacing between
        # samples, so phage and host co-abundance is expected to be
        # synchronous at the sampling resolution used here. A stronger
        # coupling and lower noise than the exploratory first pass makes the
        # co-abundance signal detectable against 24 samples.
        coupling = rng.uniform(0.85, 1.05)
        phage_baseline = rng.normal(-1, 0.3)
        noise = rng.normal(0, 0.25, n_samples)
        phage_log_abund[:, p] = phage_baseline + coupling * host_series + noise

    phage_abund = np.exp(phage_log_abund)
    phage_rel = phage_abund / phage_abund.sum(axis=1, keepdims=True)
    phage_taxa_names = [f"Phage_{i+1:02d}" for i in range(n_phage_taxa)]

    sample_names = [f"S{i+1:02d}_age{a:.2f}" for i, a in enumerate(age)]

    bact_df = pd.DataFrame(bact_rel, index=sample_names, columns=bact_taxa_names)
    phage_df = pd.DataFrame(phage_rel, index=sample_names, columns=phage_taxa_names)

    metadata = pd.DataFrame({
        "sample": sample_names,
        "relative_age": age,  # 0 = present/actively thawing, 1 = oldest (deep Pleistocene)
        "era": np.where(age < 0.4, "Holocene/thaw-front",
                 np.where(age < 0.75, "Late Pleistocene", "Deep Pleistocene")),
    }).set_index("sample")

    true_host_links = pd.DataFrame({
        "phage": phage_taxa_names,
        "true_host": [bact_taxa_names[h] for h in host_idx],
    })

    return bact_df, phage_df, metadata, true_host_links


if __name__ == "__main__":
    bact_df, phage_df, metadata, true_links = simulate_chronosequence()
    bact_df.to_csv("results/microbiome_abundance.csv")
    phage_df.to_csv("results/virome_abundance.csv")
    metadata.to_csv("results/sample_metadata.csv")
    true_links.to_csv("results/true_host_links_for_validation.csv", index=False)
    print(f"Simulated {bact_df.shape[0]} samples, {bact_df.shape[1]} bacterial taxa, "
          f"{phage_df.shape[1]} phage taxa.")
    print("Saved abundance tables + metadata + ground-truth links to results/")
