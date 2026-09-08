# MITACS Portfolio — Bioinformatics Demo Projects

Two self-contained, runnable demo pipelines built to demonstrate bioinformatics and
computational biology skills relevant to the two project briefs below. Both projects use
**simulated data with realistic, deliberately-embedded biological structure**, since no
real datasets from either study were available at prototyping time — every script is
fully documented on exactly what is simulated, why, and what the corresponding real
analysis would additionally require.

## 1. `project1_genomic_imputation/`
*Comparing and benchmarking genomic imputation algorithms in a non-model system*

A from-scratch NumPy implementation of the Li & Stephens haplotype-copying HMM (the
statistical core of BEAGLE/GLIMPSE/STITCH), benchmarked for imputation accuracy across a
grid of reference panel sizes and sequencing depths, using a coalescent-style simulator
parameterized for low recombination/diversity (as in outbreeding non-model species like
white clover).

**Skills demonstrated:** population/statistical genetics, HMM implementation, simulation
design, benchmarking methodology, scientific visualization.

## 2. `project2_permafrost_virome/`
*Mapping viromes to microbiomes in thawing permafrost*

A simulated microbiome + virome chronosequence (Pleistocene → Holocene/thaw front) with
embedded community turnover and kill-the-winner phage-host dynamics, analyzed with alpha/beta
diversity metrics, Bray-Curtis ordination, and co-abundance-based virus-host linkage
prediction (validated against ground truth).

**Skills demonstrated:** microbial ecology/metagenomics analysis, community ecology
statistics, network/linkage inference, scientific visualization.

## Running everything

Each project folder is independent. From inside either folder:

```bash
pip install numpy pandas matplotlib scipy scikit-learn
python3 <simulate_script>.py
python3 <analysis_script>.py
# project 1 only:
python3 plot_results.py
```

All outputs (CSVs + figures) land in that project's `results/` subfolder and are already
included pre-generated in this bundle.
