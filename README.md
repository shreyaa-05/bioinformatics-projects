# Bioinformatics Projects

A collection of computational bioinformatics projects exploring genomic data analysis, statistical modelling, microbiome/virome analysis, and biological data interpretation.

## Projects

### 1. Genomic Imputation Benchmark

**Folder:** `bioinformatics/module1_genomic_imputation/`

This project investigates how **reference panel size** and **sequencing depth** affect genomic imputation accuracy.

#### Objective
To simulate haplotype data and evaluate how different reference panel sizes and sequencing depths influence the accuracy of recovering missing genetic variants.

#### Methodology
- Simulated reference haplotypes and held-out target haplotypes
- Modelled low-coverage sequencing by randomly masking variants
- Used a simplified **Li & Stephens-inspired haplotype-copying approach**
- Evaluated multiple reference panel sizes and sequencing depths
- Measured:
  - Genotype concordance
  - Dosage R²
- Generated performance plots and heatmaps

#### Technologies
- Python
- NumPy
- Pandas
- Matplotlib

#### Key Concepts
- Haplotype modelling
- Reference panels
- Genomic imputation
- Hidden-state modelling
- Haplotype copying
- Genotype concordance
- Dosage R²

#### Terminal Output



<img width="743" height="711" alt="Genomic_Imputation_Output" src="https://github.com/user-attachments/assets/8959b5a3-2010-41dc-80ff-ce451f4a057e" />


---

### 2. Permafrost Virome–Microbiome Analysis

**Folder:** `bioinformatics/module2_permafrost_virome/`

This project simulates paired **microbiome and virome abundance data** across a permafrost age chronosequence and investigates microbial diversity, community composition, and virus–host relationships.

#### Objective
To explore how microbial communities may change across permafrost age and to predict potential bacteriophage–bacteria host relationships from abundance patterns.

#### Methodology
- Simulated bacterial and phage abundance tables
- Modelled samples across a permafrost age gradient
- Calculated:
  - Shannon diversity
  - Simpson diversity
  - Species richness
- Used **Bray–Curtis dissimilarity** to measure community differences
- Applied **multidimensional scaling (MDS)** for community ordination
- Used **Spearman correlation** to identify potential virus–host relationships
- Evaluated host prediction using:
  - Precision
  - Recall
- Generated correlation heatmaps and diversity plots

#### Technologies
- Python
- NumPy
- Pandas
- SciPy
- scikit-learn
- Matplotlib

#### Key Concepts
- Microbiome analysis
- Virome analysis
- Alpha diversity
- Beta diversity
- Bray–Curtis dissimilarity
- Ordination
- Virus–host prediction
- Co-abundance analysis
- Statistical evaluation

---

## Repository Structure

```text
bioinformatics-projects/
│
├── README.md
│
└── bioinformatics/
    │
    ├── README.md
    │
    ├── module1_genomic_imputation/
    │   ├── data/
    │   ├── results/
    │   ├── simulate_data.py
    │   ├── run_imputation_benchmark.py
    │   ├── plot_results.py
    │   └── README.md
    │
    └── module2_permafrost_virome/
        ├── results/
        ├── simulate_metagenome.py
        ├── analyze_virome_microbiome.py
        └── README.md
