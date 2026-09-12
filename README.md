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

<img width="708" height="505" alt="WhatsApp Image 2026-09-12 at 10 08 38 PM" src="https://github.com/user-attachments/assets/ab104fd0-b249-45ce-8359-ab0c706af62c" />

<img width="909" height="369" alt="WhatsApp Image 2026-09-12 at 10 10 52 PM" src="https://github.com/user-attachments/assets/422f2d35-f745-4b2f-9b73-6d59365891be" />
<img width="704" height="502" alt="WhatsApp Image 2026-09-12 at 10 09 21 PM" src="https://github.com/user-attachments/assets/fa7e3c6e-62cc-476f-9026-bbfabdbaf0b7" />



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

#### Terminal Output

<img width="902" height="184" alt="Permafrost_Virome_Output" src="https://github.com/user-attachments/assets/06e69868-faa5-4daa-adaa-f3c247974215" />

<img width="909" height="369" alt="WhatsApp Image 2026-09-12 at 10 10 52 PM" src="https://github.com/user-attachments/assets/3abd4a29-1ea5-4cca-849f-b90faac66a8c" />



<img width="909" height="368" alt="WhatsApp Image 2026-09-12 at 10 11 39 PM" src="https://github.com/user-attachments/assets/d5ef5fd6-21f5-4665-bb35-78b7f9a2e5fb" />




<img width="672" height="549" alt="WhatsApp Image 2026-09-12 at 10 12 40 PM" src="https://github.com/user-attachments/assets/f7cf3a35-86fa-4327-be7a-f3da5f1ab4d6" />


