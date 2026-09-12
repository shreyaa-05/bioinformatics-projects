# Bioinformatics Projects — Genomic and Microbial Data Analysis

This repository contains computational bioinformatics projects developed to explore genomic imputation, microbiome and virome analysis, statistical modelling, and biological data interpretation.

The projects use simulated biological datasets to investigate practical computational problems while demonstrating the application of statistical and algorithmic methods in bioinformatics.

---

# Project 1 — Genomic Imputation Benchmark

This project investigates how **reference panel size** and **sequencing depth** influence genomic imputation accuracy.

The study simulates haplotype data, introduces missing genetic information to represent low-coverage sequencing, and evaluates how effectively the missing variants can be recovered using information from a reference panel.

The analysis progressively evaluates multiple reference panel sizes and sequencing depths and compares their effect on imputation performance.

## Technical Implementation

The computational framework is implemented in Python using numerical, data-analysis, and visualization libraries.

The main components of the project are:

- **Haplotype simulation:** Generates reference haplotypes and held-out target haplotypes representing genomic samples.
- **Reference panel:** Provides known haplotype patterns that can be used to infer missing variants in target samples.
- **Low-coverage simulation:** Randomly masks variants to simulate incomplete sequencing observations at different sequencing depths.
- **Haplotype-copying model:** Uses a simplified Li & Stephens-inspired approach to model the target haplotype as being reconstructed from reference haplotypes.
- **HMM concepts:** Hidden states represent candidate reference haplotypes, while transitions allow switching between reference haplotypes along the genome.
- **Imputation:** Missing variants are inferred from the observed target variants and reference-panel information.
- **Benchmarking:** The imputation method is evaluated across different reference panel sizes and sequencing depths.
- **Performance analysis:** Genotype concordance and dosage R² are used to quantify imputation accuracy.

## Methodology

The project follows a computational benchmarking workflow:

- **Reference haplotype generation**  
  A simulated reference panel containing known haplotypes is generated.

- **Target sample generation**  
  A separate set of target haplotypes is held out from the reference panel so that the imputation results can be compared with known ground truth.

- **Sequencing-depth simulation**  
  Different sequencing depths are simulated by varying the proportion of genetic variants that remain observed.

- **Haplotype matching and copying**  
  Target haplotypes are compared with reference haplotypes to identify the most compatible genetic patterns.

- **HMM-based state modelling**  
  The model treats the reference haplotype being copied as a hidden state and allows transitions between reference haplotypes.

- **Variant imputation**  
  Unobserved variants are estimated using the selected reference haplotype information.

- **Performance evaluation**  
  The resulting imputed variants are compared with the known simulated truth.

## Algorithms Used

- **Li & Stephens-inspired haplotype copying** — models a target haplotype as a mosaic of segments copied from reference haplotypes.
- **Hidden Markov Model (HMM) concepts** — represents the reference haplotype being copied as a hidden state and models transitions between states.
- **Probabilistic genotype imputation** — estimates missing genetic variants using observed variants and reference-panel information.
- **Genotype concordance** — measures the proportion of imputed genotype calls that agree with the true genotype.
- **Dosage R²** — measures the squared correlation between predicted allele dosage and the true allele dosage.

> **Note:** This is a simplified, from-scratch implementation inspired by the Li & Stephens haplotype-copying framework. External imputation software such as BEAGLE, GLIMPSE, and STITCH was not directly implemented.

## Technologies Used

- Python
- NumPy
- Pandas
- Matplotlib

## Experimental Setup

The benchmark evaluates:

- **Reference panel sizes:** 20, 50, 100, and 200 haplotypes
- **Sequencing depths:** 0.1×, 0.25×, 0.5×, 1×, 2×, and 4×
- **Reference panel:** 200 simulated haplotypes
- **Target samples:** 15 held-out haplotypes
- **Variants:** 800 SNP positions

## Results

The benchmark shows that larger reference panels generally provide better imputation accuracy.

With the largest reference panel of 200 haplotypes, the benchmark achieved:

- **Maximum genotype concordance:** 1.000
- **Maximum dosage R²:** 0.997

The results demonstrate the importance of reference-panel representation when recovering missing genetic information.

### Terminal Output

<img width="743" height="711" alt="Genomic_Imputation_Output" src="https://github.com/user-attachments/assets/8959b5a3-2010-41dc-80ff-ce451f4a057e" />


### Results Figures

<img width="708" height="505" alt="WhatsApp Image 2026-09-12 at 10 08 38 PM" src="https://github.com/user-attachments/assets/ab104fd0-b249-45ce-8359-ab0c706af62c" />


<img width="909" height="369" alt="WhatsApp Image 2026-09-12 at 10 10 52 PM" src="https://github.com/user-attachments/assets/422f2d35-f745-4b2f-9b73-6d59365891be" />

<img width="704" height="502" alt="WhatsApp Image 2026-09-12 at 10 09 21 PM" src="https://github.com/user-attachments/assets/fa7e3c6e-62cc-476f-9026-bbfabdbaf0b7" />

---
# Project 2 — Permafrost Virome–Microbiome Analysis

This project simulates paired **microbiome and virome abundance data** across a permafrost age chronosequence.

The objective is to investigate changes in microbial diversity and community composition across permafrost age and to identify potential bacteriophage–bacteria host relationships using abundance-based statistical analysis.

The project combines ecological diversity metrics, community ordination, correlation analysis, and quantitative evaluation of host predictions.

## Technical Implementation

The computational framework is implemented in Python using scientific computing, statistical analysis, machine-learning, and visualization libraries.

The main components of the project are:

- **Microbiome simulation:** Generates bacterial abundance profiles across simulated permafrost samples.
- **Virome simulation:** Generates phage abundance profiles linked to bacterial abundance patterns.
- **Age-gradient modelling:** Represents microbial communities across different stages of permafrost age.
- **Alpha-diversity analysis:** Calculates Shannon diversity, Simpson diversity, and species richness.
- **Beta-diversity analysis:** Uses Bray–Curtis dissimilarity to quantify differences between microbial communities.
- **Community ordination:** Uses Multidimensional Scaling (MDS) to visualize relationships between samples.
- **Virus–host association:** Uses Spearman rank correlation to identify bacterial taxa whose abundance patterns are associated with individual phages.
- **Host prediction:** Ranks candidate bacterial hosts based on correlation strength and applies a correlation threshold.
- **Validation:** Compares predicted host relationships with the simulated ground-truth host relationships.
- **Performance evaluation:** Calculates precision and recall for the host predictions.

## Methodology

The project follows a computational microbiome and virome analysis workflow:

- **Microbial community simulation**  
  Bacterial and phage abundance tables are generated for simulated permafrost samples.

- **Permafrost age modelling**  
  Samples are assigned positions along a relative permafrost age gradient to represent differences in microbial community structure.

- **Alpha-diversity analysis**  
  Shannon diversity, Simpson diversity, and species richness are calculated to describe within-sample microbial diversity.

- **Beta-diversity analysis**  
  Bray–Curtis dissimilarity is calculated between samples to quantify differences in community composition.

- **Ordination**  
  Multidimensional Scaling (MDS) is applied to the Bray–Curtis distance matrix to visualize sample-level community relationships.

- **Virus–host association analysis**  
  Spearman rank correlations are calculated between each phage and bacterial taxon across samples.

- **Host prediction**  
  The bacterial taxon with the strongest qualifying positive correlation is selected as the predicted host.

- **Validation**  
  Predicted host relationships are compared against the simulated true host relationships.

- **Performance evaluation**  
  Precision and recall are calculated to measure the quality of the predicted virus–host relationships.

## Algorithms Used

- **Shannon diversity** — measures microbial diversity by considering both richness and relative abundance.
- **Simpson diversity** — measures community diversity with greater emphasis on dominant taxa.
- **Species richness** — counts the number of observed taxa in a sample.
- **Bray–Curtis dissimilarity** — measures differences in community composition between samples based on abundance.
- **Multidimensional Scaling (MDS)** — projects sample-to-sample dissimilarities into a lower-dimensional space for visualization.
- **Spearman rank correlation** — measures monotonic relationships between phage and bacterial abundance patterns.
- **Co-abundance-based host prediction** — identifies potential phage hosts by ranking bacteria according to their abundance correlation with each phage.
- **Precision and recall** — evaluate the accuracy and coverage of the predicted host relationships.

## Technologies Used

- Python
- NumPy
- Pandas
- SciPy
- scikit-learn
- Matplotlib

## Experimental Setup

The simulated dataset contains:

- **24 samples**
- **40 bacterial taxa**
- **25 phage taxa**
- A simulated permafrost age gradient
- Known simulated phage–bacteria host relationships for validation

## Results

The virus–host analysis identified potential host relationships for **23 out of 25 phages** using a Spearman correlation threshold of ρ ≥ 0.5.

The predictions achieved:

- **Precision:** 0.96
- **Recall:** 0.88

These results demonstrate how abundance-based statistical relationships can be used to identify potential virus–host associations in simulated microbial communities.

### Terminal Output


<img width="902" height="184" alt="Permafrost_Virome_Output" src="https://github.com/user-attachments/assets/06e69868-faa5-4daa-adaa-f3c247974215" />


### Results Figures
<img width="909" height="369" alt="WhatsApp Image 2026-09-12 at 10 10 52 PM" src="https://github.com/user-attachments/assets/3abd4a29-1ea5-4cca-849f-b90faac66a8c" />


<img width="909" height="368" alt="WhatsApp Image 2026-09-12 at 10 11 39 PM" src="https://github.com/user-attachments/assets/d5ef5fd6-21f5-4665-bb35-78b7f9a2e5fb" />


<img width="672" height="549" alt="WhatsApp Image 2026-09-12 at 10 12 40 PM" src="https://github.com/user-attachments/assets/f7cf3a35-86fa-4327-be7a-f3da5f1ab4d6" />


