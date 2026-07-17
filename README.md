# S_brain
Single brain MR



# Step 1: Extract Significant eQTL Instruments for Mendelian Randomization

## Why this step?

The SingleBrain `*_full_assoc.tsv.gz` files are approximately **4.5 GB** each. Loading these directly into R is inefficient and will likely exceed available memory.

For Mendelian Randomization (MR), we only need **genome-wide significant eQTLs** to serve as genetic instruments.

Additionally:

- Merge with the GWAS data will be performed using **rsID**
- Do **not** merge by chromosome/position because:
  - PGC GWAS datasets are generally **GRCh37 (hg19)**
  - SingleBrain eQTL data are **GRCh38 (hg38)**
- SingleBrain uses uppercase alleles while some GWAS files use lowercase alleles. This can be standardized later during harmonization.

---

# Move to the SingleBrain directory

```bash
cd /scratch/negishi/hasan128/S_brain/SingleBrain_eQTL_data/primary
```

---

# Verify the column order

Before filtering, inspect the header.

```bash
zcat Ast_eqtl_full_assoc.tsv.gz | head -1 | tr '\t' '\n' | nl
```

Confirm that:

- Column **29** = `Random_FDR`

If the column differs, adjust the `awk` command accordingly.

---

# Filter Astrocyte eQTLs

Keep:

- Header
- Variants with

```
Random_FDR < 0.05
```

Command:

```bash
zcat Ast_eqtl_full_assoc.tsv.gz \
| awk -F'\t' 'NR==1 || ($29 != "NA" && $29 < 0.05)' \
| gzip > Ast_eqtl_sig_instruments.tsv.gz
```

---

# Verify the output

Check the number of retained variants.

```bash
zcat Ast_eqtl_sig_instruments.tsv.gz | wc -l
```

Inspect the first few rows.

```bash
zcat Ast_eqtl_sig_instruments.tsv.gz | head
```

Check the file size.

```bash
ls -lh Ast_eqtl_sig_instruments.tsv.gz
```

---

# Process all primary cell types

Once the Astrocyte file works correctly, run the same filtering for every primary cell type.

```bash
for file in *_eqtl_full_assoc.tsv.gz
do
    echo "Processing $file..."

    prefix="${file%_full_assoc.tsv.gz}"

    zcat "$file" \
    | awk -F'\t' 'NR==1 || ($29 != "NA" && $29 < 0.05)' \
    | gzip > "${prefix}_sig_instruments.tsv.gz"

done
```

---

# Expected output files

```
Ast_eqtl_sig_instruments.tsv.gz
End_eqtl_sig_instruments.tsv.gz
Ext_eqtl_sig_instruments.tsv.gz
IN_eqtl_sig_instruments.tsv.gz
MG_eqtl_sig_instruments.tsv.gz
OD_eqtl_sig_instruments.tsv.gz
OPC_eqtl_sig_instruments.tsv.gz
```

---

# Verify all outputs

```bash
ls -lh *_sig_instruments.tsv.gz
```

Count the retained variants in each file.

```bash
for file in *_sig_instruments.tsv.gz
do
    echo "$file"
    zcat "$file" | wc -l
done
```

---

# Why this filtering?

The original files contain association statistics for **every tested SNP–gene pair**.

MR only requires **strong genetic instruments**, so retaining only significant eQTLs:

- greatly reduces file size
- reduces memory usage
- speeds up downstream LD clumping
- speeds up harmonization with GWAS
- makes MR analyses substantially more efficient

---

# Notes

- Merge using **rsID** (`variant_id` in SingleBrain and `SNP` in PGC GWAS).
- Do **not** merge using chromosome or genomic position because the datasets use different genome builds (GRCh38 vs. GRCh37).
- Allele case (uppercase/lowercase) will be standardized during harmonization.

---

# Next Step

After generating the significant instrument files, the next step will be:

**Step 2 — LD Clumping of Significant eQTL Instruments**

This ensures the retained SNPs are independent genetic instruments before Mendelian Randomization.











# Step 2: LD Clumping of Significant eQTL Instruments

## Objective

After extracting genome-wide significant eQTLs, the next step is to identify **independent genetic instruments** for Mendelian Randomization (MR).

LD clumping removes highly correlated SNPs so that only approximately independent variants remain.

---

## Why LD Clumping?

Mendelian Randomization assumes that instrumental variables are independent.

Using correlated SNPs can:

- Inflate MR estimates
- Bias causal effect estimates
- Violate instrumental variable assumptions

A commonly used LD threshold for MR is:

| Parameter | Value |
|-----------|------:|
| LD reference | 1000 Genomes European |
| `r²` threshold | **0.01** |
| Window size | **1000 kb** |

---

# Step 1: Request an Interactive Compute Node

Large PLINK jobs should **not** be run on the login node.

Request an interactive compute node:

```bash
salloc \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=4 \
    --mem=16G \
    --time=02:00:00 \
    --partition=shared
```

Wait until your prompt changes from something like

```
login05
```

to

```
a001
```

or another compute node name.

---

# Step 2: Move to the Working Directory

```bash
cd /scratch/negishi/hasan128/S_brain/SingleBrain_eQTL_data/primary
```

---

# Step 3: Decompress the Filtered eQTL File

PLINK requires an uncompressed text file.

```bash
zcat Ast_eqtl_sig_instruments.tsv.gz > Ast_working.txt
```

Verify that the file was created.

```bash
ls -lh Ast_working.txt
```

Inspect the first few lines.

```bash
head Ast_working.txt
```

---

# Step 4: Run LD Clumping with PLINK

Use the filtered Astrocyte eQTL file as input.

```bash
plink \
    --bfile /depot/ppaschou/data/NEW_CKM/qc/g1000_eur \
    --clump Ast_working.txt \
    --clump-snp-field variant_id \
    --clump-field Random_P \
    --clump-p1 0.05 \
    --clump-r2 0.01 \
    --clump-kb 1000 \
    --out Ast_EUR_clumped
```

### Parameter Explanation

| Option | Description |
|---------|-------------|
| `--bfile` | 1000 Genomes European reference panel |
| `--clump` | Input association file |
| `--clump-snp-field variant_id` | SNP identifier column |
| `--clump-field Random_P` | P-value column used for clumping |
| `--clump-p1 0.05` | Include SNPs with P ≤ 0.05 |
| `--clump-r2 0.01` | LD threshold |
| `--clump-kb 1000` | 1 Mb clumping window |
| `--out` | Output file prefix |

---

# Step 5: Verify the Results

PLINK should generate several output files.

```bash
ls -lh Ast_EUR_clumped*
```

The most important output is:

```
Ast_EUR_clumped.clumped
```

Inspect the first few rows.

```bash
head Ast_EUR_clumped.clumped
```

Count the number of independent instruments.

```bash
wc -l Ast_EUR_clumped.clumped
```

---

# Expected Output

```
Ast_EUR_clumped.clumped
Ast_EUR_clumped.log
Ast_EUR_clumped.nosex
```

The `.clumped` file contains the **independent SNP instruments** that will be used for downstream Mendelian Randomization analyses.

---

# Next Step

Repeat LD clumping for the remaining primary brain cell types:

- Astrocytes (Ast)
- Endothelial cells (End)
- Excitatory neurons (Ext)
- Inhibitory neurons (IN)
- Microglia (MG)
- Oligodendrocytes (OD)
- Oligodendrocyte precursor cells (OPC)

Each will produce its own set of independent eQTL instruments.

---

# Next Analysis

After LD clumping is complete for all cell types, proceed to:

**Step 3 — Harmonize the independent eQTL instruments with the PGC GWAS summary statistics using rsID before performing Mendelian Randomization.**







# Step 3: Create the Final Exposure Dataset for Mendelian Randomization

## Objective

After LD clumping, PLINK identifies the **lead (independent) SNPs** within each LD block. However, the `.clumped` file only contains summary information about these SNPs.

To prepare the exposure dataset for Mendelian Randomization (MR), we need to:

1. Extract the list of independent SNPs from the PLINK output.
2. Match those SNPs back to the original eQTL association results.
3. Retain all required MR information, including:
   - SNP ID
   - Gene
   - Effect allele
   - Other allele
   - Effect size (β)
   - Standard error (SE)
   - P-value
   - Additional annotation

The final output will serve as the **exposure dataset** for MR.

---

## Step 1: Extract Independent SNP IDs

PLINK stores the lead SNP for each clump in the **3rd column** of the `.clumped` file.

Extract these SNP IDs:

```bash
awk 'NR>1 && NF>0 {print $3}' Ast_EUR_clumped.clumped > Ast_clumped_snps_list.txt
```

---

## Verify the SNP List

Count the number of independent SNPs.

```bash
wc -l Ast_clumped_snps_list.txt
```

Inspect the first few SNPs.

```bash
head Ast_clumped_snps_list.txt
```

---

## Step 2: Create the Final MR Exposure Dataset

Filter the original significant eQTL file so that only the independent SNPs remain.

```bash
awk '
NR==FNR {
    keep[$1]
    next
}
FNR==1 || ($2 in keep)
' Ast_clumped_snps_list.txt Ast_working.txt > Ast_MR_exposure.tsv
```

This command:

- reads the independent SNP list
- scans the original filtered eQTL file
- retains:
  - the header
  - only rows whose `variant_id` appears in the clumped SNP list

---

## Verify the Exposure Dataset

Check the number of rows.

```bash
wc -l Ast_MR_exposure.tsv
```

Inspect the first few lines.

```bash
head Ast_MR_exposure.tsv
```

Check the file size.

```bash
ls -lh Ast_MR_exposure.tsv
```

---

## Expected Output

```
Ast_MR_exposure.tsv
```

This file contains the independent eQTL instruments together with all required association statistics for Mendelian Randomization.

Typical retained information includes:

- Gene (feature)
- SNP (variant_id)
- Chromosome
- Position
- Reference allele
- Alternate allele
- Effect allele
- Effect size (β)
- Standard error (SE)
- Z statistic
- Random-effects P-value
- Random-effects FDR

---

## Optional Cleanup

After confirming that the exposure dataset was created successfully, temporary files can be removed to save storage.

```bash
rm Ast_working.txt
rm Ast_clumped_snps_list.txt
```

---

## Files Generated So Far

```
Ast_eqtl_sig_instruments.tsv.gz
Ast_working.txt                (optional temporary file)
Ast_EUR_clumped.clumped
Ast_EUR_clumped.log
Ast_clumped_snps_list.txt      (optional temporary file)
Ast_MR_exposure.tsv            ← Final exposure dataset
```

---

## Pipeline Progress

```text
SingleBrain Full eQTL
          │
          ▼
Filter Significant eQTLs (FDR < 0.05)
          │
          ▼
LD Clumping (PLINK)
          │
          ▼
Extract Independent SNPs
          │
          ▼
Create Final MR Exposure Dataset
          │
          ▼
Ready for MR Harmonization
```

---

# Next Step

The next stage is to prepare the **outcome dataset** by loading one of the psychiatric GWAS summary statistics (e.g., ADHD, schizophrenia, bipolar disorder, depression) and harmonizing it with the exposure dataset using shared **rsIDs** before performing Mendelian Randomization.
































# Astrocyte–Schizophrenia Mendelian Randomization and Colocalization Pipeline

## Overview

This workflow performs the following analyses:

1. Formats astrocyte-specific SingleBrain eQTL data as the exposure dataset
2. Formats schizophrenia GWAS summary statistics as the outcome dataset
3. Harmonizes exposure and outcome alleles
4. Runs Mendelian Randomization using `TwoSampleMR`
5. Applies False Discovery Rate correction
6. Maps Ensembl gene IDs to HGNC gene symbols
7. Runs gene-level colocalization using `coloc`
8. Produces a final prioritized target file containing MR and colocalization results

---

# Directory Structure

The workflow assumes the following directory structure:

```text
/scratch/negishi/hasan128/S_brain/
├── Ast_SCZ_Master_Prioritized_Targets.csv
├── Ast_SCZ_MR_results.csv
├── run_Ast_SCZ_MR.R
├── run_automated_coloc.R
├── submit_mr.sh
├── logs/
├── psychiatric/
│   └── scz.tsv
└── SingleBrain_eQTL_data/
    └── primary/
        ├── Ast_MR_exposure.tsv
        └── Ast_eqtl_full_assoc.tsv.gz
```

---

# Part 1: Astrocyte–Schizophrenia Mendelian Randomization

## R Script

Create the following file:

```text
/scratch/negishi/hasan128/S_brain/run_Ast_SCZ_MR.R
```

Add the following code:

```r
library(TwoSampleMR)
library(data.table)
library(dplyr)

# ============================================================
# 1. File paths
# ============================================================

base_dir <- "/scratch/negishi/hasan128/S_brain"

exposure_file <- file.path(
  base_dir,
  "SingleBrain_eQTL_data",
  "primary",
  "Ast_MR_exposure.tsv"
)

outcome_file <- file.path(
  base_dir,
  "psychiatric",
  "scz.tsv"
)

output_file <- file.path(
  base_dir,
  "Ast_SCZ_MR_results.csv"
)

# ============================================================
# 2. Load and clean the astrocyte exposure data
# ============================================================

message("Loading astrocyte eQTL exposure data...")

exp_raw <- fread(exposure_file)

required_exposure_columns <- c(
  "feature",
  "variant_id",
  "Allele",
  "ref",
  "alt",
  "fixed_beta",
  "fixed_sd",
  "Random_P",
  "chr",
  "pos"
)

missing_exposure_columns <- setdiff(
  required_exposure_columns,
  names(exp_raw)
)

if (length(missing_exposure_columns) > 0) {
  stop(
    paste(
      "Missing exposure columns:",
      paste(missing_exposure_columns, collapse = ", ")
    )
  )
}

# Standardize allele formatting
exp_raw[, Allele := toupper(Allele)]
exp_raw[, ref := toupper(ref)]
exp_raw[, alt := toupper(alt)]

# Derive the non-effect allele
exp_raw[, other_allele := fifelse(
  Allele == alt,
  ref,
  alt
)]

# Remove invalid and duplicate SNP-gene associations
exp_raw <- exp_raw[
  !is.na(feature) &
  !is.na(variant_id) &
  !is.na(fixed_beta) &
  !is.na(fixed_sd) &
  fixed_sd > 0 &
  !is.na(Random_P)
]

exp_raw <- unique(
  exp_raw,
  by = c("feature", "variant_id")
)

message(
  "Exposure rows retained: ",
  format(nrow(exp_raw), big.mark = ",")
)

# ============================================================
# 3. Format the exposure dataset
# ============================================================

ast_exp <- format_data(
  dat = as.data.frame(exp_raw),
  type = "exposure",
  snp_col = "variant_id",
  phenotype_col = "feature",
  effect_allele_col = "Allele",
  other_allele_col = "other_allele",
  beta_col = "fixed_beta",
  se_col = "fixed_sd",
  pval_col = "Random_P",
  chr_col = "chr",
  pos_col = "pos"
)

# ============================================================
# 4. Load and filter the schizophrenia GWAS
# ============================================================

message("Loading schizophrenia GWAS data...")

scz_raw <- fread(outcome_file)

required_outcome_columns <- c(
  "SNP",
  "BETA",
  "SE",
  "A1",
  "A2",
  "P",
  "CHR",
  "POS"
)

missing_outcome_columns <- setdiff(
  required_outcome_columns,
  names(scz_raw)
)

if (length(missing_outcome_columns) > 0) {
  stop(
    paste(
      "Missing outcome columns:",
      paste(missing_outcome_columns, collapse = ", ")
    )
  )
}

# Standardize allele formatting
scz_raw[, A1 := toupper(A1)]
scz_raw[, A2 := toupper(A2)]

# Keep only SNPs present in the exposure dataset
scz_filtered <- scz_raw[
  SNP %in% ast_exp$SNP &
  !is.na(BETA) &
  !is.na(SE) &
  SE > 0 &
  !is.na(P)
]

scz_filtered <- unique(
  scz_filtered,
  by = "SNP"
)

message(
  "Outcome SNPs retained: ",
  format(nrow(scz_filtered), big.mark = ",")
)

# ============================================================
# 5. Format the outcome dataset
# ============================================================

scz_out <- format_data(
  dat = as.data.frame(scz_filtered),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "POS"
)

# ============================================================
# 6. Harmonize exposure and outcome data
# ============================================================

message("Harmonizing exposure and outcome data...")

dat_harmonised <- harmonise_data(
  exposure_dat = ast_exp,
  outcome_dat = scz_out,
  action = 2
)

dat_harmonised <- dat_harmonised[
  dat_harmonised$mr_keep,
]

message(
  "Harmonized rows retained for MR: ",
  format(nrow(dat_harmonised), big.mark = ",")
)

if (nrow(dat_harmonised) == 0) {
  stop("No variants remained after harmonization.")
}

# ============================================================
# 7. Run Mendelian Randomization
# ============================================================

message("Running Mendelian Randomization...")

mr_results <- mr(dat_harmonised)

# ============================================================
# 8. Save the results
# ============================================================

fwrite(
  mr_results,
  output_file
)

message("MR results saved to: ", output_file)
```

---

# Part 2: SLURM Submission Script

Create the following file:

```text
/scratch/negishi/hasan128/S_brain/submit_mr.sh
```

Add the following content:

```bash
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -A pdrineas
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=Ast_SCZ_MR
#SBATCH --output=/scratch/negishi/hasan128/S_brain/logs/Ast_SCZ_MR_%j.out
#SBATCH --error=/scratch/negishi/hasan128/S_brain/logs/Ast_SCZ_MR_%j.err

set -euo pipefail

BASE_DIR="/scratch/negishi/hasan128/S_brain"
R_SCRIPT="${BASE_DIR}/run_Ast_SCZ_MR.R"

mkdir -p "${BASE_DIR}/logs"

module --force purge
module load biocontainers
module load gcc
module load r

cd "${BASE_DIR}"

echo "Job started: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "R version:"
R --version

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

Rscript "${R_SCRIPT}"

echo "Job completed: $(date)"
```

Make the file executable:

```bash
chmod +x submit_mr.sh
```

Submit the job:

```bash
cd /scratch/negishi/hasan128/S_brain
sbatch submit_mr.sh
```

Check the queue:

```bash
squeue -u hasan128
```

Inspect the logs:

```bash
ls -lh /scratch/negishi/hasan128/S_brain/logs
```

```bash
cat /scratch/negishi/hasan128/S_brain/logs/Ast_SCZ_MR_<JOB_ID>.out
```

```bash
cat /scratch/negishi/hasan128/S_brain/logs/Ast_SCZ_MR_<JOB_ID>.err
```

---

# Part 3: Install `TwoSampleMR`

Open an interactive R session:

```bash
module load r
R
```

Install `remotes`:

```r
install.packages(
  "remotes",
  repos = "https://cloud.r-project.org"
)
```

Install `TwoSampleMR` from GitHub:

```r
remotes::install_github(
  "MRCIEU/TwoSampleMR",
  upgrade = "never"
)
```

Verify the installation:

```r
library(TwoSampleMR)
packageVersion("TwoSampleMR")
```

Exit R:

```r
q()
```

---

# Part 4: Filter Significant MR Results

Create a file named:

```text
filter_Ast_SCZ_MR_results.R
```

Add the following code:

```r
library(data.table)
library(dplyr)

base_dir <- "/scratch/negishi/hasan128/S_brain"

input_file <- file.path(
  base_dir,
  "Ast_SCZ_MR_results.csv"
)

output_file <- file.path(
  base_dir,
  "Ast_SCZ_Significant_Targets.csv"
)

results <- fread(input_file)

primary_results <- results %>%
  filter(
    method %in% c(
      "Inverse variance weighted",
      "Wald ratio"
    )
  )

primary_results <- primary_results %>%
  mutate(
    fdr = p.adjust(
      pval,
      method = "fdr"
    )
  ) %>%
  arrange(fdr)

significant_hits <- primary_results %>%
  filter(fdr < 0.05)

message(
  "Significant targets at FDR < 0.05: ",
  nrow(significant_hits)
)

print(head(significant_hits, 20))

fwrite(
  significant_hits,
  output_file
)

message(
  "Significant results saved to: ",
  output_file
)
```

Run it:

```bash
module load r
Rscript filter_Ast_SCZ_MR_results.R
```

---

# Part 5: Map Ensembl IDs to HGNC Gene Symbols

Create the file:

```text
annotate_Ast_SCZ_targets.R
```

Add the following code:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install(
    "biomaRt",
    ask = FALSE,
    update = FALSE
  )
}

library(biomaRt)
library(data.table)
library(dplyr)

base_dir <- "/scratch/negishi/hasan128/S_brain"

input_file <- file.path(
  base_dir,
  "Ast_SCZ_Significant_Targets.csv"
)

output_file <- file.path(
  base_dir,
  "Ast_SCZ_Master_Prioritized_Targets.csv"
)

significant_hits <- fread(input_file)

significant_hits[, ensembl_clean := sub(
  "\\..*$",
  "",
  exposure
)]

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

gene_annotation <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = unique(significant_hits$ensembl_clean),
  mart = mart
)

gene_annotation <- unique(
  as.data.table(gene_annotation),
  by = "ensembl_gene_id"
)

annotated_results <- merge(
  significant_hits,
  gene_annotation,
  by.x = "ensembl_clean",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

annotated_results <- annotated_results %>%
  arrange(fdr)

fwrite(
  annotated_results,
  output_file
)

print(
  head(
    annotated_results[
      ,
      .(
        exposure,
        hgnc_symbol,
        b,
        pval,
        fdr,
        description
      )
    ],
    20
  )
)

message(
  "Annotated target file saved to: ",
  output_file
)
```

Run it:

```bash
module load r
Rscript annotate_Ast_SCZ_targets.R
```

---

# Part 6: Install `coloc`

Open R:

```bash
module load r
R
```

Install the package:

```r
install.packages(
  "coloc",
  repos = "https://cloud.r-project.org"
)
```

Verify:

```r
library(coloc)
packageVersion("coloc")
```

Exit R:

```r
q()
```

---

# Part 7: Automated Colocalization Pipeline

Create the following file:

```text
/scratch/negishi/hasan128/S_brain/run_automated_coloc.R
```

Add the following code:

```r
library(data.table)
library(dplyr)
library(coloc)

# ============================================================
# 1. Configuration
# ============================================================

base_dir <- "/scratch/negishi/hasan128/S_brain"

targets_file <- file.path(
  base_dir,
  "Ast_SCZ_Master_Prioritized_Targets.csv"
)

gwas_file <- file.path(
  base_dir,
  "psychiatric",
  "scz.tsv"
)

eqtl_file <- file.path(
  base_dir,
  "SingleBrain_eQTL_data",
  "primary",
  "Ast_eqtl_full_assoc.tsv.gz"
)

output_file <- file.path(
  base_dir,
  "Ast_SCZ_Final_Coloc_Validated_Targets.csv"
)

singlebrain_n <- 983

# Replace this value with the verified schizophrenia case fraction
# when exact case and control counts are available.
scz_case_fraction <- 0.35

# ============================================================
# 2. Load prioritized targets
# ============================================================

message("Loading prioritized MR targets...")

targets <- fread(targets_file)

if (!"exposure" %in% names(targets)) {
  stop("The target file does not contain an 'exposure' column.")
}

target_genes <- unique(
  targets$exposure[
    !is.na(targets$exposure)
  ]
)

message(
  "Target genes loaded: ",
  length(target_genes)
)

# ============================================================
# 3. Load schizophrenia GWAS data
# ============================================================

message("Loading schizophrenia GWAS...")

gwas <- fread(
  gwas_file,
  select = c(
    "SNP",
    "BETA",
    "SE",
    "N"
  )
)

gwas <- gwas[
  !is.na(SNP) &
  !is.na(BETA) &
  !is.na(SE) &
  SE > 0
]

# Remove duplicated rsIDs
gwas <- unique(
  gwas,
  by = "SNP"
)

message(
  "GWAS variants retained: ",
  format(nrow(gwas), big.mark = ",")
)

# ============================================================
# 4. Load required eQTL columns
# ============================================================

message("Loading required astrocyte eQTL columns...")

eqtl_all <- fread(
  eqtl_file,
  select = c(
    "feature",
    "variant_id",
    "fixed_beta",
    "fixed_sd"
  )
)

eqtl_targets <- eqtl_all[
  feature %in% target_genes
]

rm(eqtl_all)
invisible(gc())

eqtl_targets <- eqtl_targets[
  !is.na(feature) &
  !is.na(variant_id) &
  !is.na(fixed_beta) &
  !is.na(fixed_sd) &
  fixed_sd > 0
]

eqtl_targets <- unique(
  eqtl_targets,
  by = c(
    "feature",
    "variant_id"
  )
)

message(
  "Target eQTL rows retained: ",
  format(nrow(eqtl_targets), big.mark = ",")
)

# ============================================================
# 5. Run colocalization for every target gene
# ============================================================

coloc_results_list <- vector(
  mode = "list",
  length = length(target_genes)
)

names(coloc_results_list) <- target_genes

message(
  "Starting colocalization for ",
  length(target_genes),
  " genes..."
)

for (gene in target_genes) {

  message("Processing: ", gene)

  eqtl_gene <- eqtl_targets[
    feature == gene
  ]

  if (nrow(eqtl_gene) < 2) {
    message(
      "Skipped ",
      gene,
      ": fewer than two valid eQTL variants."
    )
    next
  }

  merged_data <- merge(
    eqtl_gene,
    gwas,
    by.x = "variant_id",
    by.y = "SNP",
    all = FALSE
  )

  merged_data <- merged_data %>%
    filter(
      !is.na(variant_id),
      !is.na(fixed_beta),
      !is.na(fixed_sd),
      fixed_sd > 0,
      !is.na(BETA),
      !is.na(SE),
      SE > 0
    ) %>%
    distinct(
      variant_id,
      .keep_all = TRUE
    )

  if (nrow(merged_data) < 2) {
    message(
      "Skipped ",
      gene,
      ": fewer than two overlapping valid variants."
    )
    next
  }

  dataset1 <- list(
    beta = merged_data$fixed_beta,
    varbeta = merged_data$fixed_sd^2,
    type = "quant",
    N = singlebrain_n,
    snp = merged_data$variant_id
  )

  dataset2 <- list(
    beta = merged_data$BETA,
    varbeta = merged_data$SE^2,
    type = "cc",
    s = scz_case_fraction,
    snp = merged_data$variant_id
  )

  # Use N only when valid sample-size information is available
  if (
    "N" %in% names(merged_data) &&
    any(is.finite(merged_data$N))
  ) {
    dataset2$N <- median(
      merged_data$N[
        is.finite(merged_data$N)
      ],
      na.rm = TRUE
    )
  }

  res <- tryCatch(
    {
      capture.output(
        coloc_object <- coloc.abf(
          dataset1 = dataset1,
          dataset2 = dataset2
        )
      )

      coloc_object
    },
    error = function(e) {
      message(
        "Coloc error for ",
        gene,
        ": ",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(res)) {
    next
  }

  summary_values <- res$summary

  coloc_results_list[[gene]] <- data.frame(
    exposure = gene,
    nsnps_coloc = as.numeric(
      summary_values["nsnps"]
    ),
    PP0 = as.numeric(
      summary_values["PP.H0.abf"]
    ),
    PP1 = as.numeric(
      summary_values["PP.H1.abf"]
    ),
    PP2 = as.numeric(
      summary_values["PP.H2.abf"]
    ),
    PP3 = as.numeric(
      summary_values["PP.H3.abf"]
    ),
    PP4 = as.numeric(
      summary_values["PP.H4.abf"]
    ),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 6. Combine and save results
# ============================================================

valid_results <- coloc_results_list[
  !vapply(
    coloc_results_list,
    is.null,
    logical(1)
  )
]

if (length(valid_results) == 0) {
  stop("No genes completed colocalization successfully.")
}

final_coloc_df <- bind_rows(
  valid_results
)

final_master <- merge(
  targets,
  final_coloc_df,
  by = "exposure",
  all.x = TRUE
)

final_master <- final_master %>%
  mutate(
    coloc_interpretation = case_when(
      PP4 >= 0.80 ~ "Strong colocalization support",
      PP4 >= 0.50 ~ "Moderate colocalization support",
      !is.na(PP4) ~ "Weak colocalization support",
      TRUE ~ "Colocalization unavailable"
    )
  ) %>%
  arrange(
    desc(PP4),
    fdr
  )

fwrite(
  final_master,
  output_file
)

message(
  "Final colocalization file saved to: ",
  output_file
)

message(
  "Genes with PP4 >= 0.80: ",
  sum(
    final_master$PP4 >= 0.80,
    na.rm = TRUE
  )
)
```

---

# Part 8: Run the Automated Colocalization Script

Interactively:

```bash
cd /scratch/negishi/hasan128/S_brain
module load r
Rscript run_automated_coloc.R
```

Expected output:

```text
Ast_SCZ_Final_Coloc_Validated_Targets.csv
```

---

# Part 9: Submit Colocalization as a SLURM Job

Create:

```text
submit_coloc.sh
```

Add:

```bash
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -A pdrineas
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=Ast_SCZ_coloc
#SBATCH --output=/scratch/negishi/hasan128/S_brain/logs/Ast_SCZ_coloc_%j.out
#SBATCH --error=/scratch/negishi/hasan128/S_brain/logs/Ast_SCZ_coloc_%j.err

set -euo pipefail

BASE_DIR="/scratch/negishi/hasan128/S_brain"
R_SCRIPT="${BASE_DIR}/run_automated_coloc.R"

mkdir -p "${BASE_DIR}/logs"

module --force purge
module load biocontainers
module load gcc
module load r

cd "${BASE_DIR}"

echo "Job started: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

Rscript "${R_SCRIPT}"

echo "Job completed: $(date)"
```

Make it executable:

```bash
chmod +x submit_coloc.sh
```

Submit it:

```bash
sbatch submit_coloc.sh
```

---

# Colocalization Posterior Probabilities

The `coloc.abf()` function returns five posterior probabilities:

| Statistic | Interpretation                                                |
| --------- | ------------------------------------------------------------- |
| `PP0`     | Neither trait is associated in the region                     |
| `PP1`     | Only the eQTL trait is associated                             |
| `PP2`     | Only schizophrenia is associated                              |
| `PP3`     | Both traits are associated but have different causal variants |
| `PP4`     | Both traits share the same causal variant                     |

A commonly used interpretation is:

|                 PP4 | Interpretation                    |
| ------------------: | --------------------------------- |
|        `PP4 ≥ 0.80` | Strong evidence of colocalization |
| `0.50 ≤ PP4 < 0.80` | Moderate evidence                 |
|        `PP4 < 0.50` | Weak evidence                     |

The joint pattern of `PP3` and `PP4` should also be examined. A high `PP3` with a low `PP4` suggests that both signals exist in the region but are likely driven by different causal variants.

---

# Important Methodological Notes

## Exposure Effect Estimates

The exposure beta, standard error, and P-value should come from the same SingleBrain model.

This workflow currently uses:

```text
fixed_beta
fixed_sd
Random_P
```

Before final analysis, confirm whether SingleBrain provides corresponding random-effects beta and standard error columns. Mixing a fixed-effect beta and standard error with a random-effects P-value is not ideal.

Use internally consistent statistics whenever possible.

---

## Duplicate SNPs

Duplicated rsIDs are removed before colocalization:

```r
distinct(
  variant_id,
  .keep_all = TRUE
)
```

However, duplicated rsIDs can represent multi-allelic or inconsistently annotated variants. For final analysis, inspect duplicated variants rather than selecting an arbitrary row without review.

---

## Minor Allele Frequency

This automated script does not assign a constant dummy MAF.

When beta and variance are supplied, `coloc.abf()` may run without MAF, depending on the dataset specification. However, real ancestry-matched allele frequencies are preferable and should be added when available.

Do not interpret results based on a fabricated uniform MAF as definitive evidence.

---

## Schizophrenia Case Fraction

The current script uses:

```r
scz_case_fraction <- 0.35
```

This value must be replaced with the exact case proportion:

```text
number of cases / total sample size
```

For example:

```r
scz_case_fraction <- cases / (cases + controls)
```

The exact value should come from the documentation for the schizophrenia GWAS dataset being analyzed.

---

## Genome Build

The SingleBrain and schizophrenia datasets may use different genomic coordinate builds.

The primary merge is therefore performed using rsIDs:

```text
SingleBrain: variant_id
GWAS: SNP
```

Using rsIDs avoids direct coordinate mismatches, but it does not fully resolve:

* outdated rsIDs
* merged dbSNP identifiers
* multi-allelic variants
* strand inconsistencies
* build-specific variant representation

For final publication-level analysis, variant normalization and allele alignment should be verified carefully.

---

## Colocalization Is Not Proven Causality

A significant MR result alone does not prove that gene expression and schizophrenia share the same causal variant.

A high `PP4` strengthens the evidence that:

1. an eQTL association exists
2. a schizophrenia association exists
3. both signals may share the same underlying causal variant

Additional sensitivity analyses are still required.

---

# Recommended Final Output Columns

The final file should contain:

```text
exposure
hgnc_symbol
description
method
nsnp
b
se
pval
fdr
nsnps_coloc
PP0
PP1
PP2
PP3
PP4
coloc_interpretation
```

Final output:

```text
Ast_SCZ_Final_Coloc_Validated_Targets.csv
```

---

# Pipeline Summary

```text
SingleBrain astrocyte eQTL data
              │
              ▼
Significant eQTL filtering
              │
              ▼
LD clumping
              │
              ▼
Exposure dataset construction
              │
              ▼
Schizophrenia GWAS matching
              │
              ▼
Allele harmonization
              │
              ▼
Mendelian Randomization
              │
              ▼
FDR correction
              │
              ▼
Ensembl-to-HGNC annotation
              │
              ▼
Gene-level colocalization
              │
              ▼
Final prioritized targets
```
