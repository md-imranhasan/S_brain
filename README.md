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
````
