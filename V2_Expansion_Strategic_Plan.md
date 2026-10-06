# Socio_psy_AD — V2 → V3 Strategic Expansion Plan

**Project:** Social Connection, Cognitive Function and Alzheimer's Disease: A Public-Data Genetic and Multi-Omics Study
**Baseline:** V2.0.0 (dataset registry + implementation guide, checked 2026-10-06)
**This document:** how to add (1) wider dataset scope, (2) meta-analysis, (3) phenotype analysis, (4) AI/ML/DL, as one framework.
**Status flags used throughout**
- **V** = already in your V2 registry (URL/accession audited by you).
- **C** = candidate from my knowledge of the literature. **I did not re-verify accession, URL, N, build or licence in this session.** Every C row must pass the same metadata gate as V2 before use.

---

## 0. The one-paragraph design

V2 asks whether *different social constructs* have *distinguishable* genetic relationships with cognition and clinical AD, then asks what regulatory layers add. Its weakness is structural: it has one file per construct, so "distinguishable" cannot be separated from "this particular file/definition/sample". V3 fixes that by building a **replicated-operationalisation design**:

1. **Datasets** supply several independent-as-possible operationalisations per construct and a ring of confounder/comparator traits (education, SES, depression, BMI, sleep, hearing…).
2. **Phenotype analysis** converts the registry into a formal *design matrix* (construct dimension × definition × ascertainment × covariate adjustment × sex/APOE stratum).
3. **Meta-analysis** pools estimates *within* a construct (what is generalisable) and *tests the design-matrix moderators* (what differs). Heterogeneity becomes a result, not a nuisance.
4. **ML/DL** is used only where classical statistics are the wrong tool: variant-level pattern discovery, locus-to-gene evidence integration with abstention (your Aim 2 benchmark), and sequence-based regulatory scoring. It generates and ranks hypotheses; the meta/statistical layer tests them on held-out data.

```
 L0 Registry & harmonisation ─► L1 Phenotype ontology (design matrix X)
          │                              │
          ▼                              ▼
 L2 Estimation: h², r_g, local r_g, S & V (sampling covariance)
          │                              │
          ▼                              ▼
 L3 Meta-analysis (uses X as moderators, V as covariance) ◄──── held-out test ────┐
          │                                                                       │
          ▼                                                                       │
 L4 ML/DL (features from L1–L3; hypotheses out) ──────────────────────────────────┘
          │
          ▼
 L5 Regional multi-omics (coloc, sQTL, pQTL) ─► L6 Drug-target annotation (bounded)
```

---

## 1. Dataset expansion

### 1.1 Selection rubric (applied to every candidate)

Score each 0–3, frozen **before** looking at any association result:

| Criterion | 3 | 0 |
|---|---|---|
| Relevance to RQ | Social construct, cognition, clinical AD, or a *named* confounder in the DAG | Tangential trait |
| Compatibility | Full genome-wide stats, signed beta/SE (or signed Z + N), documented alleles/build, EUR | Top hits only, restricted, unsigned P |
| Independence | Sample-disjoint from the files it will be compared with | Same cohort, same people |
| Power | h² z ≥ 4 expected (N, polygenicity) | Underpowered for LDSC |
| Marginal information | Answers a question no existing file answers | Duplicate of an included release |

Include if Relevance ≥ 2 **and** Compatibility = 3. Independence and Power decide the *role* (primary / sensitivity / descriptive), not inclusion.

### 1.2 The overlap problem decides the design

Almost every social trait is **UK Biobank** (Neale loneliness, Hilliard isolation, Day raw loneliness/participation, Shen). On the outcome side:

| Outcome file | UKB in it? | Consequence |
|---|---|---|
| Kunkle (IGAP Stage 1) | **No** | Clean anchor: r_g with UKB social traits is not inflated by overlap. Keep as the primary AD file. |
| Bellenguez, Jansen, Wightman, Schwartzentruber (as listed in V2) | **Yes** (UKB proxy-AD cases) | Larger N, but shared participants. Use LDSC with a **free** cross-trait intercept; report as sensitivity, never as independent replication. |
| FinnGen AD endpoints (C) | No | Only genuinely independent *registry-based* AD sample. Highest-value addition for replication. |
| Davies cognition (V) | Partly (UKB + CHARGE/COGENT etc.) | Overlap with UKB social traits; intercept handling as above. |

Rule: **a release is "independent" only if the metadata cohort table says so.** Overlap table is a deliverable (V2 §3 already requires it).

### 1.3 Tiered dataset registry

Counts: Wave 1 ≈ 22 releases, Wave 2 ≈ 20, Wave 3 ≈ 20+ (descriptive/exploratory). That is "as much as possible" while keeping the *primary* family small.

#### Tier A — Social constructs (the exposures)
| ID | Release | Construct / dimension | Role | Flag |
|---|---|---|---|---|
| A1 | Neale 2020 loneliness | Subjective appraisal (binary item, linear scale) | Primary | V |
| A2 | Day 2018 raw loneliness (archive member) | Appraisal, alt. coding | Sensitivity (same cohort) | V |
| A3 | Shen 2025 loneliness | Appraisal | Sensitivity (same cohort) | V |
| A4 | Hilliard 2026 isolation (0–2 index) | Structure (contact/household) | Primary | V |
| A5 | Shen 2025 isolation | Structure, alt. definition | Definition comparison | V |
| A6 | Neale 1031 friend/family visit frequency | Structure — single component | Component decomposition | C (field known; file to verify) |
| A7 | Neale 709 household size ("lives alone") | Structure — single component | Component decomposition | C |
| A8 | Day 2018 religious-group attendance | Participation | Primary candidate (freeze one) | V |
| A9 | Day 2018 sports-club attendance | Participation | Sensitivity | V |
| A10 | Day 2018 pub/social-club attendance | Participation | Sensitivity | V |
| A11 | Neale 2110 able to confide | Functional support (narrow) | Primary candidate | V |
| A12 | Hilliard 2026 loneliness GWAS | Appraisal, third UKB definition | Sensitivity | C (same Bristol deposit; verify) |
| A13 | Abdellaoui et al. loneliness | Appraisal, partly non-UKB cohorts | **Best independence for loneliness** — verify full public file | C |
| A14 | Day 2018 MTAG | Latent social-interaction composite | Labelled exploratory only | V (excluded from main construct) |

Why: A6/A7 let you ask *which component of the isolation index* carries any AD/cognition signal; that is more informative than one composite. A13, if open, is the only route to a loneliness estimate not entirely inside UKB.

#### Tier B — Primary outcomes
| ID | Release | Role | Flag |
|---|---|---|---|
| B1 | Kunkle 2019 clinical AD (GCST007511) | Primary AD anchor | V |
| B2 | Bellenguez 2022 (GCST90027158) | AD+proxy sensitivity; also temporal hold-out for ML (loci new vs Kunkle) | V |
| B3 | Schwartzentruber 2021 (GCST90012877) | Clinical/proxy comparison | V |
| B4 | Jansen 2019 AD (UKB proxy + IGAP) | Proxy sensitivity | C |
| B5 | Wightman 2021 (no-23andMe) | Proxy sensitivity | C |
| B6 | de Rojas 2021 (EADB) | Clinical-heavy comparison | C |
| B7 | FinnGen AD / dementia endpoints (exact endpoint names to confirm) | **Independent replication**; all-cause dementia as separate outcome | C (registration may be required) |
| B8 | Davies 2018 general cognition, public subset N=282,014 (DOI 10.7488/ds/2558) | Primary cognition | V |
| B9 | Savage 2018 intelligence (CTG portal) | Cognition, different task battery | C |
| B10 | Trampush/COGENT-CHARGE cognition | Cognition, reduced UKB overlap | C |

#### Tier C — Cognitive domains and secondary cognition
Davies reaction time N=330,069 (V), Neale 399 pairs-matching errors (V-lead), Neale 20018 prospective memory (V-lead), Neale 20016 fluid intelligence (C, *not* memory — label correctly). Educational attainment goes to Tier E (confounder), not here.

#### Tier D — Emotional comparators
| Release | Role | Flag |
|---|---|---|
| Howard 2019 depression (V2 existing) | Existing comparator | V |
| PGC2025 EUR no-23andMe/no-UKB (figshare 51487022, MD5 in registry) | **Preferred** — reduced UKB overlap | V |
| PGC2025 EUR no-23 (with UKB), Clinical-EUR | Definition/overlap sensitivities (same underlying release family) | V |
| Wray 2018 w/o UKB | Reduced-overlap comparison | V |
| Neuroticism: Nagel 2018 / Neale 20127 | Emotional comparison; contains a loneliness item → overlap flagged | V-lead |
| Anxiety PGC anx2026 no-23 (figshare 31389910) | Secondary emotional; do **not** use headline N | V-lead |
| Positive affect: Baselmans 2019 (no-23andMe) / Neale 4526 happiness (proxy) | Conditional | V-lead |
| Subjective wellbeing (Okbay 2016) | Optional | C |
| PTSD (Nievergelt 2024), schizophrenia (PGC3), bipolar (PGC3), ADHD (Demontis 2023) | Wave 3 psychiatric context | C |

#### Tier E — Confounders / mediating-candidate traits (the part V2 lacks)
These are what make "social → AD" interpretable. Without them, any r_g is uninterpretable because education, SES, depression and health all sit upstream.

| Trait | Why needed | Flag |
|---|---|---|
| Educational attainment (EA4, no-23andMe public release) | Strongest known confounder of cognition and social traits | C |
| Household income, Townsend deprivation (UKB-based) | SES / selection | C |
| BMI (Yengo 2018; Neale 21001_irnt) | Cardiometabolic | V / V |
| T2D (DIAMANTE-EUR, BMI-unadjusted) | Metabolic | V-lead |
| CAD (Aragam 2022, GCST90132314/5 — verify EUR file) | Vascular | V |
| Insomnia (Neale 1200; Lane 2019) | Sleep | V / V-lead |
| Smoking initiation / alcohol (GSCAN) | Behavioural | C |
| Physical activity (accelerometer-based) | Behavioural | C |
| Hearing loss | Plausible loneliness–dementia route | C |
| CRP / lipids (LDL, HDL) | Inflammation / lipid–AD link | C |

#### Tier F — AD endophenotypes and neuro-outcomes (Wave 2–3)
| Trait | Why | Flag |
|---|---|---|
| CSF Aβ42 / p-tau / amyloid-PET | Biological AD vs diagnosis | C |
| Hippocampal volume (ENIGMA), UKB BIG40 imaging phenotypes | Structure | C |
| White-matter hyperintensity, stroke | Vascular contribution | C |
| Lewy body dementia, FTD, Parkinson | Specificity: is the signal AD-specific? | C |

#### Tier G — Survival (V2 already lists)
Deelen 90th/99th percentile (GCST008598/9), Timmers parental lifespan (10283/3209), Zenin healthspan (Zenodo 1302861; P is −log10). Keep as secondary; healthspan overlaps component diseases, parental lifespan is not independent genotyped N.

#### Tier H — Molecular layers
V2 core: SingleBrain microglia (+ excitatory neuron comparison), BrainMeta v2 splice context, ADSP FunGen NG00184 protein context. Optional additions (all C): MetaBrain cortex eQTL, GTEx v8 brain tissues, Bryois-type brain snRNA-seq eQTL, eQTLGen and UKB-PPP (blood, labelled non-brain), a brain pQTL set from ROSMAP/ Banner-type cohorts. **Rule from V2 stands: dense full associations only; resources sharing donors are not replications.**

### 1.4 Explicit exclusions (and why)
| Excluded | Reason |
|---|---|
| 23andMe-containing full summaries | Not openly redistributable; use no-23andMe releases only |
| ELSA / individual-level data, individual-level PRS validation | Public aggregate scope of V2 |
| Multi-ancestry social/AD files as primary | V2 aims EUR only; add as post-hoc transferability note at most |
| Top-hits-only files | Cannot support LDSC or colocalisation |
| MTAG/composite as main social construct | Mixes constructs, defeats distinguishability |
| Restricted validated-scale affect data | Outside public-only core |
| All 171+ pairwise tests as one family | Multiplicity; see §3.6 |

### 1.5 Complexity vs comprehensiveness
- **Primary** inference stays at V2's ≤ 8 tests (4 social × 2 outcomes). Expansion adds *replications and moderators of those same 8*, not new headline tests.
- **Secondary** (social × cognitive-domain/emotional/health): own FDR family, fixed before results.
- **Descriptive** panel-wide matrix (~50–60 traits): shown as a heatmap/embedding with no per-cell significance claims beyond a stated FDR family.
- Marginal cost per added release is the metadata audit, not computation. So the hard cap is **audit capacity**: no more than ~6 new releases through the full gate per fortnight.

---

## 2. Phenotype analysis strategy

"Phenotype analysis" here means analysis of *genetic* phenotype structure with public summary data; no individual-level phenotypes.

### 2.1 Construct framework (the ontology that drives everything)

| Dimension | Meaning | Primary instrument | Alternates |
|---|---|---|---|
| Appraisal | Subjective loneliness | A1 | A2, A3, A12, A13 |
| Structure | Objective contact/household | A4 | A5, A6, A7 |
| Participation | Organised group activity | A8 (frozen choice) | A9, A10 |
| Function | Confiding support | A11 | none public → stated as narrow |
| Outcome-cognition | General cognition | B8 | B9, B10, Tier C domains |
| Outcome-clinical | Diagnosed AD | B1 | B3, B6, B7 |
| Outcome-proxy | AD/proxy | B2 | B4, B5 |
| Comparator-emotional | Depression/anxiety/neuroticism | PGC2025 noUKB | Howard, Wray |

### 2.2 Harmonised variable dictionary (what must be recorded per release)
| Field | Rule |
|---|---|
| `construct_dim` | One of the dimensions above |
| `definition_class` | item / index / component / composite |
| `direction` | Keep the author's direction in the raw file; add `harm_sign` ∈ {+1,−1} so that higher = *adverse social state* (more lonely, more isolated, less participation, less confiding) and higher = *better* cognition. Sign flips are logged, never silent |
| `effect_scale` | log-OR / linear-binary / SD-quantitative / hazard / signed-Z-only |
| `ascertainment` | clinical / proxy / registry / self-report |
| `ukb_overlap` | yes / no / partial (from cohort table) |
| `n_policy` | per-SNP N, case/control, effective-N method |
| `covariates` | age, sex, PCs, BMI-adjusted?, education-adjusted? |
| `build`, `eaf_source` | as in V2 §3 |

Critical: linear-binary (Neale) and log-OR (Kunkle) effects are **not** on the same scale; r_g is scale-free but liability-scale h² is not. Do not mix in effect-size tables.

### 2.3 Analyses (with stratifications and their logic)

| # | Analysis | Question | Method | Stratifier and why |
|---|---|---|---|---|
| P1 | Discriminant validity of social constructs | Are loneliness, isolation, participation, support genetically distinct? | r_g with CI vs 1 and vs 0; Genomic SEM 1- vs 2- vs 3-factor model comparison | Tests the V2 premise directly |
| P2 | Component decomposition of isolation | Which component drives signal? | r_g of A6, A7 vs outcomes; compare with A4 | Removes "index hides components" |
| P3 | Conditional/residual constructs | Social-specific genetic signal beyond depression/education? | Genomic SEM **GWAS-by-subtraction** (residual loneliness after depression; after EA) | Statistical decomposition only; not mediation. Requires adequate power — report model fit and h² of residual first |
| P4 | Sex-stratified | Differences by sex? | Neale male/female files (both exist for round-2 traits); r_g between sexes; sex-specific r_g with outcomes | Tests if aggregate hides opposing effects. Only where file exists per sex |
| P5 | APOE-excluded | Is AD overlap just APOE? | Repeat LDSC with chr19:44–46 Mb (GRCh37) removed on both sides | V2 §9 rule |
| P6 | Ascertainment | Clinical vs proxy AD | B1 vs B2/B4/B5 vs B7 | Proxy files carry parental-history bias and UKB overlap |
| P7 | Functional partitioning | Where in the genome / which cell types? | Stratified LDSC with brain cell-type and regulatory annotations (microglia, neurons, astrocytes, oligodendrocytes) | Gives molecular interpretation *before* coloc and informs L5 target selection |
| P8 | Local genetic correlation | Sharing in specific regions despite null global r_g | LAVA (or equivalent) on prespecified blocks, multiple-testing within family | Selects regions for L5 using a rule, not drug availability |
| P9 | Genetic-structure embedding | What are the latent axes of the whole panel? | PCA / bNMF on trait × variant Z matrix (see §4) | Descriptive; feeds ML |
| P10 | Comparison with published phenotypic claims | Does genetics agree with observational literature? | Compare r_g sign/magnitude to Luchetti-type hazard ratios, qualitatively | Triangulation, not validation |
| P11 | Survival and healthspan context | Is any social signal shared with longevity? | r_g with Tier G | Labelled secondary |

**Not stratified (stated limit):** age strata and disease stage — public summaries seldom provide them. Say so rather than imply it.

### 2.4 How phenotype analysis connects to the other components
- The design-matrix columns in §2.2 are the **moderators** in the meta-analysis (§3).
- P7/P8 outputs are **features and region lists** for ML and for the multi-omics pilot.
- P1 determines whether pooling operationalisations is *justified*; if r_g between two "loneliness" files is well below 1, they are not pooled, they are reported separately.

---

## 3. Meta-analysis integration

### 3.1 What it answers that single datasets cannot
1. **Generalisability:** is the social-trait ↔ outcome r_g stable across operationalisations of the same construct, or specific to one file?
2. **Moderation:** do definition class, ascertainment, UKB overlap, sex, or APOE inclusion systematically shift r_g? (This is the formal test of V2's "distinguishable".)
3. **Precision:** pooled r_g with a prediction interval, instead of eight fragile point estimates.
4. **Sample-independent replication:** FinnGen/non-UKB releases tested against UKB-containing estimates.

### 3.2 Three meta-analytic levels (non-duplicative)

| Level | Unit pooled | Why it is legitimate | Model |
|---|---|---|---|
| **M1 Construct-level** (primary) | r_g (Fisher-z) for construct c × outcome o across releases j | Same estimand, different measurements | Multivariate random effects (`metafor::rma.mv`): `z ~ N(Xβ, V + τ²I)` where **V is the sampling covariance of the r_g estimates** (see 3.3), X holds moderators |
| **M2 Variant/locus-level** | Per-SNP beta/SE | **Only for sample-disjoint pairs** (e.g., Kunkle + FinnGen AD). Overlapping releases (Bellenguez, Jansen, Wightman) are *not* pooled | Inverse-variance fixed/random effects; overlap-corrected where partial overlap is documented |
| **M3 Evidence-level** | MR estimates / QTL-mediated effects across contexts | Same causal-hypothesis estimand per gene | IVW random effects, with Q/I² on instruments; contexts (cell types) pooled with donor-overlap covariance |

**Do not** meta-analyse colocalisation PP.H4 values (probabilities are not additive effect sizes); combine them in the evidence matrix and the ML layer instead.

### 3.3 The statistical detail that matters
- **Effect size:** r_g via LDSC or HDL, Fisher-z transformed (SE via delta method, with jackknife SE). h² is not pooled across traits.
- **Non-independence:** estimates sharing a trait (Neale-loneliness×Kunkle and Hilliard-isolation×Kunkle share Kunkle) or sharing participants are correlated. Take **V from the block-jackknife pseudo-values** (Genomic SEM's `S` and `V` output gives exactly this). Ignoring it overstates precision.
- **Sample overlap:** LDSC cross-trait intercept left free; report intercept; run a sensitivity with Kunkle-only (no UKB) outcome.
- **Small k:** with k < 5 per construct, τ² and I² are unstable. Use Hartung–Knapp adjusted CIs, report prediction intervals, and optionally a Bayesian random effect with a weakly informative half-normal prior on τ. State when k is too small to estimate heterogeneity and report the range instead.
- **Moderator model:** `z_ij = β0 + β1·definition_class + β2·ascertainment + β3·ukb_overlap + β4·sex + ε`. With few estimates, fit moderators one at a time and report *all* attempted.
- **Heterogeneity reading:** Cochran Q with full V, I² (Higgins–Thompson adapted for correlated data), τ², 95% prediction interval, leave-one-release-out influence.
  - Low heterogeneity → construct is robust to operationalisation (supports pooling).
  - High heterogeneity driven by `definition_class` → constructs are not interchangeable (a substantive finding for the thesis).
  - High heterogeneity driven by `ukb_overlap` → overlap bias; down-weight UKB-contained outcomes.
- **Publication bias:** not applicable to a prespecified, fixed registry; state this instead of running funnel plots on k = 3.

### 3.4 Two-way tests that must use covariance
Comparing r_g(loneliness, AD) with r_g(isolation, AD): use the **contrast within the joint V**, not two independent SEs (V2 §9 warning). Same for APOE-in vs APOE-out.

### 3.5 MR as a labelled triangulation, not a causal claim
Bounded and prespecified: univariable MR (IVW, weighted median, MR-Egger, CAUSE/LCV for correlated pleiotropy), Steiger directionality, **reverse MR (AD → social)**, and multivariable MR adjusting for depression and education. Weak instruments for social traits are expected (few genome-wide-significant loci); report F-statistics, and treat any null as imprecise. This stays consistent with V2: genetic covariance ≠ environmental causation.

### 3.6 Multiplicity families (frozen before results)
| Family | Content | Control |
|---|---|---|
| F1 | 4 social × 2 outcomes (≤ 8) | Holm |
| F2 | Social × cognitive domains/emotional/health | FDR within family |
| F3 | Moderator tests in M1 | FDR within family |
| F4 | Local r_g blocks (P8) | FDR over blocks |
| F5 | Panel-wide matrix | Descriptive, FDR stated, no headline claims |

### 3.7 How meta-analysis layers on the rest (not duplicates it)
- It consumes L2's `S` and `V`; it does not re-estimate r_g.
- Its **moderator table is the phenotype design matrix**.
- Its heterogeneity result is an ML **label/feature** (which constructs are stable) and gates which pooled estimates go to the report.

---

## 4. Computational intelligence (ML/DL)

### 4.1 Where ML helps and where it does not

| Task | Classical stats enough? | ML/DL added value |
|---|---|---|
| Global h², r_g, hypothesis tests | **Yes** — keep classical; ML gives no valid inference here | None |
| Trait-level structure (~50–60 traits) | PCA/factor analysis | Small N; neural nets overfit. Use bNMF/PCA |
| **Variant-level** structure (≈ 1M SNP × 60 traits) | PCA only | **Yes — large N; VAE/deep embeddings are justified** |
| Locus→gene prioritisation with heterogeneous, partly missing, correlated evidence | Hand-built scores | **Yes — supervised, calibrated, with abstention** |
| Regulatory effect of non-coding variants in microglia/neurons | Annotation overlap | **Yes — sequence-to-function DL** |
| Drug repurposing hypotheses | Database lookup | Conditional (see M4) |
| Individual risk prediction | — | **Out of scope (V2 limit)** |
| PRS predictive accuracy | — | Cannot be validated without individual data; do not claim |

### 4.2 Four modules

**M-A. Pleiotropy embedding of the panel (descriptive)**
- Input: Z-score matrix across the registry, LD-pruned, aligned to harmonised direction.
- Methods: PCA, Bayesian NMF (bNMF), and a β-VAE on variant × trait Z as the DL arm.
- Output: latent "genetic axes" (e.g., a vascular/metabolic axis, an emotional axis, a cognitive axis) and each social trait's loadings.
- Role: tests whether social constructs load on *different* axes — an unsupervised check of V2's premise that complements P1.
- Validation: reproducibility across random LD-pruned subsets and across chromosomes (train on odd, test on even); compare against PCA baseline; stability of axes under leaving out each trait.

**M-B. Evidence-integration model for locus→gene (supervised, Aim 2 backbone)**
- Features per locus-gene pair: distance, eQTL/sQTL/pQTL coloc strength per layer, SMR/HEIDI, local r_g, S-LDSC cell-type annotation, PoPS-style gene features, enhancer–gene links, constraint, Open Targets locus-to-gene scores, DL regulatory delta-score (M-C), *plus missingness indicators* (missing protein ≠ negative evidence, V2 §11).
- Models: logistic regression baseline → gradient boosting (LightGBM/XGBoost) → optionally a small neural model with missing-layer masking. Output calibrated probabilities (isotonic/Platt) and **conformal abstention** (this directly delivers V2's "accuracy–coverage trade-off including abstention").
- Labels: **gold-standard gene sets independent of the loci being scored** — e.g., Mendelian/rare-variant-supported AD genes, and for cognition a distinct list. Because gold sets are small and noisy, treat as comparison sets (V2 §11 statement), report with uncertainty.
- Validation:
  - Leave-locus-out and leave-chromosome-out CV; all genes in one LD block in the same fold.
  - **Temporal hold-out:** train on loci from Kunkle/IGAP-era, test on loci new in Bellenguez.
  - Baselines: nearest gene, expression-only, L2G/PoPS-type benchmarks.
  - Ablation: add layers one at a time (expression → +splicing → +protein).
  - Negative controls: phenotype-label permutation, randomly matched control traits.
  - Metrics: AUPRC, Brier score, reliability curves, coverage–accuracy curve, wrong-gene rate.
  - Interpretability: SHAP with a stated caveat that attributions on correlated features are not causal.

**M-C. Sequence-to-function DL for credible-set variants**
- Pretrained sequence models (e.g., Enformer/Borzoi/AlphaGenome-type) scored on fine-mapped variants in microglia/neuron-relevant tracks; requires the audited hg38 lift (V2 §11).
- Use as a *feature* in M-B and as an orthogonal check on QTL colocalisation, not as stand-alone evidence.
- Validation: concordance with measured QTL direction where available; held-out chromosomes; report calibration of predicted vs observed effect sign.
- Limit: models are trained largely on non-brain-aged data; microglia state specificity is imperfect — report as a hypothesis generator.

**M-D. Target / repurposing graph (optional, go/no-go at month 24)**
- Knowledge-graph embedding (e.g., RotatE/TransE-type via a standard KG library) on Open Targets-derived edges.
- Only if M-B yields enough confident genes. Evaluate with a **time-split**: edges before a cutoff predict later clinical-trial entries; avoid leakage from the same drug–target pair.
- Output remains "target and repurposing hypotheses for validation" (V2 §12).

### 4.3 How ML outputs integrate with statistics
| Interaction | Mode |
|---|---|
| M-A axes vs P1 / M1 | **Convergent validity** — if axes agree with r_g/factor model, support; if not, report disagreement, not a winner |
| M-B probabilities vs coloc/SMR | **Different questions**: coloc asks "shared signal at this locus?"; M-B asks "which gene, given all layers?" M-B uses coloc as an input, so it cannot validate coloc |
| M-B vs meta-analysis | Meta-analytic stability labels (which constructs are robust) set which loci enter the benchmark; ML never overrides a pooled estimate |
| M-C vs QTL | Orthogonal evidence; disagreement is informative, not an error |
| ML → hypotheses → held-out | New loci/genes from ML must be tested on data (a later release, another cohort, another tissue) not used in training |

### 4.4 Guardrails
- No tuning on the same loci used for evaluation.
- Pre-register feature sets, hyperparameter grids, and the held-out partition.
- Report failures and abstentions.
- No claim of clinical prediction or individual-level treatment inference.

---

## 5. Why the integrated design beats V2

| V2 baseline | V3 |
|---|---|
| One file per construct; "distinguishable" confounded with file choice | Multiple operationalisations per construct; difference tested in a covariance-aware meta-regression |
| Overlap handled as a caveat | Overlap is a modelled moderator and a design variable (Kunkle clean anchor, FinnGen independent replication) |
| Confounders appear only in text | Education/SES/depression/health traits included; conditional constructs via GWAS-by-subtraction |
| Local selection rule vague | P7/P8 give prespecified region lists; ML cannot choose regions post hoc |
| Aim 2 evaluation relies on heuristics | Calibrated, abstaining ML benchmark with temporal hold-out and layer ablation |
| Panel risk = 171 pairs | Fixed families F1–F5, headline limited to 8 |
| No quantitative heterogeneity | Q, τ², I², prediction intervals as results |

The key argument for the thesis: the novel contribution (V2 §1) is "a reproducible comparison of distinct social constructs coupled to an evaluated approach to molecular interpretation". M1/P1 supports the first half rigorously; M-B supports the second half with a measurable benchmark.

---

## 6. Workflow and gates (36 months, inside V2's schedule)

| Period | Work | Gate |
|---|---|---|
| Wk 1–6 | V2 audit of 5 existing files; Wave-1 metadata gate (Tier A/B/D + EA) | Overlap table, N policy, harmonised sign table frozen |
| Mo 2–4 | Full QC; h²; P1, P2, P5, P6; M1 for F1 | Each primary trait h² z ≥ 4 or reported as underpowered |
| Mo 4–6 | Wave-2 additions (Tier E/F); P3, P4, P7, P8; **freeze F1–F5 and ML feature spec** | Supervisor/novelty review (V2 §13) |
| Mo 7–12 | Panel embedding M-A; MR triangulation; manuscript 1 | Heterogeneity interpretable |
| Mo 13–24 | Regional multi-omics pilot → full; M-B/M-C benchmark; manuscript 2 | Benchmark beats nearest-gene/expression-only baselines, or limits reported |
| Mo 25–36 | M-D (conditional), drug annotation, reproducibility, thesis | Go/no-go on M-D by month 24 |

**Stop rules:** drop any release failing the compatibility gate (do not substitute after seeing results); drop a module if its pre-specified baseline cannot be beaten, and report that.

---

## 7. Risks and honest limits
- **Power:** social-trait h² is low; r_g may be imprecise. Report intervals; null with wide CI ≠ no overlap.
- **UKB dominance:** most social data share UKB; non-UKB loneliness (A13) and FinnGen are priorities for that reason.
- **Small k in meta-analysis:** heterogeneity estimates unstable; reported with the caveat.
- **Gold-standard labels for ML are imperfect and few;** performance estimates carry wide uncertainty.
- **Residual constructs (GWAS-by-subtraction)** can be underpowered or ill-fitting; fit diagnostics first.
- **C-flagged datasets are unverified** in this session; accession, build, N, licence must be confirmed through the V2 metadata gate.
- **Genetic covariance is not causation, and cross-sectional cognition is not decline** (V2 stance retained).
- Novelty must still be re-checked against Almeida e Sousa, Hilliard & Zhao, and Shen et al. before the title is committed.

## 8. Immediate next actions
1. Add the harmonised-sign and overlap columns (§2.2) to the registry JSON schema.
2. Run the metadata gate on Wave-1 C-flagged rows: A13, B7 (FinnGen), B9, B10, EA4, plus Tier E essentials.
3. Freeze F1–F5 and the M1 model specification (moderators, V source, small-k policy).
4. Decide the single participation activity (A8 vs alternatives) before seeing r_g.
5. Specify the M-B label set and the temporal hold-out partition before any benchmark run.
