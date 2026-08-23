# Rlig1 tRNA-Repair Co-Enrichment Reanalysis

An independent, reproducible R/Bioconductor reanalysis of the human RNA ligase **Rlig1 (C12orf29)** affinity-enrichment mass-spectrometry (AE-MS) dataset from **Stumpf et al. (2026, *Nucleic Acids Research*)**.

The original study processed its data in Perseus and reported that two tRNA end-processing enzymes (CLP1, ANGEL2) co-enrich with Rlig1. This reanalysis reproduces the broad interactome character in a modern pipeline and **extends it systematically**: almost the entire human tRNA end-processing and repair machinery co-enriches with active Rlig1 (**19 of 22 detected enzymes; ~26-fold over-representation; Fisher's exact p = 3.9 × 10⁻¹¹**), robustly under oxidative stress.

> **Note on claims.** The data show *co-enrichment* (physical association), not direct binding or catalytic partnership. All wording is kept to "co-enriches with" / "associates with" rather than "interacts."

---

## Biological background

- **Rlig1 (C12orf29 / RLIG1)** is a human 5′→3′, ATP-dependent RNA ligase that rejoins broken RNA. It auto-AMPylates at lysine 57 (K57); the **K57A** point mutant is catalytically dead but structurally intact (an activity control).
- Prior biochemistry (Hu et al., 2024, *PNAS*) established Rlig1 as a bacteriophage-T4-like tRNA-**repair** ligase that is **not** a tRNA-splicing enzyme and is predominantly cytoplasmic. It proposed a human tRNA-repair module (CLP1 + ANGEL2 + Rlig1) analogous to the T4 system, but stated the coordination was untested.
- This reanalysis provides systematic co-enrichment evidence spanning that pathway.

---

## Experimental design (source data — Stumpf et al., 2026)

Affinity-enrichment pull-down: recombinant Rlig1 baits capture interacting proteins from HEK293 Rlig1-KO lysates.

| Factor | Levels |
|---|---|
| Bait | AMP-Rlig1 (active), K57A-Rlig1 (inactive), Bead control |
| Stress | Unstressed, Stressed (40 µM menadione) |
| Replicates | 3 biological replicates |

**3 baits × 2 stress × 3 reps = 18 samples.** Lysates were RNase-treated before pull-down. Quantified by label-free DIA (Spectronaut, protein-group level).

---

## Input files

| File | Description |
|---|---|
| `20240710_..._WithoutRNA_Paper_Report.csv` | Spectronaut v16 protein-group report (~3,271 groups). Read with `check.names = FALSE`. Missing values = `NaN` → `NA`. Intensities in `...PG.Quantity` columns. |
| `metadata.tsv` | Tab-separated. One row per sample: `sample_id`, `group` (bait), `stress`. Read with `read.delim`. All 18 sample IDs match the report columns. |

---

## Pipeline overview

```
Spectronaut report + metadata
        │
        ▼
1. Import  ──────────────  QFeatures object (metadata matched, tripwire check)
        ▼
2. Clean & log2  ────────  remove 57 contaminants → 3,214 proteins; log2 transform
        ▼
3. Impute & normalise  ──  kNN imputation (seeded); median normalisation
        ▼
4. QC  ──────────────────  PCA (replicates cluster; baits separate from beads)
        ▼
5. Differential enrichment  limma: 6 contrasts, eBayes (robust, trend), BH-FDR
        ▼
6. Targeted pathway test  ─ Fisher's exact test on the tRNA gene set
        ▼
7. Functional & network  ─ GO (clusterProfiler); STRING clustering; IntAct query
```

### Step detail

1. **Import & setup** — `readQFeatures()` builds a `QFeatures` object with the 18 `.PG.Quantity` columns as the assay and `PG.ProteinGroups` as feature names. Metadata is matched to assay columns; `stopifnot(!anyNA(idx))` halts on any mismatch.
2. **Cleaning & transformation** — remove 57 `Cont_`-flagged contaminants (trypsin, keratins, BSA, casein); `logTransform(base = 2)` to symmetrise the distribution and make fold changes interpretable (log2 = 1 → 2-fold).
3. **Imputation & normalisation** — `impute(method = "knn", MARGIN = 1)` with `set.seed(118)`; `normalize(method = "diff.median")`. Median (not quantile) normalisation preserves the genuinely different bait/mutant/bead distributions.
4. **Quality control** — `FactoMineR::PCA(scale.unit = TRUE)` + `factoextra::fviz_screeplot`; confirms replicate clustering and that variation tracks the two designed factors.
5. **Differential enrichment** — means-model design (`~ 0 + condition`); 6 contrasts (AMP/K57A vs Beads, AMP vs K57A, per stress state); `lmFit → contrasts.fit → eBayes(robust = TRUE, trend = TRUE)`; BH-FDR. Enriched = `adj.P.Val < 0.05 & logFC > 1`.
6. **Targeted pathway test** — a priori tRNA end-processing/repair gene set; 2×2 table (tRNA vs not × enriched vs not); `fisher.test(alternative = "greater")`.
7. **Functional & network analysis** — `clusterProfiler::enrichGO` (background = all quantified proteins); STRING clustering (confidence ≥ 0.700); IntAct query for documented Rlig1 interactors.

---

## The six contrasts

| Contrast | Question |
|---|---|
| `AMP_vs_Beads_unst` / `_str` | What does active Rlig1 specifically enrich? (primary interactors) |
| `K57A_vs_Beads_unst` / `_str` | What does the inactive mutant enrich? (activity-independence) |
| `AMP_vs_K57A_unst` / `_str` | Direct activity-dependence |

---

## Key results

- **Reproduction:** GO analysis recovers an RNA- and translation-centric co-enrichment set, matching the original study.
- **Central finding:** 19 of 22 tRNA end-processing/repair enzymes co-enrich with active Rlig1 (86% Fisher p = 3.9 × 10⁻¹¹), spanning almost every pathway step (5′/3′ processing, end-healing, the RTCB ligase complex, CCA maturation).
- **Robustness:** present in both stress states; also recovered by the inactive K57A mutant (13/22), indicating the association is largely activity-independent but strengthened by activity.
- **Networks:** STRING clustering recovers the known tRNA sub-complexes; Rlig1 and ANGEL2 sit unclustered (recently-recognised repair roles, sparse database annotation). IntAct lists only ~4 curated Rlig1 interactors — confirming the association is novel.

---

## Figures

| Figure | Content | Tool |
|---|---|---|
| GO dotplot | Molecular Function enrichment (AMP vs Beads) | clusterProfiler + ggplot2 |
| tRNA pathway bar chart | Per-enzyme log2 enrichment, grouped by pathway step, coloured by significance | ggplot2 |
| Heatmap | Row-Z-scored abundance across 18 samples | pheatmap |
| STRING network | Clustered tRNA machinery + Rlig1 | string-db.org |
| UpSet plot | Enriched proteins shared across the six contrasts (future-work) | UpSetR |

---

## R packages

**Core:** `QFeatures`, `limma`, `clusterProfiler`, `org.Hs.eg.db`
**QC / stats:** `FactoMineR`, `factoextra`, `stats` (`fisher.test`)
**Networks:** `STRINGdb` (or string-db.org), IntAct database
**Wrangling / plots:** `tidyverse` (dplyr, tidyr, stringr, ggplot2), `pheatmap`, `RColorBrewer`, `UpSetR`

Install Bioconductor packages with:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("QFeatures", "limma", "clusterProfiler", "org.Hs.eg.db", "STRINGdb"))
install.packages(c("tidyverse", "FactoMineR", "factoextra", "pheatmap", "RColorBrewer", "UpSetR"))
```

---

## Reproducibility notes

- `set.seed(118)` fixes the stochastic kNN imputation, so the pipeline returns identical results on every run.
- The core tRNA finding was confirmed across imputation methods (MinProb, mixed MAR/MNAR, kNN) — it does not depend on the imputation strategy.
- Enrichment threshold: `adj.P.Val < 0.05 & logFC > 1` throughout.
- GO uses the quantified proteins as background (not the whole genome), so enrichment reflects the detectable proteome.

---

## References

Stumpf FM, Glauner M, Jansen J, Marchand V, Motorin Y, Stengel F, Marx A. Human RNA ligase 1 as a novel regulator of ribosome function and translation under oxidative stress. Nucleic Acids Research. 2026;54(11):gkag528. doi:10.1093/nar/gkag528

Orellana EA, Siegal E, Gregory RI. tRNA dysregulation and disease. Nature Reviews Genetics. 2022;23(11):651–664. doi:10.1038/s41576-022-00501-9

Hu Y, Lopez VA, Xu H, Pfister JP, Song B, Servage KA, Sakurai M, Jones BT, Mendell JT, Wang T, Wu J, Lambowitz AM, Tomchick DR, Pawłowski K, Tagliabracci VS. Biochemical and structural insights into a 5′ to 3′ RNA ligase reveal a potential role in tRNA ligation. Proceedings of the National Academy of Sciences. 2024;121(42):e2408249121. doi:10.1073/pnas.2408249121

Yuan Y, Stumpf FM, Schlor LA, Schmidt OP, Saumer P, Huber LB, Frese M, Höllmüller E, Scheffner M, Stengel F, Diederichs K, Marx A. Chemoproteomic discovery of a human RNA ligase. Nature Communications. 2023;14:842. doi:10.1038/s41467-023-36451-x



---

*Reanalysis conducted in R/Bioconductor, independently of the original Perseus workflow. Data show co-enrichment (association), not proven physical interaction or catalytic function.*
