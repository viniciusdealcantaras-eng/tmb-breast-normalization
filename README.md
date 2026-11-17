# Platform normalization and a null expectation for TMB reanalysis

Reproducibility code for:

> de Alcantara Sousa LV, Cubero DIG, Del Giglio A. Platform Normalization and a Null Expectation as Preconditions for Reanalyzing Tumor Mutational Burden: A Reassessment of 3381 Breast Cancer Samples from Public Genomic Cohorts. *Int. J. Environ. Res. Public Health* (under review, manuscript ijerph-4617917).

**No patient-level data are included.** All data come from the public cBioPortal for Cancer Genomics and are downloaded by the first script. Please cite cBioPortal (Cerami et al., Cancer Discov 2012; Gao et al., Sci Signal 2013) and the original studies when using these data, and follow the terms of use of each study.

## Studies

| Name in the paper | cBioPortal study ID | Platform |
|---|---|---|
| TCGA-BRCA | `brca_tcga_pan_can_atlas_2018` | Whole exome |
| MSK-IMPACT | `breast_msk_2018` | Targeted panel (341, 410 and 468 genes) |
| MBCProject | `brca_mbcproject_2022` | Whole exome |

Data were accessed on 20 July 2026. The repository is updated over time, so a new download may differ slightly from the published numbers.

## How to run

```bash
pip install -r requirements.txt
cd code
python 01_extract_cbioportal.py   # downloads data, writes data/analytic_table.csv and data/gene_counts_pooled.csv
python run_all.py                 # scripts 02 to 09, results go to results/
```

## Scripts

| Script | What it does | Paper |
|---|---|---|
| `common.py` | Settings, assay-normalized TMB, helpers | |
| `01_extract_cbioportal.py` | Clinical data, nonsynonymous mutation counts, repository TMB, gene counts | Sections 2.2 to 2.4 |
| `02_denominator_audit.py` | Implicit denominator of the repository TMB (count / TMB, TMB x 30 integer check) | Section 3.1 |
| `03_prevalence_cohorts.py` | Prevalence of high TMB, sample-type comparison, 38 Mb exome sensitivity, footprint and cut-off sensitivity | Tables 2, 3, S5, S6 |
| `04_repeated_samples.py` | One sample per patient, 1,000 random draws, standard errors clustered by patient, cohort and subtype adjustment | Table S7 |
| `05_gene_null.py` | Fisher tests with Benjamini-Hochberg FDR, length-aware null (mean and median load ratio) | Figure 2, Table S1, Data S1 |
| `06_survival.py` | Kaplan-Meier, log-rank, Cox models within cohort, Schoenfeld test, one sample per patient | Tables S3, S4, Figure S2 |
| `07_continuous_ci.py` | Hodges-Lehmann shift with cluster-bootstrap 95% CI for continuous contrasts | Section 3.2 |
| `08_subtype.py` | TMB by receptor subtype | Table S2 |
| `09_gene_length.py` | Baseline mutation frequency versus coding-sequence length for the 21 genes in Figure 2 (canonical Ensembl transcripts, Ensembl REST API) | Figure S1 |

## Key definitions

- **Nonsynonymous variants:** missense, nonsense, nonstop, frameshift and in-frame indels, splice-site and translation-start variants.
- **Repository TMB:** `TMB_NONSYNONYMOUS` from cBioPortal. The audit shows it equals the nonsynonymous count divided by 30 Mb in all three studies.
- **Assay-normalized TMB:** exome cohorts keep the repository value. MSK-IMPACT is recomputed as nonsynonymous count / 1.06 Mb (IMPACT410). Sensitivity analyses use 0.98, 1.22 Mb and the version-specific footprint, with the assay version read from the sample identifier suffix (IM3, IM5, IM6).
- **High TMB:** 10 or more mutations/Mb, used as an operational research definition.
- **Receptor subtype:** MSK-IMPACT from ER and PR status of the sequenced sample and HER2 FISH status. TCGA-BRCA from the PAM50 classification (Luminal A and B as HR+/HER2-, HER2-enriched as HER2+, Basal-like as TNBC). Not available for MBCProject.
- **Age in MBCProject:** released in 5-year bands; the lower bound of the band is used (descriptive only).
- **MSK-IMPACT panel genes** (excluded from the gene-level analysis): `code/msk_impact_468_genes.csv`, from the FDA decision summary DEN170058.

## Notes on reproducibility

- The extraction was checked against a new download of the three studies (November 2026): sample and patient counts, nonsynonymous mutation counts per sample, receptor subtype, sample type, survival data and gene-level counts were identical to those used in the paper.

- The random draws (script 04) and the bootstrap (script 07) use a fixed seed, and the manuscript reports the values obtained with it. Another seed may change percentiles in the second decimal.
- The gene-length correlation (script 09) depends on the Ensembl release used for coding-sequence lengths.

## Software

Python 3.12 or later, pandas, NumPy, SciPy, statsmodels, lifelines, matplotlib and requests (versions in `requirements.txt`).

## License

Code under the MIT License. The data are not redistributed and remain under the terms of cBioPortal and of each original study.

## Contact

Luiz Vinicius de Alcantara Sousa, Centro Universitário FMABC, Santo André, Brazil (luiz.alcantara@fmabc.net).
