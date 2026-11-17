"""Shared settings and helpers for the TMB reanalysis.

All paths are relative to the repository root. The derived tables in data/
are produced by 01_extract_cbioportal.py and are not distributed with this
repository (see README).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ANALYTIC = DATA / "analytic_table.csv"
GENE_COUNTS = DATA / "gene_counts_pooled.csv"
PANEL_GENES = ROOT / "code" / "msk_impact_468_genes.csv"

STUDIES = {
    "brca_tcga_pan_can_atlas_2018": "TCGA-BRCA",
    "breast_msk_2018": "MSK-IMPACT",
    "brca_mbcproject_2022": "MBCProject",
}
MSK = "breast_msk_2018"

REPO_DENOM_MB = 30.0          # denominator found in the repository audit
PANEL_DENOM_MB = 1.06         # primary footprint, IMPACT410
VERSION_FOOTPRINT_MB = {"IM3": 0.98, "IM5": 1.06, "IM6": 1.22}
HIGH_TMB = 10.0               # operational research definition (mut/Mb)
SEED = 20260720


def load_analytic(denom_panel=PANEL_DENOM_MB):
    """Analytic table with assay-normalized TMB.

    Exome cohorts keep the repository TMB (count / 30 Mb). MSK-IMPACT TMB is
    recomputed as nonsynonymous count / panel footprint. Eighteen TCGA samples
    without a repository TMB have zero nonsynonymous mutations and get 0.
    """
    d = pd.read_csv(ANALYTIC)
    d["cohort"] = d["studyId"].map(STUDIES)
    d["panel"] = (d["studyId"] == MSK).astype(int)
    d["metastatic"] = (d["SAMPLE_TYPE_STD"] == "Metastatico").astype(int)
    d["msk_version"] = d["sampleId"].str.extract(r"-(IM\d)$")[0]
    d["tmb_repo"] = d["N_NONSYN_MUT"] / REPO_DENOM_MB
    d.loc[d["panel"] == 0, "tmb_repo"] = d.loc[d["panel"] == 0, "TMB"].fillna(0)
    d["tmb"] = d["tmb_repo"]
    if callable(denom_panel):
        den = denom_panel(d)
    else:
        den = denom_panel
    msk = d["panel"] == 1
    d.loc[msk, "tmb"] = d.loc[msk, "N_NONSYN_MUT"] / (den[msk] if hasattr(den, "__len__") else den)
    d["high"] = (d["tmb"] >= HIGH_TMB).astype(int)
    return d


def version_footprint(d):
    return d["msk_version"].map(VERSION_FOOTPRINT_MB)


def wilson(k, n):
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return 100 * k / n, 100 * lo, 100 * hi


def woolf_or(a, b, c, d):
    """OR and Woolf 95% CI for a 2x2 table [[a, b], [c, d]]
    (a = exposed with outcome, b = exposed without, c/d = unexposed)."""
    or_ = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return or_, np.exp(np.log(or_) - 1.96 * se), np.exp(np.log(or_) + 1.96 * se)


def one_per_patient(d):
    """Keep the sample with the lowest identifier for each patient."""
    return d.sort_values("sampleId").groupby("patientId", as_index=False).head(1)


def fmt(x, nd=2):
    return f"{x:.{nd}f}"
