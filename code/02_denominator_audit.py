"""Audit of the implicit denominator of the repository TMB (Section 2.4, 3.1).

For each sample, implicit denominator = nonsynonymous count / repository TMB.
A fixed 30 Mb denominator implies that repository TMB x 30 is an integer.
"""
import numpy as np
import pandas as pd

from common import ANALYTIC, RESULTS, STUDIES, HIGH_TMB, wilson

d = pd.read_csv(ANALYTIC)
d["cohort"] = d["studyId"].map(STUDIES)
# 18 TCGA samples had no repository TMB; their TMB was filled with 0 in the
# analytic table (TMB_SOURCE == "calculado"), so treat them as missing here.
d.loc[d["TMB_SOURCE"] == "calculado", "TMB"] = np.nan

rows = []
for c, g in d.groupby("cohort"):
    gh = g[g["TMB"].notna() & (g["TMB"] > 0)]
    x30 = gh["TMB"] * 30
    integer = np.isclose(x30, np.round(x30), atol=1e-6).mean() * 100
    implicit = gh["N_NONSYN_MUT"] / gh["TMB"]
    k = int((g["TMB"].fillna(0) >= HIGH_TMB).sum())
    p, lo, hi = wilson(k, len(g))
    rows.append({
        "cohort": c,
        "samples": len(g),
        "samples_without_repository_TMB": int(g["TMB"].isna().sum()),
        "pct_TMBx30_integer": round(integer, 1),
        "implicit_denominator_median_Mb": round(implicit.median(), 2),
        "median_repository_TMB": round(g["TMB"].median(), 2),
        "high_TMB_repository_n": k,
        "high_TMB_repository_pct": round(p, 2),
    })
out = pd.DataFrame(rows)
print(out.to_string(index=False))
miss = d[d["TMB"].isna()]
print("\nSamples without repository TMB:", len(miss),
      "| all with zero nonsynonymous mutations:", bool((miss["N_NONSYN_MUT"] == 0).all()))
k = int((d["TMB"].fillna(0) >= HIGH_TMB).sum())
print(f"Pooled prevalence with repository values: {100 * k / len(d):.2f}% ({k}/{len(d)})")
out.to_csv(RESULTS / "denominator_audit.csv", index=False)
