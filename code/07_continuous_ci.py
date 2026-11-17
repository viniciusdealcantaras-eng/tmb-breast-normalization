"""Confidence intervals for the continuous TMB contrasts (Section 2.7).

TMB is a count divided by a fixed denominator, so its distribution is
discrete and a bootstrap of the difference in medians is degenerate. The
location shift is therefore summarized by the Hodges-Lehmann estimator (median
of all pairwise differences) with a 95% percentile interval from 1,000
bootstrap resamples of patients (cluster bootstrap).
"""
import numpy as np
import pandas as pd

from common import RESULTS, SEED, load_analytic

N_BOOT = 1000


def hodges_lehmann(x1, x0):
    return float(np.median(np.subtract.outer(x1, x0)))


def cluster_boot(d, group, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    idx = d.groupby("patientId").indices
    pats = np.array(list(idx))
    x, g = d["tmb"].to_numpy(), d[group].to_numpy()
    est = hodges_lehmann(x[g == 1], x[g == 0])
    out = []
    for _ in range(n_boot):
        s = np.concatenate([idx[p] for p in rng.choice(pats, len(pats))])
        xs, gs = x[s], g[s]
        out.append(hodges_lehmann(xs[gs == 1], xs[gs == 0]))
    lo, hi = np.percentile(out, [2.5, 97.5])
    return est, lo, hi


d = load_analytic()
rows = []
for label, sub, grp in [
    ("Panel (MSK-IMPACT) vs exome cohorts", d, "panel"),
    ("Metastatic vs primary, pooled", d, "metastatic"),
    ("Metastatic vs primary, MSK-IMPACT", d[d.panel == 1], "metastatic"),
]:
    est, lo, hi = cluster_boot(sub.reset_index(drop=True), grp)
    m1 = sub[sub[grp] == 1]["tmb"].median()
    m0 = sub[sub[grp] == 0]["tmb"].median()
    rows.append({"contrast": label, "median_1": round(m1, 2), "median_0": round(m0, 2),
                 "difference_in_medians": round(m1 - m0, 2),
                 "hodges_lehmann": round(est, 2), "ci_95": f"{lo:.2f}-{hi:.2f}"})
out = pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv(RESULTS / "continuous_contrasts_ci.csv", index=False)
