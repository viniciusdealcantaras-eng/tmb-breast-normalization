"""Prevalence of high TMB, sample-type comparisons and footprint/cut-off
sensitivity analyses (Sections 3.1-3.2, Tables 2, 3, S5, S6)."""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact

from common import (RESULTS, HIGH_TMB, load_analytic, version_footprint,
                    wilson, woolf_or, one_per_patient, fmt)


def prevalence_table(d):
    groups = [(c, d[d["cohort"] == c]) for c in ["TCGA-BRCA", "MBCProject", "MSK-IMPACT"]]
    groups += [("Exome cohorts", d[d["panel"] == 0]), ("Pooled (descriptive)", d)]
    rows = []
    for name, g in groups:
        k = int(g["high"].sum())
        p, lo, hi = wilson(k, len(g))
        q1, q3 = g["tmb"].quantile([0.25, 0.75])
        rows.append({"group": name, "n": len(g), "median": round(g["tmb"].median(), 2),
                     "IQR": f"{q1:.2f}-{q3:.2f}", "max": round(g["tmb"].max(), 1),
                     "high_n": k, "pct": round(p, 2), "ci": f"{lo:.2f}-{hi:.2f}"})
    return pd.DataFrame(rows)


def clustered_or(d):
    """Logistic OR of high TMB (metastatic vs primary), SE clustered by patient."""
    r = smf.logit("high ~ metastatic", d).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": pd.factorize(d["patientId"])[0]})
    lo, hi = np.exp(r.conf_int().loc["metastatic"])
    return round(np.exp(r.params["metastatic"]), 2), f"{lo:.2f}-{hi:.2f}"


def sample_type(d):
    a = d[d["metastatic"] == 1]["high"]
    c = d[d["metastatic"] == 0]["high"]
    o = woolf_or(a.sum(), (1 - a).sum(), c.sum(), (1 - c).sum())
    p = fisher_exact([[a.sum(), (1 - a).sum()], [c.sum(), (1 - c).sum()]])[1]
    return {
        "met_n": len(a), "met_high": int(a.sum()), "met_pct": round(100 * a.mean(), 2),
        "prim_n": len(c), "prim_high": int(c.sum()), "prim_pct": round(100 * c.mean(), 2),
        "met_median": round(d[d["metastatic"] == 1]["tmb"].median(), 2),
        "prim_median": round(d[d["metastatic"] == 0]["tmb"].median(), 2),
        "OR": round(o[0], 2), "OR_ci": f"{o[1]:.2f}-{o[2]:.2f}", "fisher_p": p,
    }


d = load_analytic()

print("== Table 2. Assay-normalized TMB and prevalence of high TMB")
t2 = prevalence_table(d)
print(t2.to_string(index=False))
t2.to_csv(RESULTS / "table2_prevalence.csv", index=False)

print("\n== Repository values for MSK-IMPACT (30 Mb)")
m = d[d["panel"] == 1]
print(f"median {m['tmb_repo'].median():.2f}; high {int((m['tmb_repo'] >= HIGH_TMB).sum())}/{len(m)}"
      f" ({100 * (m['tmb_repo'] >= HIGH_TMB).mean():.2f}%)")
print(f"median at 1.06 Mb {m['tmb'].median():.2f}; high {int(m['high'].sum())}")

print("\n== Panel vs exome medians")
print(f"panel {d[d.panel == 1]['tmb'].median():.2f} | exome cohorts {d[d.panel == 0]['tmb'].median():.2f}"
      f" | TCGA-BRCA {d[d.cohort == 'TCGA-BRCA']['tmb'].median():.2f}")

print("\n== Table 3. Sample type, pooled and within MSK-IMPACT")
t3 = pd.DataFrame([{"analysis": "Pooled", **sample_type(d)},
                   {"analysis": "MSK-IMPACT only", **sample_type(m)}])
print(t3.to_string(index=False))
t3.to_csv(RESULTS / "table3_sample_type.csv", index=False)

print("\n== Exome cohorts standardized to 38 Mb (MSK-IMPACT at 1.06 Mb)")
d38 = d.copy()
d38.loc[d38.panel == 0, "tmb"] = d38.loc[d38.panel == 0, "N_NONSYN_MUT"] / 38
print(f"pooled median {d38['tmb'].median():.2f}; pooled prevalence {100 * (d38['tmb'] >= HIGH_TMB).mean():.2f}%")

print("\n== Table S5. MSK-IMPACT footprint")
rows = []
base = d[d.panel == 1][["sampleId", "high"]].rename(columns={"high": "high_base"})
for label, den in [("0.98 Mb, all samples", 0.98), ("1.06 Mb, all samples (primary)", 1.06),
                   ("1.22 Mb, all samples", 1.22), ("Version-specific", version_footprint)]:
    x = load_analytic(den)
    mx = x[x.panel == 1].merge(base, on="sampleId")
    up = int(((mx.high == 1) & (mx.high_base == 0)).sum())
    down = int(((mx.high == 0) & (mx.high_base == 1)).sum())
    p, lo, hi = wilson(int(mx.high.sum()), len(mx))
    pp = 100 * x.high.mean()
    o_or, o_ci = clustered_or(x[x.panel == 1])
    rows.append({"denominator": label, "msk_high": int(mx.high.sum()), "msk_pct": round(p, 2),
                 "msk_ci": f"{lo:.2f}-{hi:.2f}", "pooled_pct": round(pp, 2),
                 "low_to_high": up, "high_to_low": down, "msk_OR_clustered": o_or, "msk_OR_clustered_ci": o_ci})
print("Samples by assay version:", d[d.panel == 1]["msk_version"].value_counts().to_dict())
x = load_analytic(version_footprint)
for v, g in x[x.panel == 1].groupby("msk_version"):
    p, lo, hi = wilson(int(g.high.sum()), len(g))
    rows.append({"denominator": f"  {v} ({version_footprint(g).iloc[0]} Mb), n = {len(g)}",
                 "msk_high": int(g.high.sum()), "msk_pct": round(p, 2), "msk_ci": f"{lo:.2f}-{hi:.2f}"})
s5 = pd.DataFrame(rows)
print(s5.to_string(index=False))
s5.to_csv(RESULTS / "tableS5_footprint.csv", index=False)

print("\n== Table S6. Prevalence at alternative cut-offs")
rows = []
for cut in [5, 10, 16, 20]:
    r = {"cutoff": cut}
    for c in ["TCGA-BRCA", "MBCProject", "MSK-IMPACT", "Pooled"]:
        g = d if c == "Pooled" else d[d.cohort == c]
        p, lo, hi = wilson(int((g["tmb"] >= cut).sum()), len(g))
        r[c] = f"{p:.2f} ({lo:.2f}-{hi:.2f})"
    rows.append(r)
s6 = pd.DataFrame(rows)
print(s6.to_string(index=False))
s6.to_csv(RESULTS / "tableS6_cutoffs.csv", index=False)
