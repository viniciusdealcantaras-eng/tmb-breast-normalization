"""Patient-level and clustered sensitivity analyses (Section 2.8, Table S7).

Two hundred patients contributed more than one sample. The primary-versus-
metastatic comparison is repeated (i) with one sample per patient (lowest
sample identifier), (ii) over 1,000 random draws of one sample per patient and
(iii) with logistic models on all samples and standard errors clustered by
patient. Woolf intervals are used for the 2x2 tables, cluster-robust Wald
intervals for the logistic models.
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from common import RESULTS, SEED, load_analytic, one_per_patient, wilson, woolf_or

N_DRAWS = 1000


def or_2x2(d):
    a = d[d.metastatic == 1]["high"]
    c = d[d.metastatic == 0]["high"]
    return woolf_or(a.sum(), (1 - a).sum(), c.sum(), (1 - c).sum())


def logit_or(formula, d, cluster=True):
    kw = {"cov_type": "cluster", "cov_kwds": {"groups": pd.factorize(d["patientId"])[0]}} if cluster else {}
    r = smf.logit(formula, d).fit(disp=0, **kw)
    lo, hi = np.exp(r.conf_int().loc["metastatic"])
    return np.exp(r.params["metastatic"]), lo, hi


d = load_analytic()
msk = d[d.panel == 1]

multi = d.groupby("patientId").size()
multi = multi[multi > 1]
print("Patients with >1 sample:", len(multi), "| samples:", int(multi.sum()))
print(d[d.patientId.isin(multi.index)].groupby("cohort").size().to_dict())

rows = []


def add(analysis, scope, est):
    rows.append({"analysis": analysis, "scope": scope, "OR": round(est[0], 2),
                 "CI_or_percentiles": f"{est[1]:.2f}-{est[2]:.2f}"})


add("All samples, Woolf", "Pooled", or_2x2(d))
add("All samples, Woolf", "MSK-IMPACT", or_2x2(msk))
add("All samples, SE clustered by patient", "Pooled", logit_or("high ~ metastatic", d))
add("All samples, SE clustered by patient", "MSK-IMPACT", logit_or("high ~ metastatic", msk))
add("Pooled adjusted for cohort, clustered", "Pooled", logit_or("high ~ metastatic + C(cohort)", d))
add("MSK-IMPACT adjusted for receptor subtype, clustered", "MSK-IMPACT",
    logit_or("high ~ metastatic + C(SUBTYPE_STD)", msk))

opp = one_per_patient(d)
add("One sample per patient (lowest identifier)", "Pooled", or_2x2(opp))
add("One sample per patient (lowest identifier)", "MSK-IMPACT", or_2x2(opp[opp.panel == 1]))

rng = np.random.default_rng(SEED)
draws = {"Pooled": [], "MSK-IMPACT": []}
for _ in range(N_DRAWS):
    s = d.sample(frac=1, random_state=rng.integers(2**32 - 1)).groupby("patientId").head(1)
    draws["Pooled"].append(or_2x2(s)[0])
    draws["MSK-IMPACT"].append(or_2x2(s[s.panel == 1])[0])
for k, v in draws.items():
    v = np.array(v)
    add(f"{N_DRAWS} random draws of one sample per patient (median, 2.5th-97.5th pct)", k,
        (np.median(v), *np.percentile(v, [2.5, 97.5])))

s7 = pd.DataFrame(rows)
print(s7.to_string(index=False))

print("\nPrevalence of high TMB, one sample per patient")
prev = []
for c in ["TCGA-BRCA", "MBCProject", "MSK-IMPACT", "Exome cohorts", None]:
    g = opp if c is None else (opp[opp.panel == 0] if c == "Exome cohorts" else opp[opp.cohort == c])
    p, lo, hi = wilson(int(g.high.sum()), len(g))
    prev.append({"cohort": c or "Pooled", "n": len(g), "high": int(g.high.sum()),
                 "pct": round(p, 2), "ci": f"{lo:.2f}-{hi:.2f}"})
prev = pd.DataFrame(prev)
print(prev.to_string(index=False))

s7.to_csv(RESULTS / "tableS7_repeated_samples.csv", index=False)
prev.to_csv(RESULTS / "tableS7_prevalence_one_per_patient.csv", index=False)
