"""Overall survival by TMB group, one sample per patient
(Sections 2.8 and 3.4, Tables S3 and S4, Figure S2).

Cox models are fitted within each cohort because the survival clock differs
(diagnosis in TCGA-BRCA, sequencing date in MSK-IMPACT). A model stratified by
cohort is a secondary summary. MBCProject has no survival data in the public
files. Patients with missing or zero follow-up are excluded (no imputation).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test

from common import RESULTS, load_analytic, one_per_patient

d = one_per_patient(load_analytic())
d = d[d.cohort.isin(["TCGA-BRCA", "MSK-IMPACT"])].copy()
d = d.rename(columns={"OS_MONTHS_STD": "time", "OS_EVENT_STD": "event"})
n0 = d.groupby("cohort").size()
d = d[d.time.notna() & (d.time > 0) & d.event.notna()].copy()
print("Excluded (missing or zero follow-up):", (n0 - d.groupby("cohort").size()).to_dict())

d["age_sd"] = (d.AGE_STD - d.AGE_STD.mean()) / d.AGE_STD.std()
tcga = d[d.cohort == "TCGA-BRCA"]
msk = d[d.cohort == "MSK-IMPACT"]


def fit(df, covs, strata=None, label=""):
    cph = CoxPHFitter()
    cols = ["time", "event"] + covs + ([strata] if strata else [])
    cph.fit(df[cols].dropna(), "time", "event", strata=[strata] if strata else None)
    ph = proportional_hazard_test(cph, df[cols].dropna(), time_transform="rank").summary
    used = df[cols].dropna()
    out = []
    for cov in covs:
        s = cph.summary.loc[cov]
        out.append({"model": label, "covariate": cov, "patients": len(used),
                    "deaths": int(used.event.sum()), "high_TMB_patients": int(used.high.sum()),
                    "deaths_high_TMB": int(used[used.high == 1].event.sum()),
                    "HR": round(s["exp(coef)"], 2),
                    "CI": f"{s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}",
                    "p": round(s["p"], 3), "schoenfeld_p": round(ph.loc[cov, "p"], 3)})
    return out


rows = (fit(tcga, ["high"], label="TCGA-BRCA, crude")
        + fit(tcga, ["high", "age_sd"], label="TCGA-BRCA, adjusted for age (per SD)")
        + fit(msk, ["high"], label="MSK-IMPACT, crude")
        + fit(msk, ["high", "metastatic"], label="MSK-IMPACT, adjusted for sample type")
        + fit(d, ["high"], strata="cohort", label="Stratified by cohort (secondary)"))
s4 = pd.DataFrame(rows)
print(s4.to_string(index=False))
s4.to_csv(RESULTS / "tableS4_cox.csv", index=False)

print("\nKaplan-Meier, log-rank and median follow-up (reverse Kaplan-Meier)")
s3 = []
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, (c, g) in zip(axes, [("TCGA-BRCA", tcga), ("MSK-IMPACT", msk)]):
    lr = logrank_test(g[g.high == 1].time, g[g.high == 0].time,
                      g[g.high == 1].event, g[g.high == 0].event)
    rkm = KaplanMeierFitter().fit(g.time, 1 - g.event)
    for h, lab in [(0, "TMB < 10"), (1, "TMB >= 10")]:
        gg = g[g.high == h]
        km = KaplanMeierFitter().fit(gg.time, gg.event, label=f"{lab} (n = {len(gg)})")
        km.plot_survival_function(ax=ax, ci_show=False, color="black", linestyle="-" if h == 0 else "--")
        fu = KaplanMeierFitter().fit(gg.time, 1 - gg.event).median_survival_time_
        s3.append({"cohort": c, "group": lab, "patients": len(gg), "deaths": int(gg.event.sum()),
                   "median_follow_up_months": round(fu, 1),
                   "median_OS_months": round(km.median_survival_time_, 1),
                   "OS_60_months_pct": round(100 * float(km.predict(60)), 1)})
    s3[-1]["logrank_p"] = round(lr.p_value, 3)
    print(f"{c}: median follow-up, whole cohort {rkm.median_survival_time_:.1f} months")
    ax.set_title(c)
    ax.set_xlabel("Months")
axes[0].set_ylabel("Overall survival")
fig.tight_layout()
fig.savefig(RESULTS / "figureS2_km.png", dpi=300)
s3 = pd.DataFrame(s3)
print(s3.to_string(index=False))
s3.to_csv(RESULTS / "tableS3_km.csv", index=False)
