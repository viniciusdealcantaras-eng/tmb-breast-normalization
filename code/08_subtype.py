"""TMB by receptor subtype (Table S2).

Subtype: PAM50 in TCGA-BRCA, receptor status (ER/PR of the sequenced sample
and HER2 FISH) in MSK-IMPACT, unavailable in MBCProject (see 01_extract).
Kruskal-Wallis across subtypes, Mann-Whitney post hoc with Bonferroni
adjustment, and chi-square for the proportion of high TMB.
"""
from itertools import combinations

import pandas as pd
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu

from common import RESULTS, load_analytic

SUBTYPES = ["HR+/HER2-", "HER2+", "TNBC"]
d = load_analytic()
d = d[d.SUBTYPE_STD.isin(SUBTYPES)]
print("Samples with subtype:", len(d), d.groupby("cohort").size().to_dict())

rows = []
for scope, g in [("Pooled", d), ("TCGA-BRCA", d[d.cohort == "TCGA-BRCA"]),
                 ("MSK-IMPACT", d[d.cohort == "MSK-IMPACT"])]:
    groups = [g[g.SUBTYPE_STD == s].tmb for s in SUBTYPES]
    h, p = kruskal(*groups)
    tab = pd.crosstab(g.SUBTYPE_STD, g.high)
    chi, pchi, _, _ = chi2_contingency(tab, correction=False)
    print(f"\n{scope}: Kruskal-Wallis H = {h:.2f}, p = {p:.3g}; "
          f"chi-square (high TMB) = {chi:.2f}, p = {pchi:.2f}")
    for s, x in zip(SUBTYPES, groups):
        hi = int(g[g.SUBTYPE_STD == s].high.sum())
        rows.append({"scope": scope, "subtype": s, "n": len(x), "median": round(x.median(), 2),
                     "high_n": hi, "high_pct": round(100 * hi / len(x), 1),
                     "kruskal_H": round(h, 2), "kruskal_p": p, "chi2_p": round(pchi, 2)})
    for a, b in combinations(SUBTYPES, 2):
        pp = mannwhitneyu(g[g.SUBTYPE_STD == a].tmb, g[g.SUBTYPE_STD == b].tmb).pvalue
        print(f"  {a} vs {b}: Bonferroni p = {min(1, 3 * pp):.3g}")

s2 = pd.DataFrame(rows)
print()
print(s2.drop(columns=["kruskal_p"]).to_string(index=False))
s2.to_csv(RESULTS / "tableS2_subtype.csv", index=False)
