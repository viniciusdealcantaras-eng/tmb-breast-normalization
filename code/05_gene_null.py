"""Gene-level enrichment against a length-aware null expectation
(Sections 2.6 and 3.3, Figure 2, Table S1, Data S1).

Input: gene_counts_pooled.csv, the number of mutated samples per gene in the
high- and low-TMB groups defined with the repository TMB in the pooled data
(40 high, 3,341 low). For genes outside the MSK-IMPACT panel the panel samples
contribute no mutations, so these counts are also the exome-only counts
(39 high, 1,424 low): the single repository-high panel sample is not
sequenced for these genes and the panel low-TMB samples add only zeros.

Null: expected frequency in the high-TMB group = 1 - (1 - p)^r, with p the
low-TMB mutated proportion and r the mutation-load ratio (mean, or median, of
nonsynonymous mutations per exome sample, high vs low TMB).
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from common import GENE_COUNTS, PANEL_GENES, RESULTS, load_analytic

N_HIGH_POOLED, N_LOW_POOLED = 40, 3341
SELECT_MIN = 14   # genes shown in Figure 2: mutated in >= 14 of 39 high-TMB exomes

d = load_analytic()
ex = d[d.panel == 0]
hi, lo = ex[ex.high == 1], ex[ex.high == 0]
N_HIGH, N_LOW = len(hi), len(lo)
r_mean = hi.N_NONSYN_MUT.mean() / lo.N_NONSYN_MUT.mean()
r_median = hi.N_NONSYN_MUT.median() / lo.N_NONSYN_MUT.median()
print(f"Exome high/low: {N_HIGH}/{N_LOW}")
print(f"Mean mutations {hi.N_NONSYN_MUT.mean():.1f} vs {lo.N_NONSYN_MUT.mean():.1f}, r = {r_mean:.1f}")
print(f"Median mutations {hi.N_NONSYN_MUT.median():.0f} vs {lo.N_NONSYN_MUT.median():.0f}, r = {r_median:.1f}")

g = pd.read_csv(GENE_COUNTS)
panel = set(pd.read_csv(PANEL_GENES)["simbolo_usado"].dropna())
g = g[~g.GENE.isin(panel)].copy()
print("Non-panel genes tested:", len(g))


def fisher_or(a, n1, c, n0):
    b, dd = n1 - a, n0 - c
    p = fisher_exact([[a, b], [c, dd]])[1]
    or_ = (a * dd) / (b * c) if b * c > 0 else np.inf
    return or_, p


res = [fisher_or(a, N_HIGH, c, N_LOW) for a, c in zip(g.N_HIGH, g.N_LOW)]
g["OR_exome"] = [x[0] for x in res]
g["p_exome"] = [x[1] for x in res]
g["q_exome"] = multipletests(g.p_exome, method="fdr_bh")[1]
g["OR_pooled_repository"] = [fisher_or(a, N_HIGH_POOLED, c, N_LOW_POOLED)[0]
                             for a, c in zip(g.N_HIGH, g.N_LOW)]

g["freq_low"] = g.N_LOW / N_LOW
g["freq_obs_high"] = g.N_HIGH / N_HIGH
g["freq_exp_mean"] = 1 - (1 - g.freq_low) ** r_mean
g["freq_exp_median"] = 1 - (1 - g.freq_low) ** r_median
g["obs_exp_mean"] = g.freq_obs_high / g.freq_exp_mean
g["obs_exp_median"] = g.freq_obs_high / g.freq_exp_median

sig = g[g.q_exome < 0.05]
print(f"q < 0.05: {len(sig)} genes; at or below the mean-based null: {(sig.obs_exp_mean <= 1).sum()}")

for gene in ["TTN", "MUC16"]:
    x = g[g.GENE == gene].iloc[0]
    print(f"{gene}: OR pooled (repository) {x.OR_pooled_repository:.1f} -> OR exome {x.OR_exome:.1f}")

top = g[g.N_HIGH >= SELECT_MIN].sort_values("N_HIGH", ascending=False)
print(f"\nGenes mutated in >= {SELECT_MIN}/{N_HIGH} high-TMB exomes: {len(top)}")
cols = ["GENE", "N_HIGH", "N_LOW", "freq_obs_high", "freq_exp_mean", "obs_exp_mean", "obs_exp_median"]
print(top[cols].round(3).to_string(index=False))
print("Exceed mean-based null:", int((top.obs_exp_mean > 1).sum()),
      f"(ratios {top[top.obs_exp_mean > 1].obs_exp_mean.min():.2f}-{top.obs_exp_mean.max():.2f})")
print("Exceed median-based null:", int((top.obs_exp_median > 1).sum()))

top[cols].to_csv(RESULTS / "tableS1_top_genes.csv", index=False)
g.to_csv(RESULTS / "dataS1_all_genes.csv", index=False)
