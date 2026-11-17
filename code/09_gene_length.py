"""Baseline mutation frequency versus coding-sequence length (Figure S1).

For the large passenger genes shown in Figure 2 (mutated in >= 14 of the 39
high-TMB exomes), the baseline frequency (proportion of low-TMB exome samples
with a nonsynonymous mutation) is related to coding-sequence length.
Coding-sequence length is taken from the canonical Ensembl transcript,
3 x (protein length + 1), retrieved from the Ensembl REST API. The length
source is external, so the coefficient may vary slightly with the Ensembl
release.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import pearsonr, spearmanr

from common import DATA, RESULTS

top = pd.read_csv(RESULTS / "tableS1_top_genes.csv")      # run 05_gene_null.py first
cache = DATA / "cds_length_ensembl.csv"

if cache.exists():
    lengths = pd.read_csv(cache)
else:
    out = []
    for g in top.GENE:
        r = requests.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}",
                         params={"expand": 1}, headers={"Accept": "application/json"}, timeout=600)
        r.raise_for_status()
        t = [t for t in r.json().get("Transcript", []) if t.get("is_canonical") == 1]
        if t and t[0].get("Translation"):
            out.append({"GENE": g, "cds_length": 3 * (t[0]["Translation"]["length"] + 1)})
    lengths = pd.DataFrame(out)
    lengths.to_csv(cache, index=False)

m = top.merge(lengths, on="GENE")
m["freq_low"] = m["N_LOW"] / 1424
rho, p = spearmanr(m.cds_length, m.freq_low)
r_log, p_log = pearsonr(np.log(m.cds_length), m.freq_low)
print(f"Genes: {len(m)}; Spearman rho = {rho:.2f} (p = {p:.2g}); Pearson on log length r = {r_log:.2f}")
print(m[["GENE", "cds_length", "freq_low"]].round(3).to_string(index=False))

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(m.cds_length, m.freq_low * 100, color="black", s=14)
for _, x in m.iterrows():
    ax.annotate(x.GENE, (x.cds_length, x.freq_low * 100), fontsize=6, xytext=(3, 2),
                textcoords="offset points")
b = np.polyfit(np.log(m.cds_length), m.freq_low * 100, 1)
xs = np.linspace(np.log(m.cds_length).min(), np.log(m.cds_length).max(), 50)
ax.plot(np.exp(xs), np.polyval(b, xs), color="black", linewidth=0.8)
ax.set_xscale("log")
ax.set_xlabel("Coding-sequence length (bp, log scale)")
ax.set_ylabel("Mutated in low-TMB exomes (%)")
ax.set_title(f"Spearman rho = {rho:.2f}")
fig.tight_layout()
fig.savefig(RESULTS / "figureS1_length.png", dpi=300)
