"""Download the three breast cancer studies from the public cBioPortal API and
build the two derived tables used by the analysis scripts.

Outputs (written to data/, not distributed with the repository):
  analytic_table.csv      one row per sample
  gene_counts_pooled.csv  mutated samples per gene, repository high vs low TMB

The original extraction was run on 20 July 2026. The repository is updated
over time, so a new run may differ slightly from the published numbers.
"""
import time

import pandas as pd
import requests

from common import DATA, STUDIES

API = "https://www.cbioportal.org/api"
NONSYN = {
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Splice_Site", "Translation_Start_Site",
}
SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json", "Content-Type": "application/json"})


def get(path, **params):
    for attempt in range(5):
        r = SESSION.get(f"{API}{path}", params=params, timeout=300)
        if r.ok:
            return r.json()
        time.sleep(2 ** attempt)
    r.raise_for_status()


def post(path, body, **params):
    for attempt in range(5):
        r = SESSION.post(f"{API}{path}", params=params, json=body, timeout=600)
        if r.ok:
            return r.json()
        time.sleep(2 ** attempt)
    r.raise_for_status()


def clinical(study, kind, attrs):
    rows = get(f"/studies/{study}/clinical-data", clinicalDataType=kind,
               projection="SUMMARY", pageSize=10_000_000)
    df = pd.DataFrame(rows)
    df = df[df.clinicalAttributeId.isin(attrs)]
    key = "sampleId" if kind == "SAMPLE" else "patientId"
    return df.pivot_table(index=key, columns="clinicalAttributeId", values="value",
                          aggfunc="first").reset_index()


def subtype_tcga(pam50):
    return {"BRCA_LumA": "HR+/HER2-", "BRCA_LumB": "HR+/HER2-",
            "BRCA_Her2": "HER2+", "BRCA_Basal": "TNBC"}.get(pam50, "Desconhecido")


def subtype_msk(er, pr, her2_fish):
    if her2_fish == "Positive":
        return "HER2+"
    if her2_fish == "Negative":
        if "Positive" in (er, pr):
            return "HR+/HER2-"
        if er == "Negative" and pr == "Negative":
            return "TNBC"
    return "Desconhecido"


def os_event(status):
    if not isinstance(status, str):
        return None
    return 1.0 if status.startswith("1") else 0.0


tables, mutations = [], []
for study in STUDIES:
    print("Study:", study)
    samples = pd.DataFrame(get(f"/studies/{study}/samples", projection="SUMMARY", pageSize=100000))
    samples = samples[["sampleId", "patientId"]]
    s_attrs = ["TMB_NONSYNONYMOUS", "SAMPLE_TYPE", "ER_STATUS", "PR_STATUS", "HER2_FISH_STATUS"]
    p_attrs = ["OS_MONTHS", "OS_STATUS", "AGE", "AGE_AT_DIAGNOSIS", "SUBTYPE"]
    cs = clinical(study, "SAMPLE", s_attrs)
    cp = clinical(study, "PATIENT", p_attrs)
    d = samples.merge(cs, on="sampleId", how="left").merge(cp, on="patientId", how="left")
    for c in s_attrs + p_attrs:
        if c not in d:
            d[c] = None

    mut = pd.DataFrame(post(f"/molecular-profiles/{study}_mutations/mutations/fetch",
                            {"sampleListId": f"{study}_all"}, projection="DETAILED"))
    mut["GENE"] = mut["gene"].apply(lambda g: g["hugoGeneSymbol"]).str.upper()
    mut = mut[mut.mutationType.isin(NONSYN)]
    n = mut.groupby("sampleId").size().rename("N_NONSYN_MUT")
    d = d.merge(n, on="sampleId", how="left")
    d["N_NONSYN_MUT"] = d["N_NONSYN_MUT"].fillna(0).astype(int)
    mutations.append(mut[["sampleId", "GENE"]].drop_duplicates())

    d["studyId"] = study
    d["PLATFORM"] = "PANEL" if study == "breast_msk_2018" else "WES"
    d["TMB"] = pd.to_numeric(d["TMB_NONSYNONYMOUS"], errors="coerce")
    d["TMB_SOURCE"] = d["TMB"].notna().map({True: "portal", False: "calculado"})
    d["TMB"] = d["TMB"].fillna(d["N_NONSYN_MUT"] / 30.0)
    if study == "brca_tcga_pan_can_atlas_2018":
        d["SAMPLE_TYPE_STD"] = "Primario"
        d["AGE_STD"] = pd.to_numeric(d["AGE"], errors="coerce")
        d["SUBTYPE_STD"] = d["SUBTYPE"].map(subtype_tcga)
    elif study == "breast_msk_2018":
        d["SAMPLE_TYPE_STD"] = d["SAMPLE_TYPE"].map({"Primary": "Primario", "Metastasis": "Metastatico"})
        d["AGE_STD"] = None
        d["SUBTYPE_STD"] = [subtype_msk(e, p, h) for e, p, h in
                            zip(d.ER_STATUS, d.PR_STATUS, d.HER2_FISH_STATUS)]
    else:
        d["SAMPLE_TYPE_STD"] = "Metastatico"
        # age is released in 5-year bands ("41-45"); the lower bound is used
        d["AGE_STD"] = pd.to_numeric(d["AGE_AT_DIAGNOSIS"].astype(str).str.split("-").str[0], errors="coerce")
        d["SUBTYPE_STD"] = "Desconhecido"
    d["OS_MONTHS_STD"] = pd.to_numeric(d["OS_MONTHS"], errors="coerce")
    d["OS_EVENT_STD"] = d["OS_STATUS"].map(os_event)
    d["TMB_GROUP"] = (d["TMB"] >= 10).map({True: "TMB-High", False: "TMB-Low"})
    tables.append(d)

cols = ["sampleId", "patientId", "studyId", "PLATFORM", "AGE_STD", "SAMPLE_TYPE_STD",
        "SUBTYPE_STD", "OS_MONTHS_STD", "OS_EVENT_STD", "N_NONSYN_MUT", "TMB", "TMB_SOURCE",
        "TMB_GROUP"]
analytic = pd.concat(tables, ignore_index=True)[cols]
DATA.mkdir(exist_ok=True)
analytic.to_csv(DATA / "analytic_table.csv", index=False)
print("analytic_table.csv:", len(analytic), "samples,", analytic.patientId.nunique(), "patients")

# Mutated samples per gene, high vs low TMB with the repository values (pooled)
mut = pd.concat(mutations).merge(analytic[["sampleId", "TMB_GROUP"]], on="sampleId")
n_hi = int((analytic.TMB_GROUP == "TMB-High").sum())
n_lo = int((analytic.TMB_GROUP == "TMB-Low").sum())
g = (mut.groupby(["GENE", "TMB_GROUP"]).sampleId.nunique().unstack(fill_value=0)
     .rename(columns={"TMB-High": "N_HIGH", "TMB-Low": "N_LOW"}).reset_index())
g["PCT_HIGH"] = (100 * g.N_HIGH / n_hi).round(2)
g["PCT_LOW"] = (100 * g.N_LOW / n_lo).round(2)
g = g.sort_values("N_HIGH", ascending=False)
g.to_csv(DATA / "gene_counts_pooled.csv", index=False)
print(f"gene_counts_pooled.csv: {len(g)} genes (high {n_hi}, low {n_lo})")
