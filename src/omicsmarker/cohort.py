"""Cohort acquisition, alignment, quality control and RNA–protein concordance."""

import hashlib
import io
import json
from importlib.metadata import version
from urllib.request import urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA


def synthetic(out, seed):
    rng = np.random.default_rng(seed)
    y = rng.permutation(np.repeat([0, 1], 30))
    latent = rng.normal(size=(60, 300))
    latent[:, :10] += 4 * y[:, None]
    ids = [f"S{i:03}" for i in range(60)]
    genes = [f"GENE{i:03}" for i in range(300)]
    for name in ("proteomics", "transcriptomics"):
        assay = pd.DataFrame(
            latent + rng.normal(size=latent.shape), index=ids, columns=genes
        )
        if name == "proteomics":
            assay = assay.mask(rng.random(assay.shape) < 0.05)
        assay.to_parquet(out / f"{name}.parquet")
    pd.DataFrame({"Sample_Tumor_Normal": "Tumor", "ER": y}, index=ids).to_parquet(
        out / "clinical.parquet"
    )


def fetch(out):
    import cptac

    cohort = cptac.Brca()
    protein = cohort.get_proteomics(source="umich")
    protein.columns = [
        "|".join(map(str, c)) if isinstance(c, tuple) else str(c)
        for c in protein.columns
    ]
    protein.to_parquet(out / "proteomics.parquet")
    rna = (
        cohort.get_transcriptomics(source="broad")
        .T.groupby(level="Name")
        .sum(min_count=1)
        .T
    )
    rna.to_parquet(out / "transcriptomics.parquet")
    clinical = cohort.get_clinical(source="mssm")
    # Use only the MSSM index. The study-specific annotation supplies targets.
    tables = {}
    hashes = {
        "sample": "f28568f1608c9f63dfeceed0ea301e29ee2dea0efbd88bd826fa2510434896cd",
        "patient": "d9b0dee7ee144105a3fdd37b3573a76558114ec1d8e6fa7262117657517a1275",
    }
    for level, expected in hashes.items():
        url = (
            "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_cptac_2020/data_clinical_"
            + level
            + ".txt"
        )
        with urlopen(url, timeout=90) as response:
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError(f"Review updated {level} annotation before using it.")
        tables[level] = pd.read_csv(io.BytesIO(raw), sep="\t", comment="#").set_index(
            "PATIENT_ID"
        )
    keys = {
        "ER_UPDATED_CLINICAL_STATUS": "ER",
        "PR_CLINICAL_STATUS": "PR",
        "ERBB2_UPDATED_CLINICAL_STATUS": "HER2",
    }
    labels = (
        tables["sample"][["PAM50", "TMT_PLEX", "CANCER TYPE"]]
        .join(tables["patient"][list(keys)], validate="one_to_one")
        .rename(columns=keys)
    )
    labels.index = labels.index.str.replace(r"^X(?=\d)", "", regex=True)
    for t in keys.values():
        labels[t] = (
            labels[t].str.lower().map({"positive": "Positive", "negative": "Negative"})
        )
    labels["Sample_Tumor_Normal"] = labels["CANCER TYPE"].map(
        {"Breast Cancer": "Tumor"}
    )
    labels.loc[labels.index.intersection(clinical.index)].to_parquet(
        out / "clinical.parquet"
    )
    (out / "sources.json").write_text(
        json.dumps(
            {
                "cptac": version("cptac"),
                "protein": "umich",
                "rna": "broad",
                "clinical": "mssm",
                "label_checksums": hashes,
            },
            indent=2,
        )
    )


def prepare(data, out, cfg):
    protein, rna, clinical = [
        pd.read_parquet(data / (k + ".parquet"))
        for k in ("proteomics", "transcriptomics", "clinical")
    ]
    if any(not t.index.is_unique for t in (protein, rna, clinical)):
        raise ValueError("Resolve duplicate sample IDs before modeling")
    clinical = clinical.loc[clinical[cfg["sample_type"]] == cfg["tumour_value"]]
    labels = clinical[list(cfg["labels"].values())].dropna(how="all")
    labels.columns = list(cfg["labels"])
    ids = labels.index.intersection(protein.index).intersection(rna.index)
    if len(ids) < 2:
        raise ValueError("Too few matched tumours")
    labels.loc[ids].to_parquet(out / "labels.parquet")
    protein = protein.loc[ids].copy()
    protein = protein.iloc[:, np.argsort(protein.isna().sum(), kind="stable")]
    protein.columns = protein.columns.str.split("|", regex=False).str[0]
    protein = protein.loc[:, ~protein.columns.duplicated()]
    protein = protein.loc[:, protein.notna().mean() >= 0.5].copy()
    rng = np.random.default_rng(cfg["seed"])
    for sample, row in protein.iterrows():
        if row.count() < 2:
            raise ValueError(f"Insufficient protein observations for {sample}")
        missing = row.index[row.isna()]
        protein.loc[sample, missing] = rng.normal(
            row.mean() - 1.8 * row.std(), 0.3 * row.std(), len(missing)
        )
    rna = rna.loc[ids]
    if cfg["rna_scale"] == "linear":
        if (rna < 0).any().any():
            raise ValueError("Negative values in linear RNA")
        rna = np.log2(1 + rna.loc[:, (rna > 0).mean() >= 0.2])
    elif cfg["rna_scale"] != "log":
        raise ValueError("Choose linear or log RNA")
    for name, frame in [("rna", rna), ("protein", protein)]:
        if not frame.columns.is_unique or not np.isfinite(frame.to_numpy()).all():
            raise ValueError(f"Invalid {name} matrix")
        frame.to_parquet(out / f"{name}.parquet")
        pcs = PCA(n_components=2, random_state=cfg["seed"]).fit_transform(frame)
        groups = (
            clinical.loc[ids, cfg["batch"]].fillna("Unknown")
            if cfg.get("batch")
            else pd.Series("All", index=ids)
        )
        fig, ax = plt.subplots()
        for group in groups.unique():
            pts = pcs[groups.eq(group)]
            ax.scatter(pts[:, 0], pts[:, 1], label=str(group))
        ax.set(xlabel="PC1", ylabel="PC2", title=name)
        ax.legend(fontsize=7)
        fig.savefig(out / f"{name}_pca.png", bbox_inches="tight")
        plt.close(fig)


def integrate(data, out, gmt=None, pride=None):
    r, p = [pd.read_parquet(data / (k + ".parquet")) for k in ("rna", "protein")]
    if not r.index.equals(p.index):
        raise ValueError("RNA and protein sample order must match")
    cor = pd.DataFrame(
        [(g, *spearmanr(r[g], p[g])) for g in r.columns.intersection(p.columns)],
        columns=["gene", "rho", "p"],
    ).sort_values("rho")
    cor.to_csv(out / "correlation.tsv", sep="\t", index=False)
    ax = cor.rho.hist(bins=30)
    ax.set(xlabel="Spearman correlation", ylabel="Gene count")
    ax.figure.savefig(out / "correlation.png", bbox_inches="tight")
    plt.close(ax.figure)
    if gmt:
        import gseapy

        gseapy.prerank(
            rnk=cor[["gene", "rho"]].dropna(), gene_sets=str(gmt), outdir=None, seed=42
        ).res2d.to_csv(out / "hallmark.tsv", sep="\t", index=False)
    if pride:
        pr = pd.read_parquet(pride).median(axis=1).rank().rename("PRIDE")
        ranks = pd.concat([pr, p.median().rank().rename("CPTAC")], axis=1, join="inner")
        ranks.to_csv(out / "overlap.tsv", sep="\t")
        if len(ranks):
            ax = ranks.plot.scatter("PRIDE", "CPTAC")
            ax.figure.savefig(out / "overlap.png")
            plt.close(ax.figure)
