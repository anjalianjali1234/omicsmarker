import json
import re
import numpy as np
import pandas as pd


def summarize(psms, lfq):
    hits = psms[(psms.label == 1) & (psms.peptide_q <= 0.01) & (psms.protein_q <= 0.01)]
    mappings = hits[["peptide", "proteins"]].drop_duplicates()
    unique = mappings[
        ~mappings.peptide.duplicated(False)
        & ~mappings.proteins.str.contains(";", regex=False)
    ]
    q = lfq.merge(unique, on="peptide", suffixes=("_lfq", ""), validate="many_to_one")
    if "q_value" in q:
        q = q[q.q_value <= 0.01]
    files = [k for k in lfq if k in set(psms.filename)]
    if not files:
        raise ValueError("No LFQ file headers match Sage PSM filenames")
    matrix = q.groupby("proteins")[files].sum(min_count=1)
    matrix = np.log2(matrix.where(matrix > 0))
    qc = hits.groupby("filename").agg(
        psms=("peptide", "size"), peptides=("peptide", "nunique")
    )
    qc["unique_proteins"] = hits.merge(unique).groupby("filename").proteins.nunique()
    return matrix - matrix.median(), qc.reindex(files).fillna(0).astype(int).to_dict(
        "index"
    )


def raw(psms, lfq, fasta, out):
    matrix, qc = summarize(pd.read_csv(psms, sep="\t"), pd.read_csv(lfq, sep="\t"))
    symbols = {}
    with fasta.open() as f:
        for line in f:
            if line.startswith(">") and (match := re.search(r"\bGN=(\S+)", line)):
                full = line[1:].split()[0]
                for token in full.split("|") + [full]:
                    symbols[token] = match[1]
    matrix.to_parquet(out / "protein_matrix.parquet")
    matrix.loc[matrix.index.isin(symbols)].rename(index=symbols).groupby(
        level=0
    ).median().to_parquet(out / "gene_matrix.parquet")
    (out / "search_qc.json").write_text(json.dumps(qc, indent=2))
