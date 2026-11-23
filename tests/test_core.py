import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from omicsmarker.cohort import synthetic, prepare, integrate
from omicsmarker.raw import summarize
from omicsmarker.train import transform


def test_alignment_imputation_and_concordance(tmp_path):
    cfg = yaml.safe_load(Path("config/test.yaml").read_text())
    synthetic(tmp_path, cfg["seed"])
    out = tmp_path / "ready"
    out.mkdir()
    prepare(tmp_path, out, cfg)
    a = pd.read_parquet(out / "protein.parquet")
    prepare(tmp_path, out, cfg)
    pd.testing.assert_frame_equal(a, pd.read_parquet(out / "protein.parquet"))
    assert not a.isna().any().any()
    assert a.index.equals(pd.read_parquet(out / "labels.parquet").index)
    integrate(out, out)
    assert len(pd.read_csv(out / "correlation.tsv", sep="\t")) == 300


def test_held_out_rows_cannot_change_selector_or_scaler():
    x = pd.DataFrame({"g": [1, 4, 3, 2, 10], "h": [0, 0, 0, 0, 500]})
    training, evaluation = np.arange(4), [4]
    args = (training, evaluation, {"rna": 1})
    _, _, a, names = transform({"rna": x}, *args)
    x.loc[4, "h"] = -8e6
    _, _, b, _ = transform({"rna": x}, *args)
    assert a == b and names == ["rna:g"] and a["rna"]["mean"] == [2.5]


def test_sage_removes_decoys_and_shared_peptides():
    psms = pd.DataFrame(
        {
            "peptide": list("abcde"),
            "proteins": ["P", "rev_P", "Q", "R", "P;Q"],
            "label": [1, -1, 1, 1, 1],
            "peptide_q": [0, 0, 0.1, 0, 0],
            "protein_q": [0, 0, 0, 0.2, 0],
            "filename": ["file"] * 5,
        }
    )
    lfq = pd.DataFrame({"peptide": list("abcde"), "file": [4, 5, 8, 9, 10]})
    matrix, qc = summarize(psms, lfq)
    assert matrix.index.tolist() == ["P"] and matrix.loc["P", "file"] == 0
    assert qc["file"] == {"psms": 2, "peptides": 2, "unique_proteins": 1}
