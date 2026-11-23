import argparse
from pathlib import Path
import yaml
from . import cohort


def main():
    p = argparse.ArgumentParser(
        prog="omicsmarker", description="Paired tumour RNA and proteomics workflow"
    )
    subs = p.add_subparsers(dest="command", required=True)
    spec = {
        "synthetic": ["out"],
        "fetch": ["out"],
        "prepare": ["data", "config", "out"],
        "integrate": ["data", "out"],
        "train": ["data", "config", "out"],
        "report": ["analysis", "out"],
        "sage-to-protein": ["psms", "lfq", "fasta", "out"],
    }
    for command, keys in spec.items():
        sub = subs.add_parser(command)
        for key in keys:
            sub.add_argument("--" + key, type=Path, required=True)
        if command == "synthetic":
            sub.add_argument("--seed", type=int, default=42)
        if command == "integrate":
            sub.add_argument("--gmt", type=Path)
            sub.add_argument("--pride", type=Path)
        if command == "train":
            sub.add_argument("--target", required=True)
            sub.add_argument("--git-commit", required=True)
        if command == "report":
            for name in ("records", "predictions"):
                sub.add_argument("--" + name, nargs="+", type=Path, required=True)
            sub.add_argument("--qc", type=Path)
    a = p.parse_args()
    if a.command != "report":
        a.out.mkdir(parents=True, exist_ok=True)
    if hasattr(a, "config"):
        cfg = yaml.safe_load(a.config.read_text())
    if a.command == "synthetic":
        cohort.synthetic(a.out, a.seed)
    elif a.command == "fetch":
        cohort.fetch(a.out)
    elif a.command == "prepare":
        cohort.prepare(a.data, a.out, cfg)
    elif a.command == "integrate":
        cohort.integrate(a.data, a.out, a.gmt, a.pride)
    elif a.command == "train":
        from .train import train

        train(a.data, a.out, cfg, a.target, a.git_commit)
    elif a.command == "report":
        from .report import report

        report(a.analysis, a.records, a.predictions, a.out, a.qc)
    else:
        from .raw import raw

        raw(a.psms, a.lfq, a.fasta, a.out)
