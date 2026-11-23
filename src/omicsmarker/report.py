import base64
import io
import json
import matplotlib.pyplot as plt
import pandas as pd
from jinja2 import Environment
from sklearn.metrics import RocCurveDisplay, ConfusionMatrixDisplay


def image():
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    plt.close()
    return base64.b64encode(buf.getvalue()).decode()


def report(analysis, records, predictions, out, qc=None):
    runs = [json.loads(p.read_text()) for p in records]
    metrics = pd.DataFrame([m for r in runs for m in r["metrics"]])
    summary = (
        metrics.groupby(["target", "features", "model"])[["auroc", "balanced_accuracy"]]
        .agg(["mean", "std"])
        .round(3)
    )
    ranks = pd.DataFrame([m for r in runs for m in r["top_features"]])
    corr = pd.read_csv(analysis / "correlation.tsv", sep="\t")
    pred = pd.concat([pd.read_json(p) for p in predictions])
    curves = []
    bars = []
    for key, g in pred.groupby(["target", "features", "model"]):
        probs = pd.DataFrame(g.probabilities.tolist())
        _, ax = plt.subplots()
        for k in range(1 if probs.shape[1] == 2 else 0, probs.shape[1]):
            RocCurveDisplay.from_predictions(
                g.truth == k, probs[k], ax=ax, name=f"class {k}"
            )
        ax.set_title(" / ".join(key))
        curves.append(image())
        if key[0] == "PAM50":
            ConfusionMatrixDisplay.from_predictions(
                g.truth, probs.to_numpy().argmax(axis=1)
            )
            curves.append(image())
    for key, g in ranks.groupby(["target", "features", "model"]):
        ax = (
            g.nlargest(10, "mean_abs_shap")
            .sort_values("mean_abs_shap")
            .plot.barh(x="feature", y="mean_abs_shap", legend=False)
        )
        ax.set_title(" / ".join(key))
        bars.append(image())
    tpl = Environment(autoescape=True).from_string("""<!doctype html><meta charset="utf-8"><title>OmicsMarker</title>
<style>body{font:16px system-ui;max-width:76rem;margin:2rem auto;padding:0 1rem}img{max-width:100%}td,th{padding:.3rem}pre{white-space:pre-wrap}</style>
<h1>OmicsMarker research report</h1><h2>1. Data summary</h2><pre>{{counts}}\n{{qc}}</pre>
<h2>2. RNA–protein concordance</h2><p>Median Spearman correlation {{median}} across {{genes}} shared genes.</p>
<img alt="Correlation histogram" src="data:image/png;base64,{{hist}}">{{pathways|safe}}
<h2>3. Outer-fold metrics</h2>{{summary|safe}}
<h2>4. Held-out predictions</h2>{% for c in curves %}<img alt="Classification evaluation" src="data:image/png;base64,{{c}}">{% endfor %}
<h2>5. Feature importance</h2>{% for b in bars %}<img alt="Held-out SHAP importance" src="data:image/png;base64,{{b}}">{% endfor %}{{ranks|safe}}
<h2>6. Run details</h2><pre>{{details}}</pre>""")
    pathway = analysis / "hallmark.tsv"
    out.write_text(
        tpl.render(
            counts={r["target"]: r["counts"] for r in runs},
            qc=qc.read_text() if qc else "Raw search disabled",
            median=round(corr.rho.median(), 3),
            genes=len(corr),
            hist=base64.b64encode((analysis / "correlation.png").read_bytes()).decode(),
            pathways=pd.read_csv(pathway, sep="\t").head(20).to_html()
            if pathway.exists()
            else "",
            summary=summary.to_html(),
            curves=curves,
            bars=bars,
            ranks=ranks.to_html(index=False),
            details=json.dumps(
                [
                    {
                        k: r[k]
                        for k in (
                            "target",
                            "git_commit",
                            "run_id",
                            "classes",
                            "versions",
                            "excluded_classes",
                        )
                    }
                    | {"seed": r["config"]["seed"]}
                    for r in runs
                ],
                indent=2,
            ),
        )
    )
