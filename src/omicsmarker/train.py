"""Fit each feature transformation using only the training partition."""

import hashlib
import json
import os
from datetime import datetime, timezone
from importlib.metadata import version

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import ParameterGrid, StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.utils.class_weight import compute_class_weight


def transform(data, training, evaluation, limits):
    train, test, state, names = [], [], {}, []
    for assay, table in data.items():
        cols = table.iloc[training].var().nlargest(limits[assay]).index
        array = table[cols].to_numpy()
        scale = StandardScaler().fit(array[training])
        train.append(scale.transform(array[training]))
        test.append(scale.transform(array[evaluation]))
        names.extend(f"{assay}:{g}" for g in cols)
        state[assay] = {
            "columns": cols.tolist(),
            "mean": scale.mean_.tolist(),
            "scale": scale.scale_.tolist(),
        }
    return train, test, state, names


def build(widths, nclasses, units, dropout, seed):
    import tensorflow as tf

    tf.keras.backend.clear_session()
    tf.keras.utils.set_random_seed(seed)
    tf.config.experimental.enable_op_determinism()
    inp = [tf.keras.Input((n,)) for n in widths]
    reg = tf.keras.regularizers.l2(1e-4)

    def dense(x, n):
        return tf.keras.layers.Dropout(dropout)(
            tf.keras.layers.Dense(n, activation="relu", kernel_regularizer=reg)(x)
        )

    branches = [dense(x, units) for x in inp]
    merged = (
        branches[0] if len(branches) == 1 else tf.keras.layers.Concatenate()(branches)
    )
    binary = nclasses == 2
    score = tf.keras.layers.Dense(
        1 if binary else nclasses, "sigmoid" if binary else "softmax"
    )(dense(merged, units // 2))
    model = tf.keras.Model(inp, score)
    model.compile(
        optimizer="adam",
        loss="binary_crossentropy" if binary else "sparse_categorical_crossentropy",
    )
    return model


def auc(y, p):
    return float(roc_auc_score(y, p[:, 1] if p.shape[1] == 2 else p, multi_class="ovr"))


def fit(data, y, tr, te, kind, hyp, cfg):
    seed = cfg["seed"]
    if kind == "tensorflow":
        import tensorflow as tf

        tr, val = train_test_split(tr, test_size=0.2, stratify=y[tr], random_state=seed)
    x, z, state, names = transform(data, tr, te, cfg["features"])
    if kind == "xgboost":
        from xgboost import XGBClassifier

        model = XGBClassifier(
            **hyp, n_jobs=1, random_state=seed, subsample=0.8, colsample_bytree=0.5
        )
        model.fit(np.concatenate(x, axis=1), y[tr])
        prob = model.predict_proba(np.concatenate(z, axis=1))
    elif kind == "tensorflow":
        classes = np.unique(y[tr])
        model = build([v.shape[1] for v in x], len(classes), seed=seed, **hyp)
        _, valid, _, _ = transform(data, tr, val, cfg["features"])
        weights = compute_class_weight("balanced", classes=classes, y=y[tr])
        model.fit(
            x,
            y[tr],
            validation_data=(valid, y[val]),
            class_weight=dict(zip(classes, weights)),
            verbose=0,
            epochs=cfg["epochs"],
            batch_size=16,
            callbacks=[
                tf.keras.callbacks.EarlyStopping(patience=20, restore_best_weights=True)
            ],
        )
        prob = model.predict(z, verbose=0)
        if prob.shape[1] == 1:
            prob = np.column_stack([1 - prob[:, 0], prob[:, 0]])
    else:
        raise ValueError(f"Unknown model {kind}")
    return prob, model, x, z, state, names


def importance(kind, model, x, z, names, seed):
    import shap

    if kind == "xgboost":
        blocks = [shap.TreeExplainer(model).shap_values(np.concatenate(z, axis=1))]
    else:
        sel = np.random.default_rng(seed).choice(
            len(x[0]), min(100, len(x[0])), replace=False
        )
        bg = [arr[sel] for arr in x]
        ex = shap.GradientExplainer(model, bg if len(bg) > 1 else bg[0])
        values = ex.shap_values(z if len(z) > 1 else z[0], rseed=seed)
        blocks = values if isinstance(values, list) else [values]
    scores = []
    for block in blocks:
        vals = np.abs(np.asarray(block)).mean(axis=0)
        scores.extend(vals if vals.ndim == 1 else vals.mean(axis=-1))
    return pd.Series(scores, index=names)


def train(data, out, cfg, target, revision):
    matrices = {k: pd.read_parquet(data / (k + ".parquet")) for k in ("rna", "protein")}
    labels = pd.read_parquet(data / "labels.parquet")[target]
    if any(not labels.index.equals(v.index) for v in matrices.values()):
        raise ValueError("Unaligned sample order")
    labels = labels.dropna().astype(str)
    counts = labels.value_counts()
    labels = labels[labels.isin(counts[counts >= cfg["min_class"]].index)]
    enc = LabelEncoder()
    y = enc.fit_transform(labels)
    nclasses = len(enc.classes_)
    if nclasses < 2 or (target != "PAM50" and nclasses != 2):
        raise ValueError(f"Review {target} classes: {counts.to_dict()}")
    outer = StratifiedKFold(cfg["outer_folds"], shuffle=True, random_state=cfg["seed"])
    inner = StratifiedKFold(cfg["inner_folds"], shuffle=True, random_state=cfg["seed"])
    metrics, preds, top = [], [], []
    for modality in ("rna", "protein", "combined"):
        tables = {
            k: v.loc[labels.index]
            for k, v in matrices.items()
            if modality in (k, "combined")
        }
        for kind in cfg["models"]:
            scores = []
            for fold, (tr, te) in enumerate(outer.split(y, y)):
                tuning = []
                for hyp in ParameterGrid(cfg[kind]):
                    perf = [
                        auc(y[tr[b]], fit(tables, y, tr[a], tr[b], kind, hyp, cfg)[0])
                        for a, b in inner.split(tr, y[tr])
                    ]
                    tuning.append((float(np.mean(perf)), hyp))
                hyp = max(tuning, key=lambda p: p[0])[1]
                prob, model, x, z, state, names = fit(tables, y, tr, te, kind, hyp, cfg)
                prefix = f"{modality}-{kind}-{fold}"
                if kind == "xgboost":
                    model.save_model(out / (prefix + ".json"))
                else:
                    model.save(out / (prefix + ".keras"))
                (out / (prefix + "-scale.json")).write_text(json.dumps(state))
                meta = {
                    "target": target,
                    "features": modality,
                    "model": kind,
                    "fold": fold,
                }
                metrics.append(
                    {
                        **meta,
                        "auroc": auc(y[te], prob),
                        "balanced_accuracy": float(
                            balanced_accuracy_score(y[te], prob.argmax(axis=1))
                        ),
                        "parameters": hyp,
                    }
                )
                for j, i in enumerate(te):
                    preds.append(
                        {
                            **meta,
                            "sample": str(labels.index[i]),
                            "truth": int(y[i]),
                            "probabilities": prob[j].tolist(),
                        }
                    )
                scores.append(importance(kind, model, x, z, names, cfg["seed"]))
            ranked = pd.concat(scores, axis=1).fillna(0).mean(axis=1).nlargest(20)
            top.extend(
                {
                    "target": target,
                    "features": modality,
                    "model": kind,
                    "feature": g,
                    "mean_abs_shap": v,
                }
                for g, v in ranked.items()
            )
    (out / "predictions.json").write_text(json.dumps(preds))
    rec = {
        "git_commit": revision,
        "target": target,
        "classes": enc.classes_.tolist(),
        "config": cfg,
        "counts": labels.value_counts().to_dict(),
        "excluded_classes": counts[counts < cfg["min_class"]].to_dict(),
        "metrics": metrics,
        "top_features": top,
        "versions": {k: version(k) for k in ["scikit-learn", "shap", *cfg["models"]]},
        "artifacts": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in out.iterdir()
            if p.is_file() and p.name != "run_record.json"
        },
        "logged_at": datetime.now(timezone.utc).isoformat(),
    }
    uri = os.environ.get(cfg["mongo_variable"])
    rec["run_id"] = "local"
    if uri:
        from pymongo import MongoClient

        with MongoClient(uri, serverSelectionTimeoutMS=5000) as client:
            coll = client.omicsmarker.runs
            coll.create_index([("target", 1), ("logged_at", -1)])
            rec["run_id"] = str(coll.insert_one(dict(rec)).inserted_id)
    (out / "run_record.json").write_text(json.dumps(rec, indent=2, allow_nan=False))
