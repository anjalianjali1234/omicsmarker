#!/usr/bin/env bash
set -euo pipefail
d=${1:-runs/demo}
export TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1 OMP_NUM_THREADS=1
omicsmarker synthetic --out "$d/input"
omicsmarker prepare --data "$d/input" --config config/test.yaml --out "$d/ready"
omicsmarker integrate --data "$d/ready" --out "$d/integration"
omicsmarker train --data "$d/ready" --config config/test.yaml --target ER --git-commit "$(git rev-parse HEAD 2>/dev/null || echo source-archive)" --out "$d/train"
omicsmarker report --analysis "$d/integration" --records "$d/train/run_record.json" --predictions "$d/train/predictions.json" --out "$d/report.html"
