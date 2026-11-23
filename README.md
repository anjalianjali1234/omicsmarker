# OmicsMarker

OmicsMarker compares RNA and protein measurements in matched breast tumours. It trains classifiers for ER, PR, HER2 and PAM50 using RNA alone, protein alone and both assays together. A separate workflow searches independent PRIDE spectra with Sage. The overlap between the two studies is evaluated by gene symbol, never by patient identity.

## Run the example

Unzip the project, open a terminal in the folder containing this README, and use Python 3.11 or 3.12 on Linux:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[analysis]' pytest
pytest -q tests
bash tests/smoke.sh
```

Open `runs/demo/report.html` locally. A ready-to-open copy from the verification run is in [docs/example-report.html](docs/example-report.html). The example generates 60 artificial samples and 300 genes per assay with signal in ten genes. It uses one binary target and three outer folds to exercise the pipeline.

## Workflow

| Stage | Command | Result |
|---|---|---|
| Retrieve CPTAC tables and study annotations | omicsmarker fetch | Matched-source tables and source checksums |
| Select tumours and prepare assays | omicsmarker prepare | RNA, protein and label matrices plus PCA plots |
| Measure within-gene concordance | omicsmarker integrate | Spearman statistics, histogram and optional Hallmark enrichment or PRIDE overlap |
| Evaluate classifiers | omicsmarker train | Nested-fold predictions, model files, transformations, SHAP rankings and a run record |
| Summarize a Sage search | omicsmarker sage-to-protein | Protein and gene matrices and per-file FDR counts |
| Generate offline HTML | omicsmarker report | Data summary, concordance, scores, ROC plots, feature rankings and run details |

The full config uses five outer and three inner stratified folds. Each training partition fits the feature selector and scaler. The TensorFlow network uses separate assay inputs and stratified early stopping. Held-out SHAP values are averaged across folds, assigning zero to features excluded in a fold. Initial cohort-level missingness filtering precedes cross-validation. Use an independent cohort before biological interpretation. Research use only.

## Workflow engines

Build the local images and generate inputs first:

```bash
docker build -f docker/analysis.Dockerfile -t omicsmarker-analysis:dev .
docker build -f docker/proteomics.Dockerfile -t omicsmarker-proteomics:dev .
omicsmarker synthetic --out runs/demo/input
miniwdl check wdl/omicsmarker.wdl
miniwdl run wdl/omicsmarker.wdl -i wdl/inputs/test.json --dir runs/wdl
java -jar "$CROMWELL_JAR" run wdl/omicsmarker.wdl -i wdl/inputs/test.json
cwltool --validate cwl/sage_search.cwl
```

Edit full.json to use absolute paths and a Git revision. For raw search set run_raw_search to true, supply raw_files, fasta and sage_config. Use hallmark_gmt for a local gene-set file. Versions of published images can be supplied through analysis_image and proteomics_image. The standalone CWL task takes fasta and config File objects and an array of mzML File objects. Set OM_MONGO_URI in the training process to enable database tracking. Every successful target run writes run_record.json.

[Data sources and limitations](docs/data.md) · [Example metrics](docs/results.md) · [MIT license](LICENSE)
