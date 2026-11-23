# Data access and interpretation

Install the CPTAC extra, then export a local cohort with the commands below. Keep downloaded datasets outside the source tree. The exporter requests Broad RNA, UMich proteomics and MSSM clinical data through cptac 1.5.14. It sums Broad transcript TPM to genes before log transformation. The cBioPortal BRCA study sample and patient tables supply receptor calls, PAM50 and TMT plex. Both files are checksum-verified. Equivocal receptor calls remain missing. config/run.yaml uses clinical labels from this export.

```bash
pip install -e '.[data]'
omicsmarker fetch --out ~/om_data/cptac
omicsmarker prepare --data ~/om_data/cptac --config config/run.yaml --out ~/om_data/ready
omicsmarker integrate --data ~/om_data/ready --out ~/om_data/integration
```

An earlier run recorded 119 matched tumours and 9,827 shared genes with median Spearman correlation 0.398. Those measurements were not rerun for the rebuilt download. Recorded class counts were ER 80 positive / 37 negative; PR 68 / 44; HER2 13 / 82. PAM50 counts were LumA 56, Basal 27, LumB 17, Her2 14 and Normal-like 5. The project minimum of 15 per class was missed for HER2 and PAM50. Review these counts, any rare-class exclusions and batch plots for each export.

PRIDE PXD010785 contains label-free DDA breast-cancer cell-line data. The earlier project used the first eight RAW filenames in sorted order as its proposed demonstration subset of sixteen. Use the reviewed human Swiss-Prot reference (UP000005640), record its checksum and let Sage generate decoys. Verify tolerances against the instrument and compare search QC with the source study before using the results. No RAW files, clinical rows or patient-level predictions are distributed.

[PRIDE PXD010785](https://www.ebi.ac.uk/pride/archive/projects/PXD010785) · [UniProt reference](https://www.uniprot.org/proteomes/UP000005640) · [Study annotations](https://github.com/cBioPortal/datahub/tree/master/public/brca_cptac_2020) · [Sage](https://github.com/lazear/sage)
