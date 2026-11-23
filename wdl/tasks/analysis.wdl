version 1.0
task prepare {
  input { File rna
          File protein
          File clinical
          File config
          String image }
  command <<<
    set -eu
    mkdir input
    ln -s '~{rna}' input/transcriptomics.parquet
    ln -s '~{protein}' input/proteomics.parquet
    ln -s '~{clinical}' input/clinical.parquet
    omicsmarker prepare --data input --config '~{config}' --out ready
    tar czf ready.tar.gz ready
  >>>
  output { File prepared="ready.tar.gz" }
  runtime { docker:image }
}
task integrate {
  input { File prepared
          File? gmt
          File? pride
          String image }
  command <<<
    set -eu
    tar xzf '~{prepared}'
    omicsmarker integrate --data ready --out result ~{if defined(gmt) then "--gmt" else ""} ~{default="" gmt} ~{if defined(pride) then "--pride" else ""} ~{default="" pride}
    tar czf integration.tar.gz result
  >>>
  output { File integrated="integration.tar.gz" }
  runtime { docker:image }
}
task train {
  input { File prepared
          File config
          String target
          String revision
          String image }
  command <<<
    set -eu
    export OMP_NUM_THREADS=1 TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1
    tar xzf '~{prepared}'
    omicsmarker train --data ready --config '~{config}' --target '~{target}' --git-commit '~{revision}' --out trained
    tar czf models.tar.gz trained
  >>>
  output { File record="trained/run_record.json"
           File predictions="trained/predictions.json"
           File models="models.tar.gz" }
  runtime { docker:image
            memory:"8 GB" }
}
task report {
  input { File integrated
          File? qc
          Array[File] records
          Array[File] predictions
          String image }
  command <<<
    set -eu
    tar xzf '~{integrated}'
    omicsmarker report --analysis result --records '~{sep="' '" records}' --predictions '~{sep="' '" predictions}' ~{if defined(qc) then "--qc" else ""} ~{default="" qc} --out report.html
  >>>
  output { File html="report.html" }
  runtime { docker:image }
}
