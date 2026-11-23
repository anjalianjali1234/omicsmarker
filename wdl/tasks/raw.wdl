version 1.0
task convert {
  input { File raw
          String image }
  command <<<
    ThermoRawFileParser.sh -i='~{raw}' -o=. -f=2
  >>>
  output { File mzml = glob("*.mzML")[0] }
  runtime { docker: image }
}
task search {
  input { Array[File] mzml
          File fasta
          File settings
          String image }
  command <<<
    set -eu
    sage -f '~{fasta}' -o result '~{settings}' '~{sep="' '" mzml}'
  >>>
  output { File psms="result/results.sage.tsv"
           File lfq="result/lfq.tsv" }
  runtime { docker:image
            cpu:4
            memory:"16 GB" }
}
task quantify {
  input { File psms
          File lfq
          File fasta
          String image }
  command <<<
    omicsmarker sage-to-protein --psms '~{psms}' --lfq '~{lfq}' --fasta '~{fasta}' --out result
  >>>
  output { File matrix="result/protein_matrix.parquet"
           File genes="result/gene_matrix.parquet"
           File qc="result/search_qc.json" }
  runtime { docker:image }
}
