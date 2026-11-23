cwlVersion: v1.2
class: CommandLineTool
baseCommand: sage
requirements:
  DockerRequirement: {dockerPull: 'omicsmarker-proteomics:dev'}
  ResourceRequirement: {coresMin: 4, ramMin: 16000}
arguments:
  - {position: 2, prefix: -o, valueFrom: result}
inputs:
  fasta: {type: File, inputBinding: {prefix: -f, position: 1}}
  config: {type: File, inputBinding: {position: 3}}
  mzml: {type: 'File[]', inputBinding: {position: 4}}
outputs:
  psms: {type: File, outputBinding: {glob: result/results.sage.tsv}}
  lfq: {type: File, outputBinding: {glob: result/lfq.tsv}}
