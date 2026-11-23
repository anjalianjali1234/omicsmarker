version 1.0
import "tasks/raw.wdl" as raw
import "tasks/analysis.wdl" as a
workflow omicsmarker {
  input {
    File proteomics_table
    File transcriptomics_table
    File clinical_table
    File run_config
    String git_commit
    Array[String]+ targets=["ER","PR","HER2","PAM50"]
    Boolean run_raw_search=false
    Array[File] raw_files=[]
    File? fasta
    File? sage_config
    File? hallmark_gmt
    String analysis_image="omicsmarker-analysis:dev"
    String proteomics_image="omicsmarker-proteomics:dev"
  }
  if (run_raw_search) {
    scatter (f in raw_files) {
      call raw.convert { input: raw=f,image=proteomics_image }
    }
    call raw.search { input: mzml=convert.mzml,fasta=select_first([fasta]),settings=select_first([sage_config]),image=proteomics_image }
    call raw.quantify { input: psms=search.psms,lfq=search.lfq,fasta=select_first([fasta]),image=analysis_image }
  }
  call a.prepare { input: rna=transcriptomics_table,protein=proteomics_table,clinical=clinical_table,config=run_config,image=analysis_image }
  call a.integrate { input: prepared=prepare.prepared,gmt=hallmark_gmt,pride=quantify.genes,image=analysis_image }
  scatter (t in targets) {
    call a.train { input: prepared=prepare.prepared,config=run_config,target=t,revision=git_commit,image=analysis_image }
  }
  call a.report { input: integrated=integrate.integrated,qc=quantify.qc,records=train.record,predictions=train.predictions,image=analysis_image }
  output { File html=report.html
           File integration=integrate.integrated
           Array[File] records=train.record
           Array[File] predictions=train.predictions
           Array[File] models=train.models
           File? raw_proteins=quantify.matrix
           File? search_qc=quantify.qc }
}
