FROM mambaorg/micromamba:2.0.5
RUN micromamba install -y -n base -c conda-forge -c bioconda sage-proteomics=0.14.7 thermorawfileparser=1.4.5 && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
