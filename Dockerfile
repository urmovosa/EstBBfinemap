FROM nfcore/base
LABEL authors="urmovosa@ut.ee" \
      description="Docker image for running SuSiE fine-mapping analyses on GWAS summary statistics from biobank"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/FineMapPipeline/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN apt-get update && apt-get install -y gcc
RUN apt-get install -y build-essential
RUN apt-get install libz-dev
RUN git clone https://github.com/statgen/emeraLD.git
WORKDIR $HOME/emeraLD
RUN make
