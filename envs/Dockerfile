FROM nfcore/base:1.9
LABEL authors="Maxime Borry" \
    description="Docker image containing python software requirements for the maxibor/kraken-nf pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name kraken-nf > kraken-nf.yml
ENV PATH /opt/conda/envs/kraken-nf/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name kraken-nf > kraken-nf.yml