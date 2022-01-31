# Don't upgrade nfcore/base, it creates "Kernel too old" error for singularity (because of the debian image)
FROM nfcore/base:dev

LABEL author="alper.kucukural@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/barcodseq pipeline"

RUN apt-get update && apt-get install -y gcc libltdl7 libtbb-dev

COPY environment.yml /
RUN conda env create -n myenv -f /environment.yml && conda clean -a
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-barcodeseq-1.0/bin:$PATH

RUN echo "conda activate myenv" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

COPY install_packages.R /
RUN Rscript /install_packages.R
