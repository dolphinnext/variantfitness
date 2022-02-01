# Don't upgrade nfcore/base, it creates "Kernel too old" error for singularity (because of the debian image)
FROM nfcore/base:1.7 

LABEL author="alper.kucukural@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/barcodseq pipeline"

RUN apt-get update && apt-get install -y gcc libltdl7 libtbb-dev libcairo2-dev
RUN conda update -n base -c defaults conda
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda install -c bioconda mageck
RUN echo ". /opt/conda/etc/profile.d/conda.sh" > /root/.bashrc
RUN echo "conda activate dolphinnext" >> /root/.bashrc

RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext/bin:$PATH

COPY install_packages.R /
RUN Rscript /install_packages.R