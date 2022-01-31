# Don't upgrade nfcore/base, it creates "Kernel too old" error for singularity (because of the debian image)
FROM nfcore/base:latest

LABEL author="alper.kucukural@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/barcodseq pipeline"

RUN apt-get update && apt-get install -y gcc libltdl7 libtbb-dev

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-barcodeseq-1.0/bin:$PATH


# Now install R and littler, and create a link for littler in /usr/local/bin
# Default CRAN repo is now set by R itself, and littler knows about it too
# r-cran-docopt is not currently in c2d4u so we install from source
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
                 littler \
 		 r-base \
 		 r-base-dev \
 		 r-recommended \
		 libfreetype6-dev libfontconfig1-dev libudunits2-dev libssl-dev \
		 libcurl4-openssl-dev libxml2-dev \

  	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
 	&& install.r docopt \
 	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
 	&& rm -rf /var/lib/apt/lists/*

COPY install_packages.R /
RUN Rscript /install_packages.R
