FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.9, R version 4, and necessary dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    wget \
    dirmngr \
    gnupg \
    make \
    cmake \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    liblapack-dev \
    libblas-dev \
    gfortran \
    build-essential \
    ca-certificates && \
    \
    # Add key and repository for R version
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    r-base && \
    \
    # Add PPA for Python 3.9
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    python3.9 \
    python3-pip && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.9 1 && \
    update-alternatives --config python3 && \
    \
    # Cleanup to reduce image size
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* 

# install necessary libraries in python
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Install remotes package for version management
RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Install specific versions of R packages
RUN Rscript -e "remotes::install_version('ggridges', version = '0.5.6', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_version('dplyr', version = '1.1.4', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_version('ggplot2', version = '3.5.1', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_version('reshape2', version = '1.4.4', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_version('moments', version = '0.14.1', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_version('ggpubr', version = '0.6.0', repos = 'https://cloud.r-project.org')" 

# To install current version of R packages
# RUN Rscript -e "install.packages(c('ggridges', 'ggplot2', 'dplyr', 'reshape2', 'moments','ggpubr'), repos='https://cloud.r-project.org')"


# Verify installations
 RUN python3 --version && R --version

# Set the working directory to the user's home directory
WORKDIR /home
RUN mkdir transcriptional_noise
COPY . /home/transcriptional_noise/

# Default command to keep the container running
CMD ["bash"]
