# Use Ubuntu 20.04 base image
FROM ubuntu:20.04

# Set environment variable to avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.9, R version 4, and necessary dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    wget \
    dirmngr \
    gnupg \
    ca-certificates && \
    \
    # Add CRAN GPG key and repository for R version 4
    # wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee /usr/share/keyrings/cran-archive-keyring.gpg > /dev/null && \
    # echo "deb [signed-by=/usr/share/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" | tee /etc/apt/sources.list.d/cran.list && \
    # apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    # gpg --keyserver pgp.mit.edu --recv-key 381BA480 && \
    # gpg -a --export 381BA480 | apt-key add - && \
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

    
# Install R packages
RUN Rscript -e "install.packages(c('ggridges', 'ggplot2', 'dplyr', 'reshape2', 'moments'), repos='https://cloud.r-project.org')"


# Verify installations
 RUN python3 --version && R --version

# Set working directory
WORKDIR /app

# Default command to keep the container running
CMD ["bash"]
