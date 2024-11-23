This repository contains all the codes necessary to reproduce the modeling result contained in the article entitles '...'. Docker file contains a minimal script necessary to create an isolated environment for testing the scripts.

**Docker Image for R, Python**

This repository provides a Docker image based on Ubuntu 20.04, pre-configured with:

a) R (version 4.x) with libraries necessary for the analysis
    
b) Python (version 3.9) with libraries necessary for the analysis
    
c) Support for installing R packages such as ggpubr and Python libraries
    
d) Essential system dependencies for compiling R and Python packages

**Getting started**

**1) Install Docker**

Instructions to install Docker is available here: https://docs.docker.com/get-started/

**2) Build locally or pull the image**

To build it locally, download this repository and navigate to this folder and run 

    docker build -t <image_name> .

To pull the image from docker hub, run

    docker pull arnabrkmv/transcriptional_noise:latest

**3) Run the container**

To start the container interactively:

    docker run -it -p 8787:8787 my-r-python-rstudio-image

For persistent data storage, mount a local directory:

    docker run -it -p 8787:8787 -v /path/to/local/dir:/home/rstudio my-r-python-rstudio-image



