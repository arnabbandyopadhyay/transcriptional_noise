This repository contains all the codes necessary to reproduce the modeling result contained in the article entitles '...'. Docker file contains a minimal script necessary to create an isolated environment for testing the scripts.

**Docker Image for R, Python**

This repository provides a Docker image based on Ubuntu 20.04, pre-configured with:

a) R (version 4.x) with libraries necessary for the analysis
    
b) Python (version 3.9) with libraries necessary for the analysis
    
c) Support for installing R packages such as ggpubr and Python libraries
    
d) Essential system dependencies for compiling R and Python packages

# Getting started

## Environment setup

**1) Install Docker**

Instructions to install Docker is available here: https://docs.docker.com/get-started/

**2) Build locally or pull the image**
 Like
To build it locally, download this repository and navigate to this folder and run 

    docker build -t <image_name> .

To pull the image from docker hub, run

    docker pull arnabrkmv/transcriptional_noise:latest

**3) Run the container**

To start the container interactively:

    docker run -it transcriptional_noise

**4) Mount a directory**

For persistent data storage, mount a local directory:

    docker run -it -v /path/to/local/dir:/mnt/ transcriptional_noise

## Description of folders

Four folders with all the necessary codes are already copied within the container in /home/transcriptional_noise folder. Working diectory is /home/. Folders have self explainatory names, for example noGrowth_noMemory_Switching folder contains scripts that only considers stochastic switching from low glpD state to high glpD state. Likewise, Growth_Memory_Switching folder contains scripts that consider growth disparities between high and low glpD expressing cells, inheritance memory and stochastic switching phenomenon. All folders contain high_ino and low_ino folder, representing high inoculum and low_inoculum setting as described in the manuscript. A sample result of 20 simulations considering ~10,000 are contained within each folder. To generate results in the manuscript, the model needs to run for ~40,000 cells. This can be done by modifying the max_cells value in the params.py file. 

## Launching simulation





