
# Epigenetic cellular memory in Pseudomonas aeruginosa generates phenotypic variation in response to host environments

Elisabeth Vatareck 1,2⁕, Tim Rick 3⁕, Nicolas Oswaldo Gomez 1⁕, Arnab Bandyopadhyay 4*, Janina Kramer 1 , Dmytro Strunin 3 , Jelena Erdmann 1,2 , Oliver Hartmann 1,2 , Kathrin Alpers 1 , Christian Boedeker 5 , Christian Sieben 6,7 , Gang Zhao 4,8 , Jürgen Tomasch 1 , and Susanne Häussler 1,2,3,9¶

1 Department of Molecular Bacteriology, Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany 

2 Institute of Molecular Bacteriology, TWINCORE Centre for Experimental and Clinical Infection Research, a joint venture of the Hannover Medical School (MHH) and the HZI, Hannover, Germany 

3 Department of Clinical Microbiology, Copenhagen University Hospital – Rigshospitalet 2100 Copenhagen, Denmark

4 Department of Systems Immunology, Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany 

5 Leibniz Institute DSMZ - Deutsche Sammlung von Mikroorganismen und Zellkulturen, Braunschweig, Germany 

6 Nanoscale infection Biology, Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany 

7 Institute of Genetics, Technische Universität Braunschweig, Braunschweig, Germany 

8 Current address: Clinical Pharmacology and Quantitative Pharmacology, Clinical Pharmacology & Safety Sciences, R&D, Astrazeneca, Cambridge, UK 

9 Cluster of Excellence RESIST (EXC 2155), Hannover Medical School, Hannover, 30625, Germany

⁕ contributed equally

¶ corresponding author

---

This repository contains all the codes necessary to reproduce the modeling result contained in the article entitles 'Epigenetic cellular memory in Pseudomonas aeruginosa generates phenotypic variation in response to host environments'. Docker file contains a minimal script necessary to create an isolated environment for testing the scripts.

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

Building a container will take ~30-60 min depending on the machine configuration and connection speed. Alternatively, the image can be pulled from docker hub. To pull the image from docker hub, run

    docker pull arnabrkmv/transcriptional_noise:latest

**3) Run the container**

To start the container interactively:

    docker run -it transcriptional_noise

**4) Mount a directory**

For persistent data storage, mount a local directory:

    docker run -it -v /path/to/local/dir:/mnt/ transcriptional_noise

## Description of folders

Four folders with all the necessary codes are already copied within the container in /home/transcriptional_noise folder. Working diectory is /home/. Folders have self explainatory names, for example noGrowth_noMemory_Switching folder contains scripts that only considers stochastic switching from low glpD state to high glpD state. Likewise, Growth_Memory_Switching folder contains scripts that consider growth disparities between high and low glpD expressing cells, inheritance memory and stochastic switching phenomenon. All folders contain high_ino and low_ino folder, representing high inoculum and low_inoculum setting as described in the manuscript. Each folder contains cell.py, gillespie.py, params.py and main.py file. Cells are treated as agents and are controlled through cell.py file. All parameters are in the param.py file. Stochastic gillespie simulation framework for the glp operon is described in the gillespie.py file. Each cell agent contains the glp operon and the state of different componenents of it is therefore determined through stochastic simulation. All agents and functions are called through main.py file. A sample result of 20 simulations started with either low or high inoculum conditions with stopping criteria ~10,000 cells are contained within each folder. To generate results in the manuscript, the model needs to run for ~40,000 cells. This can be done by modifying the max_cells value in the params.py file. 

## Launching simulation

Before launching simulation, it is recommended to set the max_cells value in params.py file to a desired value. All required libraries are pre-installed in the container environment. Further the code uses multiprocessing library and multiple simulations can be launched in parallel. To set this to a desired value, modify the last couple of lines in main.py file. Generally, a container has no resource constraints and can use as much of a given resource as the host's kernel scheduler allows. For restricted resource allocation, check the guidelines at: https://docs.docker.com/engine/containers/resource_constraints/

Finally, to launch a simulation, navigate to a desired folder and from terminal run:

    nohup python main.py &

## Computing statistics

Rplot.R provides a template for generating population distributions and for comparing statistics. All required libraries are pre-installed in the container environment. Edit the locations of the directory in the script as necessary, modify my_comparisons list accordingly. It generates two plots: sample_plot.pdf that plots the distributions and stat_sample.pdf that plots computed statistics.  

    Rscript Rplot.R





