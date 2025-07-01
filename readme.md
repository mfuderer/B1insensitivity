# Software scripts for B1-Insensitivity
## Introduction
This repository is a collection of julia-scripts that relate to the results as published in "Reducing sensitivity to B1+ in transient-state quantitative-MRI" by M. Fuderer
This collection consists of the following elements:
- Generating the RF flip-angle sequences
- 1-D simulation of the effect of B1+ variation 
- Generation of figures from the reconstructed data

## Installing and using Julia
This repository contains Julia scripts for generating and analyzing data. These scripts rely on external Julia packages to perform various tasks. If you're new to Julia, don't worry—this guide will help you set everything up.

### Requirements: To run the scripts, you need:
1. Julia installed on your computer (download it from https://julialang.org).
2. An internet connection to download the required packages.

### Setting Up the Environment
Julia uses "environments" to manage packages. This repository includes a Project.toml and Manifest.toml file, which specify all the packages and their versions needed to run the scripts. Follow these steps to set up the environment:

1. Open Julia: Start Julia from your terminal or Julia IDE (e.g., VS Code or Jupyter Notebook).
2. Activate the Environment
Navigate to the folder containing this repository. In Julia, type:    cd("path/to/this/repository")

Then activate the environment:

using Pkg
Pkg.activate(".")

This tells Julia to use the environment defined in this folder.

3. Install the Packages
Once the environment is activated, install the required packages by running:

Pkg.instantiate()

This will download and install all the packages specified in the Project.toml and Manifest.toml files.

### Running the Scripts
After setting up the environment, you can run the scripts. For example:

include("sequenceGenerationScript.jl")

Make sure to follow any specific instructions in the script comments.

### Notes
If you encounter any issues, ensure that the environment is activated (Pkg.activate(".") and Pkg.instantiate()) before running the scripts.
The scripts may take some time to run, depending on the task and your computer's performance.

## Generation of flip-angle sequences 
This is done by running the sequenceGenerationScript in folder SequenceGeneration
Note that, within the script, one has to manually cycle the variable "case" from 1 to 4. 
The processing time is vastly different between case=1 (Amplitude-only noise-optimized), which may take a minute, and case=4 (Amplitude+Phase, B1-opt), which may take up to a day.
The results should appear in the folder RFsequences (well, they may actually be there already, but if you want to reproduce these, you can move them to another location and generate them again using the script).
Importantly, this script has BLAKJac as a dependency.

## 1-D simulation of the effect of B1+ variation
This is done by running the script Simulation1D in folder Simulation
This script has BlochSimulators and MRSTATToolbox as a dependency

## Generation of figures from the reconstructed data
This is done by running the scripts FigureGenerationPhantom and FigureGenerationVolunteer in folder ScanAnalysis
IMPORTANT: there are prerequisites:
1. The data has to be downloaded separately! Motivation: it is a bit bulky (zipped size 174MB, unzipped about 1.5GB). The data is to be found in https://doi.org/10.5281/zenodo.14916494. It has to be unzipped into your favorite folder and this should result in three files, Phantom.jld2, Volunteer1.jld2 and Volunteer2.jld2. The name of that folder has to be edited into script FigureGenerationSetup, as `figurePars["dataFolder"] = "(the reference to that folder)"`
2. The HD-BET tool has to be installed from github.com/MIC-DKFZ/HD-BET; the file location of the shell script 'bet' has to be edited into FigureGenerationSetup.jl
3. The FSL tool has to be installed from https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/FAST.html; the file location of the shell script 'fast' has to be edited into FigureGenerationSetup.jl"
4. This script requires a CUDA capable GPU.

























