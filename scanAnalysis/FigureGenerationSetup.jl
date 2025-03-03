## Code common to both figure generation scripts

using Statistics
using JLD2
using Printf
using QMRIColors
using PyPlot

# A bunch of parameters passed to and from functions; filled in here with values that are common for phantom and volunteer experiments
figurePars = Dict()
figurePars["seqDescription"] = [ "Amplitude-only \n noise-opt", 
                                 "Amplitude-only \n B1-opt", 
                                 "Amplitude+Phase \n noise-opt", 
                                 "Amplitude+Phase \n B1-opt"]
figurePars["recDescription"] = ["PreCalibrated", "MisCalibrated", "NonCalibrated"]
figurePars["maxIt"]         = 21
figurePars["comparisonSet"] = [1,2]
figurePars["sliceRange"]    = 1:20 
figurePars["gels"]          = 1:12

nRecons = length(figurePars["recDescription"])
nSeqs   = length(figurePars["seqDescription"])

# To be inserted: reference to ...
#   - data, should be unzipped from 10.5281/zenodo.14916494 into a folder of your choice
#   - BET tool, see https://github.com/MIC-DKFZ/HD-BET
#   - FAST tool, see https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/FAST.html
figurePars["dataFolder"] = 
    "/smb/user/mfuderer/DS-Data/Radiotherapie/Research/Project/MRSTAT/experiments/miha/Scandata_miha/Pub_regroup/reConvert/"
figurePars["betTool"] = "/usr/local/fsl/share/fsl/bin/bet"
figurePars["fastTool"] = "/usr/local/fsl/share/fsl/bin/fast"

@assert isfile(figurePars["dataFolder"]*"Phantom.jld2") 
"Before proceeding, unzip files from 10.5281/zenodo.14916494 and edit the file location into FigureGenerationSetup.jl"

@assert isfile(figurePars["betTool"]) 
"Before proceeding, install HD-BET from github.com/MIC-DKFZ/HD-BET and edit the file location into FigureGenerationSetup.jl"

@assert isfile(figurePars["fastTool"]) 
"Before proceeding, install FSL from https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/FAST.html and edit the file location into FigureGenerationSetup.jl"

# Plotting style settings
PyPlot.rc("font", family="serif")
PyPlot.rc("font", size=14)
