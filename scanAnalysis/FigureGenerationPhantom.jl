## This script generates figures for the phantom data

include("FigureGenerationSetup.jl")
include("figuresB1phantom.jl")

figurePars["dispMax"]       = [2.0,0.5,1.0]
figurePars["diffuCorrFile"] = "scanAnalysis/mapCol240812A.jld2"
figurePars["filename"] = "Phantom"
figuresReadCompactData!(figurePars)
phantomVialMeans!(figurePars)
figureSequence2by2(figurePars)
figuresPhantomB1map(figurePars)
PyPlot.rc("font", size=9)
figuresPhantomBias(figurePars, true); # with diffusion correction
PyPlot.rc("font", size=14)
figureAllIm(figurePars)
figureDifferences(figurePars)
