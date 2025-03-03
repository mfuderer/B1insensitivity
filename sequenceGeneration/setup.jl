## Define packages
# using LinearAlgebra # package for linear algebra operations (dot, )
# using DelimitedFiles # For input of .txt files containing RF-patterns
# using FFTW
using FileIO
using ComputationalResources
# using CUDA

# # Before adding PyPlot, do
# #     ENV["PYTHON"] = "" in the Julia REPL.
# #     Add PyPlot in Pkg mode
# #     Pkg.build("PyCall")
# #     In vscode, deactivate settings > Extensions > Julia > Use Plot Plane
# #     Restart Julia
# using PyPlot # same plotting backend as commonly used in python
# using PyCall

cpu  = ComputationalResources.CPU1()
# using BlochSimulators
# using Colors
using Statistics
using BLAKJac
# using MRSTATToolbox
using Optim
using Printf
using PyPlot
# include("recon_options.jl")

# include("../numerical_phantom/RF_Shapes.jl")
# include("../numerical_phantom/k_shapes.jl")
# # src/MRSTAT/src/Utils/RF_Shapes.jl
# if (!@isdefined(gelphantoms)) include("gelphantoms.jl") end
# if (!@isdefined(tissues)) include("tissues.jl") end
# @pyimport matplotlib.animation as anim




function RF_from_file(nTR, ny)
    RF=complex.(zeros(nTR))
    #@assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    fn = recon_options["rfFile"]
    fPath = string(recon_options["rfFolder"],fn,".jld2")
    @show fPath
    vars = FileIO.load(fPath)
    RFdeg = vars["RFdeg"]
    RF_train = vec(complex.(RFdeg))
    return RF_train
end
