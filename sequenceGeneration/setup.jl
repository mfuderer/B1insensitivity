## Define packages
using FileIO
using ComputationalResources
cpu  = ComputationalResources.CPU1()
using Statistics
using BLAKJac
using Optim
using Printf
using PyPlot

function RF_from_file(nTR, ny)
    RF=complex.(zeros(nTR))
    fn = recon_options["rfFile"]
    fPath = string(recon_options["rfFolder"],fn,".jld2")
    @show fPath
    vars = FileIO.load(fPath)
    RFdeg = vars["RFdeg"]
    RF_train = vec(complex.(RFdeg))
    return RF_train
end
