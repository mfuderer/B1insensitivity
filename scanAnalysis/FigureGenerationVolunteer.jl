## This script generates figures and variance analysis for the volunteer data

include("FigureGenerationSetup.jl")
include("figuresB1phantom.jl")
include("scanAnalysisB1Brain.jl")
using NIfTI

meanRoiValuesCollection = []
devRatiosCollection=[]
analyzedType = 2      # T2 map
for volunteer in 1:2

    figurePars["dispMax"]       = [2.0,0.18,1.0]
    PyPlot.rc("font", family="serif")
    PyPlot.rc("font", size=14)
    percDev = zeros(2,2,4) # smoothState, condensedRecon (PreMis vs Non), sequence

    figurePars["filename"] = "Volunteer$volunteer"
    # figuresReadData!(figurePars)
    figuresReadCompactData!(figurePars)
    mapSet   = figurePars["mapSet"]   

    figureAllIm(figurePars)
    figureDifferences(figurePars)

    tissueName                = ["GM", "WM"]
    threshold_settings        = [0.5, 0.7] # apparently it needs to be much lower for GM than for WM 

    for tissueLabel in 1:2
        threshold                  = threshold_settings[tissueLabel]

        devRatios = zeros(length(figurePars["recDescription"]),length(figurePars["seqDescription"]))
        meanRoiValues= zeros(length(figurePars["recDescription"]),length(figurePars["seqDescription"]))
        preMisRatios = zeros(2,                                   length(figurePars["seqDescription"])) # over two smooth states

        # checking results with common ROI     
        meanRoi, stdRoi = BrainSegment(1, 1, figurePars, tissueLabel,threshold) # taking the pre-calibrated no-phase roi for all   
        wmroi         = figurePars["wmroi"]           # FH, AP, LR  
        sliceRange    = figurePars["sliceRange"]   

        @printf("Tissue %1d \n", tissueLabel)
        for case in 1:nSeqs
            seqDesc = figurePars["seqDescription"][case] 
            for r in eachindex(figurePars["recDescription"])
                recDesc = figurePars["recDescription"][r]
                t2map   = mapSet[r][case,sliceRange,analyzedType,:,:] # FH, AP, LR
                meanRoi = mean(t2map[wmroi])
                stdRoi  = std(t2map[wmroi]);
                @printf("For %s recon of %s sequence, mean T%1d is %.1f ms, dev is %.1f ms, ratio %.1f %% \n", 
                                recDesc, seqDesc, analyzedType, 1000*meanRoi, 1000*stdRoi, 100*stdRoi/meanRoi)
                devRatios[r,case] = 100*stdRoi/meanRoi
                meanRoiValues[r,case] = meanRoi
            end
            preMisRatios[1,case] = 100*(meanRoiValues[2,case] - meanRoiValues[1,case]) / meanRoiValues[1,case]
        end
        percDev[1,1,:] = 0.5 .* devRatios[1,:] .+ 0.5 .* devRatios[2,:] # collapsed the 'Pre' and 'Mis' deviations
        percDev[1,2,:] =        devRatios[3,:]                          # assumed to refer to the 'Non' deviations

        wmroi = figurePars["wmroi"]                     # FH, AP, LR
        B1map = figurePars["dreamSet"][1][2,:,:,:]      # AP, LR, FH
        # figure(); imshow(B1map[:,:,10], cmap="gray", vmin=0.8, vmax=1.2)
        @show std(permutedims(B1map,(3,1,2))[wmroi])    # FH, AP, LR
        B1mapRoi = zeros(size(B1map))

        b1map = figurePars["dreamSet"][1][2,:,:,:]  # AP, LR, FH
        # figure(); imshow(b1map[:,:,10], cmap="gray", vmin=0.8, vmax=1.2)
        wmroiTrans    = [CartesianIndex(p[2],p[3],p[1]) for p in wmroi]      # AP, LR, FH
        b1mapRoi = ones(size(b1map))
        for p in wmroiTrans
            b1mapRoi[p] = b1map[p]
        end 
        #figure(); imshow(b1mapRoi[:,:,10], cmap="gray", vmin=0.8, vmax=1.2)
        meanROIb1 = mean(b1map[wmroiTrans])
        stdROIb1  = std(b1map[wmroiTrans])
        @printf("mean B1 is %.1f, dev is %.1f, ratio %.1f %% \n", meanROIb1, stdROIb1, 100*stdROIb1/meanROIb1)

        # Smoothed over ROI
        for case in 1:nSeqs
            seqDesc = figurePars["seqDescription"][case] 
            for r in eachindex(figurePars["recDescription"])
                recDesc = figurePars["recDescription"][r]
                smoothRoi = SmoothOverMask(case, r, figurePars, analyzedType)
                meanRoi = mean(smoothRoi[wmroi])
                stdRoi  = std(smoothRoi[wmroi]);
                @printf("For %s recon of %s sequence, mean T%1d is %.1f ms, dev is %.1f ms, ratio %.1f %% \n", 
                                recDesc, seqDesc, analyzedType, 1000*meanRoi, 1000*stdRoi, 100*stdRoi/meanRoi)
                devRatios[r,case] = 100*stdRoi/meanRoi
                meanRoiValues[r,case] = meanRoi
                # if (case,r)==(3,3)
                #     figure(); imshow(smoothRoi[10,:,:])
                # end
            end
            preMisRatios[2,case] = 100*(meanRoiValues[2,case] - meanRoiValues[1,case]) / meanRoiValues[1,case]
        end
        percDev[2,1,:] = 0.5 .* devRatios[1,:] .+ 0.5 .* devRatios[2,:] # collapsed the 'Pre' and 'Mis' deviations
        percDev[2,2,:] =        devRatios[3,:]                          # assumed to refer to the 'Non' deviations

        # plotting of a smoothed deviation graph
        figure()
        ylim(0,30)
        xticks([0,1,2],figurePars["recDescription"])
        for case in 1:4
            plot(devRatios[:,case], label=figurePars["seqDescription"][case])
        end
        ylabel("relative standard dev [%]")
        legend()
        title("Tissue $(tissueName[tissueLabel]), volunteer $volunteer, smoothed")

        r=1; case=1   
        recDesc = figurePars["recDescription"][r]
        smoothRoi = SmoothOverMask(case, r, figurePars, 2)
        # figure(); imshow(smoothRoi[10,:,:])

        # does the smoothing affect the estimate of B1-variation over ROI?
        B1mapTrans = permutedims(B1map,(3,1,2))
        @show std(B1mapTrans[wmroi])    
        justRoi   = zeros(size(B1mapTrans))
        mapRoi    = zeros(size(B1mapTrans)) # ; mapRoi .= mean(map[wmroi])
        smoothRoi = copy(mapRoi)
        for p in wmroi
            justRoi[p] = 1.0
            mapRoi[p]  = B1mapTrans[p]  
        end 

        kernel = Kernel.gaussian((1,1,1))
        convolvedMapRoi  = imfilter(mapRoi, kernel)
        convolvedJustRoi = imfilter(justRoi, kernel)
        for p in wmroi
            smoothRoi[p] = convolvedMapRoi[p] / convolvedJustRoi[p]
        end 
        @show std(smoothRoi[wmroi])    

        varCorSmooth = 0.975  # approximate reduction of noise-variance by smoothing-denoising
        varCorB1     = 1.0 #0.75   # guesstimate of B1-variance-reduction by B1-correction 
        devB1        = 9.9  # measured B1-deviation over ROI
        #b1factors = [1.9, 1.85, 1.2, -0.1]
        b1factors = preMisRatios[2,:]./5.0

        noiseVarianceEst = (percDev[1,:,:].^2 .- percDev[2,:,:].^2)./varCorSmooth
        b1NonVarViaFactor = (devB1.*b1factors[:]).^2
        b1NonVarViaVarB1  = (percDev[:,2,:].^2 .- percDev[:,1,:].^2)./varCorB1
        totVarNon        = percDev[:,2,:].^2

        println()
        @printf("Percent difference in value between Mis and Pre, for unsmoothed and smoothed: \n")
        for i in 1:2
            for value in preMisRatios[i,:]
                @printf("%6.2f ", value)
            end
            println()
        end

        analysisTypes = [noiseVarianceEst, b1NonVarViaFactor, b1NonVarViaVarB1, totVarNon]
        analysisNames = ["noiseVarianceEst", "b1NonVarViaFactor", "b1NonVarViaVarB1", "totVarNon"]
        println()
        for (i,thisArray) in enumerate([noiseVarianceEst, b1NonVarViaFactor, b1NonVarViaVarB1, totVarNon])
            @printf("%s: \n",analysisNames[i])
            for row in axes(thisArray,1)
                for value in thisArray[row,:]
                    @printf("%6.2f ", value)
                end
                println()
            end
            println()
        end

        @printf("Condensed expected noise variance, expected B1 variance, measured total variance: \n")
        noiseVarianceCondensed = mean(noiseVarianceEst, dims=1) 
        b1VarianceCondensed    = [(0.1*b1NonVarViaFactor[seq]+b1NonVarViaVarB1[1,seq]+b1NonVarViaVarB1[2,seq])/2.1 for seq in axes(percDev,3)]
        for thisArray in [noiseVarianceCondensed, b1VarianceCondensed, totVarNon[1,:]]
            for value in thisArray
                @printf("%6.2f ", value)
            end
            println()
        end
        remainingVariance = totVarNon[1,:] .- noiseVarianceCondensed[1,:] .- b1VarianceCondensed
        @printf("Yet unexplained variance: \n")
        for value in remainingVariance
            @printf("%6.2f ", value)
        end
        println()

        push!(meanRoiValuesCollection, meanRoiValues)
        push!(devRatiosCollection, devRatios)
    end
end

devRatiosMean = zeros(size(devRatiosCollection[1]))
for i in 1:4 # over volunteers*tissue_types (GM, WM)
    devRatiosMean = devRatiosMean .+ devRatiosCollection[i].*(1/4.0)
end
meanRoiGM = (meanRoiValuesCollection[1] .+ meanRoiValuesCollection[3]) ./ 2.0
meanRoiWM = (meanRoiValuesCollection[2] .+ meanRoiValuesCollection[4]) ./ 2.0

fig,ax = subplots(1,1,figsize=(6,6))
ax.set_ylim(0,30)
xticks([0,1,2],figurePars["recDescription"])
for case in 1:4
    ax.plot(devRatiosMean[:,case], label=figurePars["seqDescription"][case], marker=".",markersize=20)
    # labelTissue = (case==1) ? "Grey matter" : ""
    # ax[2].plot(meanRoiGM[:,case], label=labelTissue)
    # labelTissue = (case==1) ? "White matter" : ""
    # ax[2].plot(meanRoiWM[:,case], "--", label=labelTissue)
end
ax.set_ylabel("relative standard dev [%]")
ax.legend()

@show meanRoiGM[1,:]
@show meanRoiWM[1,:]
# Mis-Pre 
@show (meanRoiGM[2,:].-meanRoiGM[1,:])./meanRoiGM[1,:]./0.05
@show (meanRoiWM[2,:].-meanRoiWM[1,:])./meanRoiWM[1,:]./0.05
