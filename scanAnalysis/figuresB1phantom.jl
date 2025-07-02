
sx = 224; sy = 224; nSweeps = 6
nZoom = 1


function ShowTubes(ρ, tubes)
    figure()
    subplot(121)
        # overlay of ROIs on proton density map
        ρ_max = maximum(ρ);
        for i in eachindex(tubes)
            ρ[tubes[i]] .= 2*ρ_max;
        end
        imshow(ρ)
        title("Mask per tube")

    subplot(122)
        # numbered
        x = zeros(size(ρ)...)
        for (i,r) in enumerate(tubes)
            x[r] .= i
        end
        imshow(x)
        colorbar()
end

function figuresReadCompactData!(figurePars)
    filename = figurePars["filename"]
    folder  = figurePars["dataFolder"]
    #fn_full = folder*filename*".jld2" 
    fn_full = joinpath(folder,(filename*".jld2")) 
    @load fn_full mapSet dreamSet RFdeg RFphaseDeg
    figurePars["mapSet"]          = mapSet
    figurePars["dreamSet"]        = dreamSet
    figurePars["RFdeg"]           = RFdeg
    figurePars["RFphaseDeg"]      = RFphaseDeg
end


function phantomVialMeans!(figurePars)
    mapSet        = figurePars["mapSet"]
    sliceRange    = figurePars["sliceRange"]         # The range of slices that are available 
    nRecons       = length(figurePars["recDescription"])  
    nSeqs         = length(figurePars["seqDescription"])     

    # Measured by Oscar using MIR-MSE measurements
    goldstandard_3T = (
            #         1    2    3    4    5    6    7    8    9    10   11   12    13    14    15    16    17    18
        T1 = Float64[216, 323, 316, 500, 488, 480, 646, 629, 814, 801, 965, 1645, 1047, 1120, 1266, 1431, 1394, 1555],
        T2 = Float64[ 44,  61,  96,  45,  78, 128,  78, 109,  96, 127, 108,  303,  179,  132,  167,  161,  145,  140]
    )

    symbolicPos = [(0,1), (0,2), (1,0), (1,1), (1,2), (1,3), (2,0), (2,1), (2,2), (2,3), (3,1), (3,2), (2.6,2.6), (1.5,1.5)];
    tubeLabel   = [  14,     9,    11,     1,     2,     4,    18,   10,    13,     16,      5,     7,          0,        0]
    origin = (67,33.5)
    distance = 43.8
    angulation = 0.24
    radius = 7 
    refSlice = 10
    usableSliceRange = 8:14 
    centerPos = [origin.+(distance*x*cos(angulation)-distance*y*sin(angulation), distance*x*sin(angulation)+distance*y*cos(angulation)) for (x,y) in symbolicPos]
    tube_centers = [CartesianIndex(round(Int64,x),round(Int64,y)) for (x,y) in centerPos]

    refRec = (length(sliceRange) > 1) ? refSlice : 1
    ρ = mapSet[1][1,refRec,3,:,:]

    d(x::CartesianIndex, y::CartesianIndex) = √( (x[1]-y[1])^2 + (x[2]-y[2])^2 )
    allidx = CartesianIndices(ρ)

    tubes = [allidx[ d.((center,), allidx) .< radius] for center ∈ tube_centers];
    ShowTubes(ρ, tubes)

    # Take mean values over ROIs
    indexRange= (length(sliceRange) > 1) ? sliceRange : (1:1)
    meanMeans = zeros(nRecons, nSeqs,2,length(tubes))
    devDevs   = zeros(nRecons, nSeqs,2,length(tubes)) # disused
    for recType in 1:nRecons
        for case in 1:nSeqs
            for m in 1:2
                meanIm = mapSet[recType][case,indexRange,m,:,:] 
                for (i,tube) in enumerate(tubes)
                    meanMeans[recType,case,m,i] = mean(meanIm[usableSliceRange,tube])
                end        
            end
        end
    end

    figurePars["meanMeans"]       = meanMeans
    figurePars["goldstandard_3T"] = goldstandard_3T
    figurePars["tubes"]           = tubes
    figurePars["tubeLabel"]       = tubeLabel
end

function figuresPhantomB1map(figurePars)
    nSeqs         = length(figurePars["seqDescription"])     
    dreamSet      = figurePars["dreamSet"]          
    sliceRange    = figurePars["sliceRange"]
    midSlice      = Int64(round((sliceRange[1]+sliceRange[end])/2))
    recPick       = midSlice

    fig_B1s, ax_B1s  = subplots(1,3,figsize=(10,3))
    fig_B1s.subplots_adjust(wspace=-0.15,hspace=0.0)
    fig_B1s.suptitle("B1 maps for Pre-calibrated, Mis-calibrated and Non-calibrated reconstructions", fontsize=12)

    # plot Dream maps 
    B1map = dreamSet[1][1,:,:,midSlice]
    maskOnes = ones(size(B1map)) .* (B1map .> 0.01)

    vminB1 = 0.75; vmaxB1=1.25
    pcm = ax_B1s[1].imshow(B1map,  vmin=vminB1,vmax=vmaxB1)
    ax_B1s[2].imshow(1.05 .* B1map,vmin=vminB1,vmax=vmaxB1)
    ax_B1s[3].imshow(maskOnes,     vmin=vminB1,vmax=vmaxB1)

    cbar_ax = fig_B1s.add_axes([0.9, 0.15, 0.05, 0.7])
    clb = fig_B1s.colorbar(pcm, cax=cbar_ax)
    clb.ax.tick_params(labelsize=14) 
    for i in 1:3
        ax_B1s[i].set_xticks([]); 
        ax_B1s[i].set_yticks([]);
        ax_B1s[i].set_yticklabels([]);  
    end
end

function figureSequence2by2(figurePars)
    RFdeg         = figurePars["RFdeg"]           
    RFphaseDeg    = figurePars["RFphaseDeg"] 

    fig_seqs,ax_seqs = subplots(2,2,figsize=(8,8))
    fig_seqs.subplots_adjust(wspace=0.0,hspace=0.0,bottom=0.2)

    for r in 1:2
        for c in 1:2
            case=2*(r-1)+c
            ax_seqs[r,c].set_ylim(-10.0,90.0)
            ax_seqs[r,c].plot(RFdeg[case,:],label="amplitude")
            anglesdd = zeros(size(RFdeg)[2])
            for i in 1:size(RFdeg)[2]-2
                anglesdd[i] = ((-RFphaseDeg[case,i]+2*RFphaseDeg[case,i+1]-RFphaseDeg[case,i+2]+720.0+270.0)) % 180.0 - 90.0
            end
            ax_seqs[r,c].plot(anglesdd,label="ϕ''")  
            ax_seqs[r,c].set_xlabel("pulse number")

        end
        ax_seqs[r,1].set_yticks([])
        ax_seqs[r,2].yaxis.set_label_position("right")
        ax_seqs[r,2].yaxis.tick_right()
        ax_seqs[r,2].set_ylabel("amplitude [deg] \n or ϕ'' [deg/TR²]")
    end
    for c in 1:2
        ax_seqs[1,c].set_xticks([])
    end

    rowLabels = ["Amplitude-only", "Amplitude+phase"]
    colLabels = ["noise-optimized","B1-insensitivity-optimized"]
    for i in 1:2
        fig_seqs.text(0.08, 0.75 - (i - 1) * 0.4, rowLabels[i], va="center", ha="center", rotation="vertical", fontsize=14)
    end
    for i in 1:2
        fig_seqs.text(0.25 + (i - 1) * 0.5, 0.90, colLabels[i], va="center", ha="center", fontsize=14)
    end
end

function figuresPhantomBias(figurePars, diffusion_corrected=false, type=2) #T2 by default
    mapSet = figurePars["mapSet"]  
    meanMeans   = figurePars["meanMeans"]       
    goldstandard_3T=figurePars["goldstandard_3T"]
    tubes       = figurePars["tubes"]
    tubeLabel   = figurePars["tubeLabel"] 
    nSeqs         = length(figurePars["seqDescription"])     
    gels          = figurePars["gels"]               # The set of vials to be used
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    seqDescription= figurePars["seqDescription"]     # vector of description of sequence types
    comparisonSet = figurePars["comparisonSet"]      # vector of two numbers, referring to indices of folder names; typically [1,2]
    diffuCorrFile = figurePars["diffuCorrFile"]      # location/name of file containing diffusion correction data
    dctxt = ""
    if diffusion_corrected
        @load diffuCorrFile mapCollection
        dctxt = " (diffusion-corrected)"
    end
    PyPlot.rc("font", size=9)

    lines_color_cycle = [p["color"] for p in plt.rcParams["axes.prop_cycle"]]

    # Prepare Bias list and relative bias 
    gtT = goldstandard_3T[type]
    xxx = ([gtT[tubeLabel[gel]] for gel in gels])
    biasTn    = zeros(length(recDescription),nSeqs,length(gels))
    biasTrel  = zeros(length(recDescription),nSeqs,length(gels))
    for recType in eachindex(recDescription)
        for case in 1:nSeqs
            biasTn[recType,case,:] = 1000.0 .* meanMeans[recType,case,type,gels].-xxx
            if diffusion_corrected
                diffCorr = [mapCollection[case][type,tubeLabel[gel]] for gel in gels]
                biasTn[recType,case,:] = biasTn[recType,case,:] .+ 1000.0 .* diffCorr
            end 
            biasTrel[recType,case,:] = biasTn[recType,case,:]./xxx
        end
    end

    fig,ax = subplots(1,3,figsize=(12,4))
    markers = ["o","x","+"]
    reconNames = recDescription

    m = type
    localSeqDescription = seqDescription
    localSeqDescription = ["A-only noise-opt", "A-only B1-opt", "A+P noise-opt", "A+P B1-opt"]
    for case in 1:nSeqs
        label = localSeqDescription[case] 
        color = lines_color_cycle[case]
        yyy = 1000.0 .* (meanMeans[comparisonSet[2],case,m,gels].-meanMeans[comparisonSet[1],case,m,gels])
        # make it relative: 
        yyyrel = 100.0 .* yyy./xxx
        @show label, mean(yyyrel)

        ax[1].scatter(xxx,yyyrel,label=label, color=color)
        # slope = mean(yyy) / mean(xxx)
        # @show label, slope
        # xl = minimum(xxx); xh = maximum(xxx)
        # ax[1].plot([xl,xh],slope.*[xl,xh], color=color)        
    end
    #ax[1].set_ylim(-5.0,23.0)
    ax[1].set_xlabel("Gold standard T$type [ms]")
    ytext = "T$type estimation difference [%] \n due to 5% offset in B1"*dctxt
    ax[1].set_ylabel(ytext)
    ax[1].set_title("Miscalibration effect on T$type [%]", fontsize=9)
    ax[1].legend(fontsize=9)

    for recType in eachindex(recDescription)
        for case in 1:nSeqs
            label = ""
            label = case==1 ? recDescription[recType] : "" 
            color = lines_color_cycle[case]
            marker = markers[recType]
            yyy = 100.0 .* biasTn[recType,case,:] ./ xxx
            ax[2].scatter(xxx,yyy,label=label, color=color, marker=marker)
        end
    end
    ax[2].set_xlabel("Gold standard T$type [ms]")
    ax[2].set_title("T$type bias [%]", fontsize=9)
    ax[2].legend(fontsize=9)

    # Alessandro's suggested deviation-bar graphs
    spacer = 4 + 2
    for recType in eachindex(recDescription)
        for case in 1:nSeqs
            pos = spacer*(recType-1) + (case-1)
            ave = 100.0 * mean(biasTrel[recType,case,:])
            dev = 100.0 * std(biasTrel[recType,case,:])
            color = lines_color_cycle[case]
            ax[3].bar(pos,ave,color=color)
            ax[3].errorbar(pos,ave,dev,linewidth=2.0,capsize=8.0,color="black")
        end
    end
    ax[3].set_xticks([]); 
    ax[3].plot([-0.8,(3-1)*spacer+4.0],[0.0,0.0],color="black")
    ax[3].yaxis.set_label_position("right")
    ax[3].yaxis.tick_right()
    ax[3].set_ylabel("T$type bias in %", fontsize=9)  
    for i in 1:3 
        ax[i].set_ylim(-22.0,53.0)
    end
    yPositions = [24,32,39] # [9,12,9]
    for recType in eachindex(recDescription)
        ax[3].text(spacer*(recType-1)-1.4,yPositions[recType],recDescription[recType], fontsize=8)
    end
end

function figureAllIm(figurePars)
    nSeqs         = length(figurePars["seqDescription"])     
    nRecons       = length(figurePars["recDescription"])       
    mapSet        = figurePars["mapSet"]     
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    seqDescription= figurePars["seqDescription"]     # vector of description of sequence types
    dispMax       = figurePars["dispMax"]
    sliceRange    = figurePars["sliceRange"]
    midSlice      = Int64(round((sliceRange[1]+sliceRange[end])/2))
    recPick       = length(sliceRange) > 1 ? midSlice : 1

    # (all of T1 maps and T2 maps for scan*recon combinations)

    wRat = ones(nSeqs)
    ddd=Dict("width_ratios" => wRat)

    for m in 1:2
        fig_maps,ax_maps = subplots(nRecons,nSeqs,figsize=(13.5,9.5),gridspec_kw=ddd)
        fig_maps.subplots_adjust(wspace=-0.05,hspace=0.0,right=0.8, top=0.9, bottom=0.0, left=0.0)
        for case in 1:nSeqs
            for r in 1:nRecons
                row = r
                # case zoom [dummy] type iter x y
                meanIm = mapSet[r][case,recPick,m,:,:] 

                rgb_vec, imClip = relaxationColorMap("T$m", meanIm, 0.0, dispMax[m])  # call to resource, generating a colormap 
                cmap = PyPlot.ColorMap("relaxationColor", rgb_vec, length(rgb_vec), 1.0) 
                pcm = ax_maps[row,case].imshow(imClip, vmax=dispMax[m], cmap=cmap);

                if (case==1 && r==1)
                    fig_maps.subplots_adjust(right=0.8)
                    cbar_ax = fig_maps.add_axes([0.85, 0.15, 0.05, 0.7])
                    clb = fig_maps.colorbar(pcm, cax=cbar_ax, ticks=[0, dispMax[m]])
                    clb.ax.tick_params(labelsize=20) 
                    dispMaxVal = 1000*dispMax[m]
                    clb.ax.set_yticklabels(["0ms", "$dispMaxVal ms"], fontsize=14)
                end

                ax_maps[row,case].set_xticks([]); 
                ax_maps[row,case].set_yticks([]);
                ax_maps[row,case].set_yticklabels([]);  
            end
            ax_maps[1,case].set_title(seqDescription[case])
        end 
        for r in 1:nRecons
            ax_maps[r,1].text(0.02,0.3, recDescription[r], transform=ax_maps[r,1].transAxes, rotation="vertical", color="white")
        end
    end
end   


function figureDifferences(figurePars)
    nSeqs         = length(figurePars["seqDescription"])     
    mapSet        = figurePars["mapSet"]     
    seqDescription= figurePars["seqDescription"]     # vector of description of sequence types
    recDescription= figurePars["recDescription"]     # vector of description of recon types
    dispMax       = figurePars["dispMax"]
    sliceRange    = figurePars["sliceRange"]
    midSlice      = Int64(round((sliceRange[1]+sliceRange[end])/2))
    recPick       = length(sliceRange) > 1 ? midSlice : 1

    for m in 1:2
        # For the maps, differences between the first reconstruction mode and the rest  
        if length(recDescription) > 1   # only meaningful if there is a 'rest'
            wRat = ones(nSeqs)
            ddd=Dict("width_ratios" => wRat)

            fig_rows = min(2,length(recDescription)-1)
            fig_maps,ax_maps = subplots(2,nSeqs,figsize=(17,8),gridspec_kw=ddd)  
            fig_maps.subplots_adjust(wspace=-0.05,hspace=-0.1,right=0.9, top=0.9, bottom=0.0, left=0.0)

            for case in 1:nSeqs
                meanImRef = mapSet[1][case,recPick,m,:,:] 
                for r in 2:length(recDescription)
                    row = r-1
                    # case zoom [dummy] type iter x y
                    meanIm = mapSet[r][case,recPick,m,:,:] 

                    vmax = [0.12,0.012][m]; vmin = -vmax;

                    pcm = ax_maps[row,case].imshow(meanIm .- meanImRef, vmin=vmin, vmax=vmax, cmap="RdBu");

                    if (case==1 && r==2)
                        fig_maps.subplots_adjust(right=0.8)
                        cbar_ax = fig_maps.add_axes([0.85, 0.15, 0.05, 0.7])
                        clb = fig_maps.colorbar(pcm, cax=cbar_ax, ticks = [vmin, 0, vmax])
                        clb.ax.tick_params(labelsize=20) 
                        clb.ax.set_yticklabels(["$(1000*vmin) ms", "0 ms", "$(1000*vmax) ms"], fontsize=14)
                    end

                    ax_maps[row,case].set_xticks([]); 
                    ax_maps[row,case].set_yticks([]);
                    ax_maps[row,case].set_yticklabels([]);  
                end
                ax_maps[1,case].set_title(seqDescription[case], fontsize=16)    
            end
            for r in 2:length(recDescription)
                vtext = recDescription[r]*" - "*recDescription[1] 
                ax_maps[r-1,1].text(0.02,0.15, vtext, transform = ax_maps[r-1,1].transAxes, rotation="vertical", color="black")
            end
        end
    end
end
