
function map_B1_sensitivities(fn_base,nShownRealizations)
    recon_options["handleB1"] = "sensitivity" 
    recon_options["B1metric"] = "derivative_at_1"         
    saved_H = Dict()

    defaultT1T2mix = copy(recon_options["T1T2set"])
    appliedT1T2mix = [(1.3,0.2)] 
    t1c = 1.0; t2c = 0.08;    

    b1_range  = 0.72:0.02:1.28
    lt1_range = -1.4:0.1:1.4
    lt2_range = -1.4:0.1:1.4
    b1sensitivities = zeros(length(lt1_range),length(lt2_range),2)
    noisesMap       = zeros(length(lt1_range),length(lt2_range),2)

    (fig,ax)=(subplots(Int(ceil(nShownRealizations/3)),3,figsize=(9,3)))
    for i in 1:nShownRealizations
        fn = "$fn_base($i)"
        recon_options["rfFile"]  = fn
        RFdeg = RF_from_file(recon_options["nTR"], recon_options["nky"])
        BLAKJac.PlotAmplAndPhaseDeriv(RFdeg, ax[i])

        noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
        @show fn, noises
    end

    pcm=0
    (fig,ax) = subplots(nShownRealizations,5,figsize=(18,18))

    for rrr in 1:nShownRealizations
        fn = "$fn_base($rrr)"
        recon_options["rfFile"]  = fn
        RFdeg = RF_from_file(recon_options["nTR"], recon_options["nky"])
        BLAKJac.PlotAmplAndPhaseDeriv(RFdeg, ax[rrr,1])

        for (i1,lt1) in enumerate(lt1_range)
            for (i2,lt2) in enumerate(lt2_range)
                recon_options["T1T2set"] =[(t1c*exp(lt1),t2c*exp(lt2))] 
                noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H)
                b1sensitivities[i1,i2,1:2] = b1f[2:3]
                noisesMap[i1,i2,1:2]       = noises[2:3]
            end
        end

        for m in 1:2 
            pcm=ax[rrr,m+1].imshow(b1sensitivities[:,:,m], cmap="RdBu",vmin=-1.5,vmax=1.5)
            ax[1,m+1].set_title("T$m-sensitivity to B1+")
        end
        colorbar(pcm, ax=ax[rrr,3])

        for m in 1:2 
            pcm=ax[rrr,m+3].imshow(noisesMap[:,:,m],vmin=2.0,vmax=6.0)
            ax[1,m+3].set_title("T$m noise factor")
        end
        colorbar(pcm, ax=ax[rrr,5])
        for loc in 2:5
            ax[rrr,loc].set_xticks([])
            ax[rrr,loc].set_yticks([])
            ax[nShownRealizations,loc].set_xticks([0,length(lt1_range)-1],[@sprintf("%.0f ms", 1000*t2c*exp(lt2_range[1])), @sprintf("%.0f ms", 1000*t2c*exp(lt2_range[end]))])
            ax[nShownRealizations,loc].set_xlabel("T2")
        end
        ax[rrr,2].set_yticks([0,length(lt1_range)-1],[@sprintf("%.1f s", t1c*exp(lt1_range[1])), @sprintf("%.1f s", t1c*exp(lt1_range[end]))])
        ax[rrr,2].set_ylabel("T1")

        t1SetRescaled = [(log(t1/t1c)+1.4  )*10+1 for (t1,t2) in defaultT1T2mix]
        t2SetRescaled = [(log(t2/t2c)+1.4  )*10+1 for (t1,t2) in defaultT1T2mix]
        ax[1,2].scatter(t2SetRescaled,t1SetRescaled,color="orange")
    end

    # Plot over (B1,T2)
    len1 = length(b1_range)
    (fig,ax) = subplots(nShownRealizations,5,figsize=(18,18))
    for rrr in 1:nShownRealizations
        fn = "$fn_base($rrr)"
        recon_options["rfFile"]  = fn
        RFdeg = RF_from_file(recon_options["nTR"], recon_options["nky"])
        ax[rrr,1].set_ylim(-10,100)
        BLAKJac.PlotAmplAndPhaseDeriv(RFdeg, ax[rrr,1])

        (t1app,dummy) = appliedT1T2mix[1]
        for (i1,b1) in enumerate(b1_range)
            for (i2,lt2) in enumerate(lt2_range)
                recon_options["T1T2set"] =[(t1app,t2c*exp(lt2))] 
                noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg.*b1, trajectorySet, recon_options, saved_H)
                b1sensitivities[i1,i2,1:2] = b1f[2:3]
                noisesMap[i1,i2,1:2]       = noises[2:3]
            end
        end

        for m in 1:2 
            pcm=ax[rrr,m+1].imshow(b1sensitivities[:,:,m], cmap="RdBu",vmin=-1.5,vmax=1.5)
            ax[1,m+1].set_title("T$m-sensitivity to B1+")
        end
        colorbar(pcm, ax=ax[rrr,3])

        for m in 1:2 
            pcm=ax[rrr,m+3].imshow(noisesMap[:,:,m],vmin=2.0,vmax=7.0)
            ax[1,m+3].set_title("T$m noise factor")
        end
        colorbar(pcm, ax=ax[rrr,5])
        for loc in 2:5
            ax[rrr,loc].set_xticks([])
            ax[rrr,loc].set_yticks([])
            ax[end,loc].set_xticks([0,length(lt1_range)-1],[@sprintf("%.0f ms", 1000*t2c*exp(lt2_range[1])), @sprintf("%.0f ms", 1000*t2c*exp(lt2_range[end]))])
            ax[end,loc].set_xlabel("T2")
        end
        ax[rrr,2].set_yticks([0,len1-1],[@sprintf("%.1f", b1_range[1]), @sprintf("%.1f", b1_range[end])])
        ax[rrr,2].text(-5,18,"B1")

        b1centers = [(0+1.4)*10+1 for (t1,t2) in defaultT1T2mix]
        t2SetRescaled = [(log(t2/t2c)+1.4  )*10+1 for (t1,t2) in defaultT1T2mix]
        ax[1,2].scatter(t2SetRescaled,b1centers,color="orange")
    end
    recon_options["T1T2set"] = copy(defaultT1T2mix)
end
