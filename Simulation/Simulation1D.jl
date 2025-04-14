# Script for simulation of the effect of B1 on T1 and T2 estimation 
# 

# This function generates a linearly varying profile over the central portion of a 1-dimensional range 
# (the outer two quarters of the range being filled in such as to avoid discontinuities)
function sawtooth(segs, center, spread, useLog=true)
    saw = zeros(segs)
    qs = segs÷4
    saw[1:qs]      = spread .* collect(range(start=0, stop=-1, length=qs))
    saw[qs+1:3*qs] = spread .* collect(range(start=-1, stop=1, length=2*qs))
    saw[3*qs+1:4*qs]=spread .* collect(range(start=1, stop=0, length=qs))
    if useLog
        saw = center .* exp.(saw)
    else
        saw = center .* max.((saw.+1.0),0.0)
    end
    return saw
end

include("setup.jl")
include("load_data_phantom.jl")

# PyPlot.rc("font", family="serif")
# PyPlot.rc("font", size=14)
PythonPlot.rc("font", family="serif")
PythonPlot.rc("font", size=14)

recon_options = Dict()

recon_options["recon_folder"] = "tmp";
recon_options["recon_subfolder"] = "tmp";
recon_options["numphantom"] = true;
recon_options["numphantom_rf_shape"] = "from_file" 
recon_options["rfFolder"]= pwd()*"/RFsequences/" # name of the folder to read the RF pattern from

recon_options["maxRho"]   = 1.0
recon_options["simulationT1center"]  = 0.7 # 1.3
recon_options["simulationT2center"]  = 0.07 # 0.2
recon_options["simulationB1center"]  = 1.0
recon_options["simulationVariable"] = "B1" # "TR/TD" # "measured from phantom"
recon_options["simulationSpread"] = 0.3 # used 0.8 for T1 and T2, 1.3 for δ # 0.3 for B1
recon_options["RhoVariation"]   = 0.0 

recon_options["numphantom_type"] = "line" 
recon_options["numphantom_size"] = (1,224)
nkx,nky = recon_options["numphantom_size"]
recon_options["numphantom_sequence"] = "Spoiled" 
recon_options["numphantom_trajectory"] = "Cartesian" 
recon_options["numphantom_noise"] = false
recon_options["simulation_parameters"] = T₁T₂B₁ρˣρʸ
recon_options["reconstruction_parameters"] = T₁T₂ρˣρʸ

recon_options["TR"]      = 0.01         # in seconds     

recon_options["slice"] = 1;
recon_options["coils"] = 1; # should be a range
recon_options["maxstate"] = 64; 
recon_options["lsqr_its"] = 10;
recon_options["trf_max_iter_steihaug"] = 30;
recon_options["trf_max_iter"] = 20;
recon_options["slice_discretization"] = 1
recon_options["slice_profile_correction"] = "shinnar_leroux" # small_tip_angle or shinnar_leroux
recon_options["slice_thickness_multiplier"] = 1 # compute slice profiles from 0 to multiplier * nominal slice thickness

recon_options["scaling_factor"] = 1.0;         

#adpt(x) = MRSTAT.Reconstruction.adapt_to(resource, x)
noplot(x; figtitle="") = println("no plotting")

spread = recon_options["simulationSpread"]
simT1 = recon_options["simulationT1center"] .* ones(nkx,nky)
simT2 = recon_options["simulationT2center"] .* ones(nkx,nky)
simB1 = recon_options["simulationB1center"] .* ones(nkx,nky)
for x in 1:nkx
    if recon_options["simulationVariable"]=="T1"
        simT1[x,:] = sawtooth(nky, recon_options["simulationT1center"], spread)
    elseif recon_options["simulationVariable"]=="T2"
        simT2[x,:] = sawtooth(nky, recon_options["simulationT2center"], spread)
    elseif recon_options["simulationVariable"]=="B1"
        simB1[x,:]  = sawtooth(nky, recon_options["simulationB1center"], spread, false)
    elseif recon_options["simulationVariable"]=="measured from phantom";
        ; # no overwrites
    else @assert false
    end
end

simRho = recon_options["maxRho"] .* ones(nkx,nky)
simMask = ones(nkx,nky)
simRho .*= simMask .|> Complex

B₁map = ones(nkx,nky)

description         = ["noise-optimized, no phase",    
                                        "B1-optimized,no phase",    
                                                          "noise-optimized, with phase",     
                                                                             "B1-optimized, with phase"]  
#cases               = ["20240702V(5)",   "20240703W(5)",  "20240701S(5)",    "20240703X(3)"]
cases               = ["NoPhase_noise(5)", "NoPhase_B1opt(5)", "Phase_noise(5)", "Phase_B1opt(3)"]

nR                  = [1,                1,               1,                  1]

sweeps              = 6
startstate          = -1
mbiasT1 = zeros(length(cases),nky)
mbiasT2 = zeros(length(cases),nky)
mT2     = zeros(length(cases),nky)

# lines_color_cycle = [p["color"] for p in plt.rcParams["axes.prop_cycle"]]
(figd,axd)=subplots(1,figsize=(8,5))

# MRSTATToolbox can use options loaded from file. 
#    Oscar set them to something that looked useful for pure simulation stuff but please double-check
# options = load_options_from_file(pwd()*"/mrstat_options_in_silico.toml")      

slopes = zeros(length(cases))
for (caseIndex,case) in enumerate(cases)
    for r in 1:nR[caseIndex]
        recon_options["rfFile"]  = nR[caseIndex]==1 ? case : case*"($r)"
        recon_options["nTR"]      = nky*sweeps
        recon_options["startstate"] = startstate

        sequence, coordinates, coilmaps, trajectory, mask, phantom = load_data_phantom(recon_options);

        # Rename sim{T1,T2,B1} to {T1,T2,B1} and simRho to {PDx,{Dy} to be able ot use @parameters from BlochSimulators
        T1, T2, B1 = simT1, simT2, simB1
        PDx, PDy = real(simRho), imag(simRho)
        parameters_with_B1 = @parameters T1 T2 B1 PDx PDy
        recon_options["simulation_parameters"] = T₁T₂B₁ρˣρʸ

        # Set all the things to single precision
        sequence = f32(sequence) |> gpu
        trajectory = f32(trajectory) |> gpu
        coilmaps = f32(reshape(only.(coilmaps), length(coilmaps),1)) |> gpu 
        trajectory = f32(trajectory) |> gpu 
        parameters_with_B1 = f32(parameters_with_B1) |> vec |> gpu
        coordinates = StructVector(f32(coordinates)) |> gpu |> vec
        transmit_field = simB1 |> f32 |> vec

        @show typeof(coilmaps), size(coilmaps), size(parameters_with_B1)

        # Rather than generating echos, applying phase encoding and then calling the Mv function, we can directly call simulate_signal from BlochSimulators. 
        # Running this on the CPU (multi-threaded) for now since it doesn't take long anyway. 
        # For 2D experiments it's a different story.

        raw_data = simulate_signal(CUDALibs(), sequence, parameters_with_B1, trajectory, coordinates, coilmaps)
       
        # All arguments should be of the same floating point precision, let's just do Float32 using BlochSimulators.f32's function        
        # The arrays also need to have a specific number of dimensions, add missing dimensions with size 1
        # For raw data: get data back to CPU, add "locations" and "repetitions" dimensions and convert to Float32
        # raw_data = reshape(raw_data, size(raw_data)..., 1, 1)
        
        # # For coordinates: Store as StructArray{<:Coordinates}, add "x", "z" and "locations" dimensions
        # coordinates = reshape(coordinates, 1, length(coordinates), 1, 1)
        
        # For coil sensitivity maps: Add "x", "z" and "locations" dimensions
        @assert size(coilmaps, 2) == 1
        # coil_sensitivities = reshape(coilmaps, 1, length(coilmaps), 1, 1, 1)
 
        # # For transmit field: Add "x", "z" and "locations" dimensions
        # transmit_field = reshape(complex.(parameters_with_B1.B₁), 1, length(coordinates), 1, 1)
 
        recon_options["simulation_parameters"] = T₁T₂ρˣρʸ
 
        # Call the MRSTAT reconstruction function
        # Note that this should be run on a machine with GPU 
        # ctx = MRSTATToolbox.MRSTAT.mrstat_recon_in_silico(sequence, trajectory, raw_data, coordinates, coil_sensitivities, transmit_field, options);
        ctx = MRSTAT.mrstat_recon(raw_data, sequence, coordinates, coilmaps, trajectory, transmit_field);

        final_optimpars = reshape(ctx.x[:,end], :, 4)

        # Extract the reconstructed T₁ and T₂ maps from the returned ctx object
        recon_res = (
            # ctx.reconstructed.T₁ has multiple named dimensions and resides on the 
            # with parent we get rid of the dimension names, with collect we get it back to 
            # the cpu and with vec we lose the dimensions
            T₁ = exp.(final_optimpars[:,1]),
            T₂ = exp.(final_optimpars[:,2]),
        )

        rT1 = reshape(recon_res.T₁, 1, nky)
        rT2 = reshape(recon_res.T₂, 1, nky)
        sT1 = reshape(simT1, 1, nky)
        sT2 = reshape(simT2, 1, nky)
        sB1 = reshape(simB1, 1, nky)

        mbiasT1[caseIndex,:] = mean(sT1.-rT1, dims=1) |> vec
        mbiasT2[caseIndex,:] = mean(sT2.-rT2, dims=1) |> vec
        mT2[caseIndex,:]     = mean(  sT2,    dims=1) |> vec

        sVar = copy(sT1)
        if recon_options["simulationVariable"]=="T1"; sVar = sT1;
        elseif recon_options["simulationVariable"]=="T2"; sVar = sT2;
        elseif recon_options["simulationVariable"]=="B1"; sVar = sB1;
        elseif recon_options["simulationVariable"]=="measured from phantom"
            from_measured =  measuredS[caseIndex]*sB1[1]/recon_options["TR"]*(sT2[1]*(sT2[1]+0.12))
            sVar .= from_measured;
        else @assert false
        end

        centerRange=nky÷2-50 : nky÷2+50
        color = lines_color_cycle[caseIndex]
        data = rT2[centerRange] .- sT2[centerRange] 
        xcoord = sVar[centerRange]  
        # data = rT1[centerRange] .- sT1[centerRange]   
        label = description[caseIndex]
        @show label, r, data[1]   
        # If r equals 1, add a label to the plot
        if r == 1
            axd.plot(sVar[centerRange], data, color=color, label=label)
        else
            axd.plot(sVar[centerRange], data, color=color)
        end
        if (nR[caseIndex]>1); axd.text(sVar[centerRange][1], data[1], string(r), color=color); end;
        slopes[caseIndex] = (data[51+10]-data[51-10])/(xcoord[51+10]-xcoord[51-10])
    end
end
axd.set_xlabel("value of $(recon_options["simulationVariable"])")
axd.set_ylabel("bias on T2 [s]")
axd.legend()
txtt1 = "T1=$(recon_options["simulationT1center"])s, "
txtt2 = "T2=$(recon_options["simulationT2center"])s, "
txtB1 = "B1=$(recon_options["simulationB1center"])"
if recon_options["simulationVariable"]=="T1"; txtt1="T1=variable, ";
elseif recon_options["simulationVariable"]=="T2"; txtt2="T2=variable, ";
elseif recon_options["simulationVariable"]=="B1"; txtB1="B1=variable";
elseif recon_options["simulationVariable"]=="measured from phantom"
    ; # no overwrite of texts
else @assert false
end
# axd.text(0.2,0.9,txtt1*txtt2*txtB1, transform=axd.transAxes) # In publication, this is for the legend

if recon_options["simulationVariable"]=="measured from phantom"
    axd.plot([0.0,0.02],[0.0,0.02],color="black")
end
@show slopes
