# Provide Dict containing UMCU MRI data format
function line_phantom(n)

    T₁ = 0.8*ones(n) |> collect
    T₂ = 0.06*ones(n)
    ρ = (1.0 * ones(n) ) .|> complex 
    ρˣ = real.(ρ)
    ρʸ = imag.(ρ)
    q = StructArray( map(BlochSimulators.T₁T₂ρˣρʸ , T₁, T₂, ρˣ, ρʸ) )
    return q
end

function phantom_2D(m,n)

    T₁ = 0.8*ones(m,n) |> collect
    T₂ = 0.06*ones(m,n)
    ρ = (1.0 * ones(m,n) ) .|> complex 
    ρˣ = real.(ρ)
    ρʸ = imag.(ρ)
    q = StructArray( map(BlochSimulators.T₁T₂ρˣρʸ , T₁, T₂, ρˣ, ρʸ) )
    return q
end

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



    # sl = shepp_logan(N, SheppLoganBrainWeb()) |> rotr90

    # regions = map( val -> findall(sl .== val), unique(sl))

    # T₁ = zeros(N,N)
    # T₂ = zeros(N,N)
    # ρˣ = zeros(N,N)
    # ρʸ = zeros(N,N)

    # for r in regions
    #     T₁[r] .= rand((0.3:0.01:2.5))
    #     T₂[r] .= rand((0.03:0.001:0.2))
    #     ρˣ[r] .= rand(0.5:0.02:1.5)
    #     ρʸ[r] .= rand(0.5:0.02:1.5)
    # end

    # T₂[ T₂ .> T₁ ] .= 0.5 * T₁[ T₂ .> T₁ ]

    # return map(T₁T₂ρˣρʸ, T₁, T₂, ρˣ, ρʸ) |> StructArray



function load_data_phantom(recon_options)

    vx, vy = recon_options["numphantom_size"]
    contrast = recon_options["numphantom_sequence"]
    trajectory = recon_options["numphantom_trajectory"]

    # Make phantom
    # if recon_options["numphantom_type"] == "brainweb"
    #     phantom = MRSTAT.MRITools.brainweb(vx,vy);
    # # elseif recon_options["numphantom_type"] == "shepp_logan"
    #     # phantom = MRSTAT.MRITools.shepp_logan(vx,vy);
    # elseif recon_options["numphantom_type"] == "checkerboard"
    #     @assert vx == vy
    #     blocksize = 8
    #     @assert rem(vx,blocksize) == 0
    #     phantom = MRSTAT.MRITools.checkerboard(vx,blocksize, B₁map = true);
    # elseif recon_options["numphantom_type"] == "tubes"
    #     phantom = MRSTAT.MRITools.tubes(vx,vy, B₁map=false);

    #     # jldsave("B1.jld2"; B1 = phantom.B₁)
    #     # phantom.B₁ .= 0.85
    # elseif recon_options["numphantom_type"] == "line"
    #     phantom = MRSTAT.MRITools.line(vy);
    # end

    # phantom = line_phantom(vy)
    phantom = phantom_2D(vx,vy)

    # not properly implemented at the moment
    multi_echo = false
    ne = 8 # nr echos for multiecho

    # Set coordinates (cm)
        nv = vx * vy;
        Δx = 0.1;
        Δy = 0.1;
        fovx = vx * Δx;
        fovy = vy * Δy;
        x =  -fovx/2 : Δx : fovx/2 - Δx;
        if vx == 1
            x = [0.0]
        end
        y =  -fovy/2 : Δy : fovy/2 - Δy;
        @show typeof(x)
        @show x
        #coordinates = vec(tuple.(x,y',1));
        coordinates = vec(Coordinates.(x,y',1.0));

    # Set coilmaps

        coilmaps = SVector.(ones(ComplexF64,nv));

        if recon_options["coils"] > 1
            # coil1 = complex.(repeat(LinRange(0.75,1.25,vx), 1,vy)) |> vec;
            # coil2 = im .* complex.(repeat(reverse(LinRange(0.75,1.25,vx)), 1,vy)) |> vec;

            # coilmaps = map(SVector{2}{ComplexF64}, coil1, coil2)
            nc = length(recon_options["coils"])
            coilmaps = [rand(SVector{nc}{ComplexF64}) for _ in 1:nv]
        end

    # Set trajectory

        nk = Int(ceil(recon_options["nTR"]/vy)); # nr of "kspaces"
        os = 2;  # factor two oversampling
        @info "Readout oversampling factor = $os"

        if trajectory == "Cartesian"
            @info "Linear Cartesian Trajectory"
            Δt_adc = 4e-6 / os;
            py_min = -vy÷2;
            py_max =  (vy-1)÷2;
            ns = vx * os;
            py = repeat(py_min:py_max, nk);
            nr = length(py)
            Δkˣ = 2π / (os*fovx);
            Δkʸ = 2π / fovy;
            k0 = [(-ns/2 * Δkˣ) + im * (py[r] * Δkʸ) for r in 1:nr];
            Δk = [Δkˣ + 0.0im for r in 1:nr];

            if multi_echo
                @info "Alternating readout order"
                k0 = [(ns/2 * Δkˣ)*(-1)^i + im * (py[p] * Δkʸ) for i in 1:ne, p in 1:length(py)];
                k0 = vec(k0)
                # determine step in kspace per sample for each readout
                Δk = [Δkˣ*(-1)^(i-1) + 0.0im for i in 1:ne, p in 1:length(py)];
                Δk = vec(Δk)
                nr *= ne
            end

            # @info "3x trajectory"
            # reps = 3
            # nr = reps*nr
            # k0 = repeat(k0, reps)
            # Δk = repeat(Δk, reps)
            # py = repeat(py, reps)

            # trajectory = MRSTAT.GradientTrajectories.CartesianTrajectory(nr,ns,Δt_adc,k0,Δk,py);
            trajectory = BlochSimulators.CartesianTrajectory2D(nr,ns,Δt_adc,k0,Δkˣ,py, os);

        elseif trajectory == "Radial"

            @info "Golden Angle Radial Trajectory"
            Δt_adc = 4e-6/os;

            ns = os*vx;
            nr = nk * vy
            φ = π/((√5+1)/2) # golden angle of ~111 degrees
            # @info "testing √2 thing"
            # Δkˣ = √2 * 2π / (os*fovx);
            # √2 * 2π / (os*fovx);
            Δkˣ = 2π / (os*fovx);
            # 2π / (os*fovx);

            # starting point in kspace for each readout
            radial_angles = collect(φ .* (0:nr-1))
            k0 = -(ns/2)*Δkˣ + 0.0im
            k0 = collect(@. exp(im*radial_angles) * k0)
            # k_start_readout = [(x = k.re, y = k.im) for k in tmp]

            Δk = Δkˣ + 0.0im
            Δk = collect(@. exp(im*radial_angles) * Δk)

            if multi_echo
                @warn "TODO"
            end

            trajectory = MRSTAT.GradientTrajectories.RadialTrajectory(nr,ns,Δt_adc,k0,Δk,radial_angles, os);
        elseif trajectory == "Spiral"
            @warn "TODO"
        end

    # Set sequence
        nTR = nr; # nr of TRs to be simulated
        γ = 26753.0;

        # sl = 1 # sweeplength
        # RF_train = sin.(LinRange(0,π,sl*vy)).^2
        # # RF_train = LinRange(0,90, nr)
        # RF_train = repeat(RF_train,1,nk÷sl)
        # Random.seed!(1234);
        # RF_train = RF_train .* (80 * rand(1,nk÷sl)) .+ 15

        # # make sure first sweep starts from zero otherwise aliasing artefacts?!
        # RF_train[1:end÷2,1] .= sin.(LinRange(0,π/2,sl*vy÷2)).^2 .* RF_train[end÷2,1]
        # RF_train[1,1] = RF_train[2,1] / 2 # start with alpha over 2 to prevent oscillations
        # RF_train = vec(complex.(RF_train))

        # rf_function = rfDictionary[recon_options["numphantom_rf_shape"]]
        RF_train = RF_from_file(nTR,vy)



        if contrast == "Balanced"
            @info "Balanced Sequence"
            @. RF_train[1:2:end] *= -1; # (0,180) phase cycling
            # plot(abs.(RF_train))
            TR = 0.00788; # s
            TE = TR/2; # s
            nRF = 15; # nr of samples to simulate RF excitation
            RFexdur = 0.001;
            Δt = (ex=RFexdur/nRF, inv = 0.01, pr = (TR - RFexdur)/2);
            γΔtGRz = (ex=0.00/nRF, inv = 0.00, pr = -0.0);
            γΔtRF = (pi/180) * (1/nRF) * SVector(ones(nRF)...); # RF waveform normalized to flip angle of 1 degree
            nz = 16
            z = SVector( vec(collect(LinRange(-1.0,1.0,nz)))... );
            # z = SVector(0.0)

            sequence = MRSTAT.BlochSimulators.bSSFP(RF_train, TR,TE,γΔtRF,Δt,γΔtGRz,z)

            if multi_echo
                TE = SVector( 0.00285, 0.00845, 0.01405, 0.01965, 0.02525, 0.03085, 0.03645, 0.04205); # s
                sequence = MRSTAT.BlochSimulators.bSSFP_ME(RF_train, TR,TE,γΔtRF,Δt,γΔtGRz,z)
            end

        elseif contrast == "Spoiled"
            @info "Spoiled Sequence"
            # @warn "Doesn't work well without phase cycling, but why?!"
            # @. RF_train[1:2:end] *= -1; # (0,180) phase cycling
            TR = recon_options["TR"]; # s
            TE = TR/2.0;
            max_state = recon_options["maxstate"];
            TI = (recon_options["startstate"]==1) ? 20.0 : 0.01;; # s
            nz = recon_options["slice_discretization"];
            z_locations =1:1

            if nz > 1
                # sliceprofiles
                GR_strength = 0.015914520263671874
                GR_refocus_area = -5.238291781787155e-6
                RF_duration = 0.0006304000020027161
                RF_waveform = [0.0, 0.00023196187512812663, 0.00048245823422532876, 0.0007509274262622404, 0.0010373694512388618, 0.0013417843091551927, 0.001663610348981867, 0.002003970872777617, 0.0023611809274543447, 0.0027352405130120503, 0.0031255879784213673, 0.0035316616726529305, 0.003952899944677374, 0.004388741143465331, 0.004837500315928704, 0.005299177462067493, 0.005772087628793602, 0.006255669165077663, 0.006748237117831578, 0.007248668184996616, 0.007755839064514046, 0.008268064803295769, 0.008784222099283052, 0.0093026259993878, 0.009821591550521914, 0.010339433799597297, 0.010855029444555214, 0.011366131881278205, 0.011871056156678172, 0.012369240619725747, 0.012857315364244734, 0.01333471873920577, 0.013799765791520754, 0.014250209917072225, 0.01468492781380145, 0.015102234528620332, 0.015500445108440772, 0.01587843625120404, 0.016234523003822034, 0.01656702041320666, 0.016875366828328554, 0.017157877296099616, 0.017413990165490476, 0.01764202048341304, 0.01784196824986731, 0.018012148511765184, 0.018152561269106665, 0.01826208321983302, 0.01834071436394425, 0.018387893050410987, 0.018403619279233233, 0.018387893050410987, 0.01834071436394425, 0.01826208321983302, 0.018152561269106665, 0.018012148511765184, 0.01784196824986731, 0.01764202048341304, 0.017413990165490476, 0.017157877296099616, 0.016875366828328554, 0.01656702041320666, 0.016234523003822034, 0.01587843625120404, 0.015500445108440772, 0.015102234528620332, 0.01468492781380145, 0.014250209917072225, 0.013799765791520754, 0.01333471873920577, 0.012857315364244734, 0.012369240619725747, 0.011871056156678172, 0.011366131881278205, 0.010855029444555214, 0.010339433799597297, 0.009821591550521914, 0.0093026259993878, 0.008784222099283052, 0.008268064803295769, 0.007755839064514046, 0.007248668184996616, 0.006748237117831578, 0.006255669165077663, 0.005772087628793602, 0.005299177462067493, 0.004837500315928704, 0.004388741143465331, 0.003952899944677374, 0.0035316616726529305, 0.0031255879784213673, 0.0027352405130120503, 0.0023611809274543447, 0.002003970872777617, 0.001663610348981867, 0.0013417843091551927, 0.0010373694512388618, 0.0007509274262622404, 0.00048245823422532876, 0.00023196187512812663, 0.0]

                z_locations = range(0.0, stop=0.01, length=nz)

                if recon_options["slice_profile_correction"] == "small_tip_angle"

                    @assert isodd(nz)

                    sliceprofile = MRSTAT.MRITools.small_tip_angle_approximation(GR_strength, GR_refocus_area, RF_waveform, RF_duration, z_locations)
                    sliceprofiles = repeat(sliceprofile', length(RF_train)) # same for each flip angle

                elseif recon_options["slice_profile_correction"] == "shinnar_leroux"

                    sliceprofiles = MRSTAT.MRITools.shinnarleroux_forward(RF_train, RF_waveform, RF_duration, GR_strength, GR_refocus_area, z_locations)
                end
            else
                # @error "nz = 1 gives a bug"
                sliceprofiles = ones(length(RF_train),1)
            end

            sequence = MRSTAT.BlochSimulators.FISP2D(RF_train, sliceprofiles, TR, TE, max_state, TI)

            if multi_echo
                @warn "todo"
                # TE = SVector( 0.00285, 0.00845, 0.01405, 0.01965, 0.02525, 0.03085, 0.03645, 0.04205); # s
                # sequence = BlochSimulators.SPGR_ME(RF_train, TR, TE, sliceprofile, max_state, TI)
            end
        end

    # Apply mask
        # mask = vec(abs.(phantom.ρ) .> 0.0); @warn "Proton density based mask is applied"
        mask = trues(prod(size(phantom.ρˣ))); @info "No mask is applied"
        coordinates = coordinates[mask]
        coilmaps    = coilmaps[mask]

    return sequence, coordinates, coilmaps, trajectory, mask, phantom
end