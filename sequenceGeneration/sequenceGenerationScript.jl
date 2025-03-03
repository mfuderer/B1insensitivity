# sequenceGenerationScript.
# Note that, within the script, one has to manually cycle the variable "case" from 1 to 4. 
# The processing time is vastly different between case=1 (Amplitude-only noise-optimized), which may take a minute, 
#     and case=4 (Amplitude+Phase, B1-opt), which may take up to a day.
# The results should appear in the folder RFsequences (well, they are actually there, but if you want to reproduce these, you can move them to another location and generate them again using the script).

case = 1 # 1 or 2, 3, 4;  to be cycled manually
caseName = ["NoPhase_noise", "NoPhase_B1opt", "Phase_noise", "Phase_B1opt"]
b1handling = ["no",          "sensitivity",   "no",          "sensitivity"]
qPhase    = [false,          false,           true,          true]

include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 224                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
recon_options["optpars"]   = Optim.Options(time_limit = 20000.0, iterations = 100000, f_tol=1.0e-5, g_tol = 1.0e-5)  
# @warn "debugging" 
# recon_options["optpars"]   = Optim.Options(time_limit = 20.0, iterations = 100000, f_tol=1.0e-5, g_tol = 1.0e-5)  

recon_options["opt_criterion"] = "noise_level" 
recon_options["account_SAR"]   = true     

recon_options["sar_limit"] = 40^2/0.01 
recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = b1handling[case] # "no" # "sensitivity" # "no" # "sensitivity" #
recon_options["lambda_B1"] = 10.0     
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_focus"] = "max"      
recon_options["opt_complex"] = false      
recon_options["opt_account_maxFlip"] = false
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] = qPhase[case] # false # true                         
recon_options["considerCyclic"] = false  
recon_options["opt_emergeCriterion"] = 500 # 2000 # 500 # 2000
ph = [] 
# ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1
recon_options["sizeSteps"] = [6]  
recon_options["B1metric"] = "multi_point_values" # "multi_point"         
nRealizations = 10  
# nRealizations = 1 ; @warn "temporary debugging line, has to be removed later" 
recon_options["rfFolder"]= "./RFsequences/"    

fn_base = caseName[case]
# fn_base = "temp"; @warn "temporary debugging line, has to be removed later"
for i in 1:nRealizations; goodseed = i
    stageText = ""
    portionRange = 0:0
    fn = recon_options["rfFolder"]*"$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, goodseed);
    FileIO.save(fn,"RFdeg",RFdeg)
    @show fn
end 

#
@show recon_options["T1T2set"]
include("BLAKJac_B1_figures.jl")
testCase = caseName[case];   Nr=10
map_B1_sensitivities(testCase,5)
@show recon_options["T1T2set"]

recon_options["handleB1"] = "sensitivity" 
recon_options["B1metric"] = "multi_point_values"   
winner=0; optscore = Inf       
figure()
for rrr in 1:Nr
    fn = "$testCase($rrr)"
    recon_options["rfFile"]  = fn
    RFdeg = RF_from_file(recon_options["nTR"], recon_options["nky"])
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options)
    mn = maximum(noises[2:3]); b1m=mean(b1f); b1max = maximum(b1f)
    mixCrit = sqrt(mn^2+(10*b1max)^2)
    @printf("Try %d: n=%.2f, T1B1=%.3f, T2B1=%.3f, criterium=%.2f \n", rrr, mn, b1f[2], b1f[3], mixCrit)
    if mixCrit < optscore
        optscore = mixCrit
        winner = rrr
    end
    scatter(mn,rrr, color="blue")   
    scatter(mixCrit,rrr, color="red")
    scatter(10*b1f[3],rrr, color="green")
end
scatter(optscore, winner, color="red", s=100)
@printf("Ant the winning sequence is ... %d \n", winner)
