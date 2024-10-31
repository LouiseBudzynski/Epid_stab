import OffsetArrays
using StatsBase
using Distributions
srcpath = "./src"
include("$srcpath/bp.jl") #functions for bp update
include("$srcpath/bp_stab.jl") #additional functions for commputing stability parameter
include("$srcpath/post_infer.jl") #main functions
include("$srcpath/observables.jl") #functions for estimating observables

function main(args)
    arg_gam=parse(Float64, args[1])
    arg_N = parse(Int64, args[2])
    arg_maxiter=parse(Int64, args[3])
    nbrun=parse(Int64, args[4])

    #set parameters
    T = 8 # discrete time: number of time-steps
    λp = 1.0# planted infection rate
    λi = λp #inferred infection rate
    γp = arg_gam #0.019 # planted autoinfection probability
    γi = γp # inferred autoinfection probability
    N = arg_N; #population size
    dilution = 0.0 #dilution of observations. dil=0 means everybody observed once
    fr = 0.0; #noise in the observation. 
    degree_dist = Dirac(3) #the distribution of the degree in the graph. Dirac means random regular.
    param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist]

    println("#param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist]=", param)
    #init pop
    M = ParametricModel(N=N,T=T,γp=γp,λp=λp,γi=γi,λi=λi,fr=fr,dilution=dilution,distribution=degree_dist);
    #iterations
    pop_dynamics(M, tot_iterations = arg_maxiter, stop_at_convergence=false);
    return 0
end
main(ARGS)