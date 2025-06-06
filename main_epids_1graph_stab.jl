import OffsetArrays
using StatsBase
using Graphs, IndexedGraphs
using Distributions,UnPack,OffsetArrays
srcpath = "./src"
include("$srcpath/post_infer_1graph.jl") 

function main(args)
    arg_gam=parse(Float64, args[1])
    arg_sig0=parse(Float64, args[2])
    arg_N = parse(Int64, args[3])
    arg_maxiter=parse(Int64, args[4])
    fixseeds=parse(Bool, args[5])
    nbrun=parse(Int64, args[6])

    #set parameters
    T = 8 # discrete time: number of time-steps
    λp = 1.0 # planted infection rate
    λi = λp # inferred infection rate
    γp = arg_gam # planted autoinfection probability
    γi = γp # inferred autoinfection probability
    σ0 = arg_sig0 # magnitude IC for perturbation
    N = arg_N; # population size
    dilution = 0.0 # dilution of observations. dil=0 means everybody observed once
    fr = 0.0; # noise in the observation. 
    maxd=3
    degree_dist = Dirac(maxd) # distribution of the degree in the graph (Dirac = random regular).
    param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist, fixseeds]

    println("#param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist, fixseeds]=", param)
    #init pop
    M=ParametricModel_1graph(N=N,T=T,γp=γp,λp=λp,fixseednumber=fixseeds,γi=γi,λi=λi,σ0=σ0,fr=fr,dilution=dilution,distribution=degree_dist,maxd=maxd);
    #iterations
    pop_dynamics_1graph_stab(M, tot_iter=arg_maxiter)
    return 0
end
main(ARGS)
