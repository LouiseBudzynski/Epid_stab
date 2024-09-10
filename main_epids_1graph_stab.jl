import OffsetArrays
using Plots
using Revise
using StatsBase
using ProgressMeter
using SparseArrays, LinearAlgebra, IndexedGraphs, Graphs
srcpath = "./src"
include("$srcpath/bp.jl") 
include("$srcpath/observables.jl") 
include("$srcpath/single_instance.jl") 
include("$srcpath/post_infer_1graph.jl") 

function main(args)
    arg_gam=parse(Float64, args[1])
    arg_N = parse(Int64, args[2])
    arg_maxiter=parse(Int64, args[3])
    if length(args) > 3
        nbrun=parse(Int64, args[4])
    end

    #set parameters
    T = 8 # discrete time: number of time-steps
    λp = 1.0# planted infection rate
    λi = λp #inferred infection rate
    γp = arg_gam #0.019 # planted autoinfection probability
    γi = γp # inferred autoinfection probability
    N = arg_N; #population size
    dilution = 0.0 #dilution of observations. dil=0 means everybody observed once
    fr = 0.0; #noise in the observation. 
    maxd=3
    degree_dist = Dirac(maxd) #the distribution of the degree in the graph. Dirac means random regular.
    param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist]

    if length(args) > 3
        prefix="fileres_epid_1graph_stab_T"*string(T)*"_LAM"*string(λp)*"_GAM"*string(γp)*"_N"*string(N)*"_iter"*string(arg_maxiter)*"_run"*string(nbrun)
    else
        prefix="fileres_epid_1graph_stab_T"*string(T)*"_LAM"*string(λp)*"_GAM"*string(γp)*"_N"*string(N)*"_iter"*string(arg_maxiter)
    end
    open(prefix, "w") do out
        open("err." * prefix, "w") do err
            redirect_stdout(out) do
                redirect_stderr(err) do
                    println("#param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist]=", param)
                    #init pop
                    M=ParametricModel_1graph(N=N,T=T,γp=γp,λp=λp,γi=γi,λi=λi,fr=fr,dilution=dilution,distribution=degree_dist,maxd=maxd);
                    #iterations
                    pop_dynamics_1graph_stab(M, tot_iter=arg_maxiter)
                end
            end
        end
    end
    return 0
end
main(ARGS)