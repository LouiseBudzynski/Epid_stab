import OffsetArrays
using Plots
using Revise
using StatsBase
using ProgressMeter
using SparseArrays, LinearAlgebra, IndexedGraphs, Graphs
using Distributions
using DelimitedFiles
srcpath = "./src"
include("$srcpath/bp.jl") #julia code which contains the bp update
include("$srcpath/post_infer.jl") #definition of the main functions and the types
include("$srcpath/observables.jl") # julia code for estimating observables

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
    degree_dist = Dirac(3) #the distribution of the degree in the graph. Dirac means random regular.
    param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist]

    if length(args) > 3
        prefix="fileres_epid_stab_T"*string(T)*"_LAM"*string(λp)*"_GAM"*string(γp)*"_N"*string(N)*"_iter"*string(arg_maxiter)*"_run"*string(nbrun)
        fileobs="obs_epid_stab_T"*string(T)*"_LAM"*string(λp)*"_GAM"*string(γp)*"_N"*string(N)*"_iter"*string(arg_maxiter)*"_run"*string(nbrun)
    else
        prefix="fileres_epid_stab_T"*string(T)*"_LAM"*string(λp)*"_GAM"*string(γp)*"_N"*string(N)*"_iter"*string(arg_maxiter)
        fileobs="obs_epid_stab_T"*string(T)*"_LAM"*string(λp)*"_GAM"*string(γp)*"_N"*string(N)*"_iter"*string(arg_maxiter)
    end
    open(prefix, "w") do out
        open("err." * prefix, "w") do err
            redirect_stdout(out) do
                redirect_stderr(err) do
                    println("#param=[T, λp, λi, γp, γi, N, dilution, fr, degree_dist]=", param)
                    #init pop
                    M = ParametricModel(N=N,T=T,γp=γp,λp=λp,γi=γi,λi=λi,fr=fr,dilution=dilution,distribution=degree_dist);
                    #print observables               
                    open(fileobs, "a") do io
                        writedlm(io, [vcat(["#param"],param)])
                         writedlm(io, ["#column 1=iter; line 1=planted; line 2=inferred; line 3=AUC"])
                    end                
                    #iterations
                    F,it = pop_dynamics_stab(M, filepr = fileobs, tot_iterations = arg_maxiter);
                end
            end
        end
    end
    return 0
end
main(ARGS)