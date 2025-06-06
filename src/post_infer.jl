"""
    ParametricModel{D,D2,Taux,M,M1,M2,O,Tλ,MO,Tsympt}

Mutable struct which represents a model of inference. It contains the planted and inference parameters, the graph distribution, the BP messages. Is the fundamental object of the code.
"""
mutable struct ParametricModel{D,D2,Taux,M,M1,M2,O,Tλ,MO,Tsympt}
    T::Int
    γp::Float64
    λp::Float64
    γi::Float64
    λi::Tλ
    σ0::Float64
    Paux::Taux
    μ::M
    mom1μ::M
    belief::M2
    ν::M1
    tmpν::M1
    mom1ν::M1
    fr::Float64
    distribution::D
    residual::D2
    Λ::O
    obs_list::MO
    obs_range::UnitRange{Int64}
    field::Float64
    p_sympt_pla::Float64
    p_test_pla::Float64
    p_sympt_inf::Tsympt
    p_test_inf::Float64
end

function ParametricModel(; N, T, γp, λp, γi=γp, λi=λp, σ0=1e-4, fr=0.0, dilution=0.0, distribution, obs_range=T:T,field=0.0,p_sympt_pla=0.0,p_test_pla=1-dilution,p_sympt_inf=p_sympt_pla,p_test_inf=p_test_pla)
    μ = fill(one(λi * p_sympt_inf) / (6*(T+2)^2), 0:T+1, 0:1, 0:T+1, 0:2, 1:N)
    mom1μ = OffsetArray(rand(-1:2:1, T+2, 2, T+2, 3, N).*σ0, 0:T+1, 0:1, 0:T+1, 0:2, 1:N)
    belief = fill(zero(λi * p_sympt_inf), 0:T+1, 0:T+1, N)
    ν = fill(zero(λi * p_sympt_inf), 0:T+1, 0:T+1, 0:T+1, 0:2)
    tmpν = fill(zero(λi * p_sympt_inf), 0:T+1, 0:T+1, 0:T+1, 0:2)
    mom1ν = fill(zero(λi * p_sympt_inf), 0:T+1, 0:T+1, 0:T+1, 0:2)
    Paux = fill(zero(λi * p_sympt_inf), 0:1, 0:2)
    Λ = OffsetArray([t <= 0 ? one(λi) : (1-λi)^t for t = -T-2:T+1], -T-2:T+1)
    ParametricModel(T, γp, λp,γi, λi, σ0,Paux, μ, mom1μ, belief, ν, tmpν, mom1ν,fr, distribution, residual(distribution), Λ,fill(false,N),obs_range,field,p_sympt_pla,p_test_pla,p_sympt_inf,p_test_inf)
end


#=function avg_err(b)
    N = size(b,3)
    T = size(b,1) - 2
    avg_bel = reshape(sum(sum(b,dims=2),dims=3) ./ (N*(T+2)),T+2) 
    err_bel = sqrt.(reshape(sum(sum(b .^ 2,dims=2),dims=3) ./ (N * (T+2)),T+2) .- avg_bel .^ 2) ./ sqrt(N)
    return avg_bel, err_bel
end=#

function avg_err(b)
    N = size(b,3)
    T = size(b,1) - 2
    avg_bel = reshape(sum(b, dims=3) ./ N, T+2, T+2) 
    err_bel = sqrt.(reshape(sum(b .^ 2, dims=3) ./ N, T+2, T+2) .- avg_bel .^ 2) ./ sqrt(N)
    return avg_bel, err_bel
end

function FatTail(support::UnitRange{Int64}, exponent,a)
    p = 1 ./ (collect(support) .^ exponent .+ a)
    return DiscreteNonParametric(collect(support), p / sum(p))
end

function sweep_stab!(M, iter, maxiter; nbstab = round(maxiter/2))
    N = popsize(M)
    F = 0.0
    Δt = 0.0
    for l = 1:N
        # update μ and mom1μ
        xi0,sij,sji,d,oi,sympt,ci,ti_obs = rand_disorder(M,M.distribution)
        neighbours = rand(1:N,d-1)
        calculate_ν!(M,neighbours,xi0,oi,sympt,ci,ti_obs)
        zψij = original_normalization(M,M.ν,sji)
        F -= 0.5*d*log(zψij)
        update_μ!(M,l,sij,sji) 
        if (iter > maxiter - nbstab)
            M.mom1ν.=0.0
            for val in neighbours
                calculate_ν_stab!(M,val,neighbours,xi0,oi,sympt,ci,ti_obs)
                zψijm = original_normalization(M,M.tmpν,sji)
                M.mom1ν .+= M.tmpν.-M.ν.*(zψijm/zψij)
            end
            M.mom1ν ./= zψij
            delt=original_L2norm(M,M.mom1ν,sji)
            Δt+=delt
            update_mom1μ!(M,l,sij,sji)
        end

        # update belief
        xi0,sij,sji,d,oi,sympt,ci,ti_obs = rand_disorder(M,M.distribution)
        M.obs_list[l] = oi 
        neighbours = rand(1:N,d)
        zψi = calculate_belief!(M,l,neighbours,xi0,oi,sympt,ci,ti_obs)
        F+=(0.5*d-1)*log(zψi)
    end
    if (iter > maxiter - nbstab)
        Δt = Δt/N
    else
        Δt=-1
    end
    return F / N, Δt
end

function sweep!(M)
    e = 1 #edge counter
    N = popsize(M)
    F_itoj = zero(M.λi * M.p_sympt_inf)
    Fψi = zero(M.λi * M.p_sympt_inf)
    for l = 1:N
        # Extraction of disorder: state of individual i: xi0, delays: sij and sji
        xi0,sij,sji,d,oi,sympt,ci,ti_obs = rand_disorder(M,M.distribution)
        M.obs_list[l] = oi #this is stored for later estimation of AUC
        neighbours = rand(1:N,d)
        for m = 1:d
            res_neigh = [neighbours[1:m-1];neighbours[m+1:end]]
            calculate_ν!(M,res_neigh,xi0,oi,sympt,ci,ti_obs)
            #from the un-normalized ν message it is possible to extract the orginal-message 
            #normalization z_i→j 
            # needed for the computation of the Bethe Free energy
            r = 1.0 / log(1-M.λp)
            sij = floor(Int,log(rand())*r) + 1
            sji = floor(Int,log(rand())*r) + 1
            zψij = original_normalization(M,M.ν,sji)
            F_itoj += log(zψij) 
            # Now we use the ν vector just calculated to extract the new μ.
            # We overwrite the μ in postition μ[:,:,:,:,l]
            update_μ!(M,e,sij,sji)  
            e = mod(e,N) + 1
        end
        zψi = calculate_belief!(M,l,neighbours,xi0,oi,sympt,ci,ti_obs)
        Fψi += (0.5 * d - 1) * log(zψi) 
    end
    return (Fψi - 0.5 * F_itoj) / N 
end


function pop_dynamics_stab(M; tot_iterations = 10, nbstab = round(tot_iterations/2), observ = false)
    N, T = popsize(M), M.T

    #uniform initial condition
    M.ν.=1/((T+1)^4)
    for l in 1:N
        (xi0, sij, sji, d, oi, sympt, ci, ti_obs) = rand_disorder(M, M.distribution)
        update_μ!(M,l,sij, sji)
    end
    if observ
        println("#1.iter 2.F 3.Δ 4.MMO[t=0] 5.AUC[t=0]")
    else
        println("#1.iter 2.F 3.Δ")
    end
    for iterations = 1:tot_iterations
        F, Δ = sweep_stab!(M, iterations, tot_iterations, nbstab=nbstab) 
        if observ
            MMO = avgOverlap(M.belief)
            if iterations >= tot_iterations-5
                ensAUC = avgAUC(M.belief,M.obs_list,count_obs=true)
                println(iterations, "\t", F, "\t", Δ, "\t", MMO[1], "\t", ensAUC[1])
            else
                println(iterations, "\t", F, "\t", Δ, "\t", MMO[1])
            end
        else
            println(iterations, "\t", F, "\t", Δ)
        end
        flush(stdout)
        if (iterations == tot_iterations - nbstab) 
            M.mom1μ .= (M.mom1μ).*(M.μ)
        end
    end
end

function pop_dynamics(M; tot_iterations = 5, tol = 1e-10, eta = 0.0, infer_lam=false, infer_gam = false,infer_sympt=false,nonlearn_iters=0,stop_at_convergence=true)
    N, T, F = popsize(M), M.T, zero(M.λi * M.p_sympt_inf)
    F_window = zeros(10)
    converged = false
    lam_window = zeros(10)
    gam_window = zeros(10)
    sympt_window = zeros(10)
    println("#1.iter 2.F")
    for iterations = 0:tot_iterations-1
        wflag = mod(iterations,10)+1
        lam_window[wflag] = M.λi |> real
        gam_window[wflag] = M.γi |> real
        sympt_window[wflag] = M.p_sympt_inf |> real
        F = sweep!(M) 
        println(iterations, "\t", F)
        flush(stdout)
        F_window[mod(iterations,length(F_window))+1] = (F |> real)
        infer_lam = check_prior(iterations, infer_lam, lam_window, eta, nonlearn_iters)
        infer_gam = check_prior(iterations, infer_gam, gam_window, eta, nonlearn_iters)
        infer_sympt = check_prior(iterations, infer_sympt, sympt_window, eta, nonlearn_iters)
        (iterations > length(F_window)) && (converged = check_convergence(F_window,2/sqrt(N),stop_at_convergence))
        if converged  & !infer_lam & !infer_gam & !infer_sympt #if we don't have to infer 
            return F_window |> real, iterations+1
        end
        (infer_lam) && (update_lam!(M,F,eta))
        (infer_gam) && (update_gam!(M))
        (infer_sympt) && (update_sympt!(M,F,eta))
    end
    return F_window |> real , tot_iterations
end



function check_prior(iterations, infer, window, eta,nonlearn)
    if (iterations >= nonlearn && infer) # if we need to infer + we we have iterated more than the nonlearn treshold
        avg = sum(window) / 10
        err = sqrt(sum(window .^ 2)/10 - avg ^ 2)
        if (err/avg <= eta) 
           return false
        end
    end
    return infer
end

function check_convergence(window,tol,stop_at_convergence)
    l = length(window)
    avg = sum(window) / l
    variance = sum(window .^ 2)/l - avg ^ 2
    err = abs(variance) > 1e-15 ? sqrt(variance) : 0.0
    if (err <= tol) 
        return stop_at_convergence
    else 
        return false
    end
end

