using Distributions,UnPack,OffsetArrays


popsize(M) = size(M.belief,3)


function obs(M, ti, τi, oi, sympt, ci, ti_obs) 
    xt_inf = (ti <= ti_obs) 
    xt_pla = (τi <= ti_obs)
    o = oi | (xt_pla & sympt) #the indiv is either observed random or because it is I and symptomatic
    p_test = xt_inf ? (M.p_sympt_inf + M.p_test_inf * (1 - M.p_sympt_inf)) : M.p_test_inf * one(M.p_sympt_inf) 
    #the one(..) above is to stabilize types
    if ci #if ci=1 we see the opposite of the value xt_pla
        ob_made = o ? (xt_inf == !xt_pla ? ((1.0 - M.fr) * p_test) : (M.fr * p_test)) : (1.0 - p_test)
    else 
        ob_made = o ? (xt_inf == xt_pla ? ((1.0 - M.fr) * p_test) : (M.fr * p_test)) : (1.0 - p_test)
    end
    return (1 - M.field) * ob_made + M.field * Int(ti == τi)
end


function calculate_ν!(M,neighbours,xi0,oi,sympt,ci,ti_obs)
    @unpack T,γi,Λ,μ,ν = M
    ν .= 0
    if xi0 == 0
        for τi = 1:T+1
            for ti = 0:T+1
                #first we check consistency between
                # the planted time τi and the inferred 
                #time ti by checking the observation constraint
                ξ = obs(M,ti,τi,oi,sympt,ci,ti_obs)
                if ξ == 0.0 #if the observation is NOT satisfied
                    continue  # ν = 0
                end
                #Since they both depend on ti only,
                # we precaclulate the prior seed probability
                # of the individual and the value of phi function 
                # which is 1 if 0<ti<T+1 and 0 if ti=0,T+1
                seed = ti == 0 ? γi : 1 - γi
                phi = ti == 0 || ti == T + 1 ? 0 : 1
                #now we calculate the four products over
                # μ functions that we need to put in the
                # expression of ν. We call them m1,..,m4
                m1, m2, m3, m4 = one(eltype(μ)),one(eltype(μ)),one(eltype(μ)),one(eltype(μ))
                # we initialize the m's to one and then we 
                # loop a product over neighbours
                for k in neighbours 
                    m1 *= μ[ti,1,τi,1,k] + μ[ti,1,τi,2,k]
                    m2 *= μ[ti,0,τi,1,k] + μ[ti,0,τi,2,k]
                    m3 *= μ[ti,1,τi,2,k]
                    m4 *= μ[ti,0,τi,2,k]
                end
                #Now we have everything to calculate ν
                for tj=0:T+1                
                    ν[ti,tj,τi,1] = ξ  * seed * (Λ[ti-tj-1] * m1 - phi * Λ[ti-tj] * m2)
                    # We use the fact that ν for σ=2 is just ν at σ=1 plus a term
                    ν[ti,tj,τi,2] = ν[ti,tj,τi,1] + ξ * (τi<T+1) * seed * (phi * Λ[ti-tj] * m4 - Λ[ti-tj-1] * m3)
                end
            end
        end
    else
        # We are now in the case in which the individual is 
        # the zero patient. In this case the computation of 
        # the ν function is a little bit different than before
        # so we separated the cases

        for tj = 0:T+1
            for ti = 0:T+1
                ξ = obs(M,ti,0,oi,sympt,ci,ti_obs)
                if ξ == 0.0  #if the observation is NOT satisfied
                    continue
                end
                #we can calculate ν now because it is constant
                # in σ and is nonzero only if τi=0

                #As before we pre-calculate ti-dependent quantities 
                seed = (ti==0 ? γi : (1-γi) )
                phi = (ti==0 || ti==T+1) ? 0 : 1
                # We perform the product over neighbours
                m1, m2 = one(eltype(μ)), one(eltype(μ))
                for k in neighbours                
                    m1 *= μ[ti,1,0,0,k] + μ[ti,1,0,1,k] + μ[ti,1,0,2,k]
                    m2 *= μ[ti,0,0,0,k] + μ[ti,0,0,1,k] + μ[ti,0,0,2,k]
                end
                #We calculate ν in the zero patient case
                ν[ti,tj,0,:] .= ξ * seed * (Λ[ti-1-tj] * m1 - phi * Λ[ti-tj] * m2)
            end
        end
    end
    if any(isnan.(ν))
        #println("NaN ν at $(M.λi), $(M.γi), $(popsize(M)), $(M.fr)")
        return
    end 
    if sum(ν) == zero(eltype(ν))
        println("sum-zero ν at $(M.λi), $(M.γi), $(popsize(M)), $(M.fr)")
        return
    end        
end


function calculate_belief!(M,l,neighbours,xi0,oi,sympt,ci,ti_obs) 
    @unpack T, belief, γi, μ = M
    belief[:,:,l] .= zero(eltype(belief))
    if xi0 == 0
        for τi = 1:T+1
            for ti = 0:T+1
                #first we check consistency between
                # the planted time τi and the inferred 
                #time ti by checking the observation constraint
                ξ = obs(M,ti,τi,oi,sympt,ci,ti_obs)
                if ξ == 0.0 #if the observation is NOT satisfied
                    continue  # ν = 0
                end
                #Since they both depend on ti only,
                # we precaclulate the prior seed probability
                # of the individual and the value of phi function 
                # which is 1 if 0<ti<T+1 and 0 if ti=0,T+1
                seed = (ti==0 ? γi : (1-γi) )
                phi = (ti==0 || ti==T+1) ? 0 : 1
                #now we calculate the four products over
                # μ functions that we need to put in the
                # expression of ν. We call them m1,..,m4
                m1, m2, m3, m4 = one(eltype(μ)),one(eltype(μ)),one(eltype(μ)),one(eltype(μ))
                # we initialize the m's to one and then we 
                # loop a product over neighbours
                for k in neighbours 
                    m1 *= μ[ti,1,τi,1,k] + μ[ti,1,τi,2,k]
                    m2 *= μ[ti,0,τi,1,k] + μ[ti,0,τi,2,k]
                    m3 *= μ[ti,1,τi,2,k]
                    m4 *= μ[ti,0,τi,2,k]
                end
                #Now we have everything to calculate ν
                    # We use the fact that ν for σ=2 is just ν at σ=1 plus a term
                belief[ti,τi,l] = ξ  * seed * ( m1 - phi * m2) + ξ * (τi<T+1) * seed * (phi *  m4 -  m3)
            end
        end
    else
        # We are now in the case in which the individual is 
        # the zero patient. In this case the computation of 
        # the ν function is a little bit different than before
        # so we separated the cases

        for ti = 0:T+1
            ξ = obs(M,ti,0,oi,sympt,ci,ti_obs)
            if ξ == 0.0  #if the observation is NOT satisfied
                continue
            end
            #we can calculate ν now because it is constant
            # in σ and is nonzero only if τi=0

            #As before we pre-calculate ti-dependent quantities 
            seed = (ti==0 ? γi : (1-γi) )
            phi = (ti==0 || ti==T+1) ? 0 : 1
            # We perform the product over neighbours
            m1, m2 = one(eltype(μ)),one(eltype(μ))
            for k in neighbours                
                m1 *= μ[ti,1,0,0,k] + μ[ti,1,0,1,k] + μ[ti,1,0,2,k]
                m2 *= μ[ti,0,0,0,k] + μ[ti,0,0,1,k] + μ[ti,0,0,2,k]
            end
            #We calculate ν in the zero patient case
            belief[ti,0,l] = ξ * seed * ( m1 - phi *  m2)
        end
    end
    S = sum(@view belief[:,:,l])
    if S == zero(eltype(belief))
        println("sum-zero belief  at $(M.λi), $(M.γi), $(M.γp)")
        return
    end    
    belief[:,:,l] ./= S
    return S
end




residual(d::Poisson) = d #residual degree of poiss distribution is poisson with same param
residual(d::Dirac) = Dirac(d.value - 1) #residual degree of rr distribution (delta) is a delta at previous vale
residual(d::DiscreteNonParametric) = DiscreteNonParametric(support(d) .- 1, (probs(d) .* support(d)) / sum(probs(d) .* support(d)))
residual(d::DiscreteUniform) = Dirac(0)

function rand_disorder(M, dist)
    @unpack γp, λp, fr, obs_range, p_sympt_pla,p_test_pla = M
    r = 1.0 / log(1-λp) #random number to generate the delays
    sij = floor(Int,log(rand())*r) + 1
    sji = floor(Int,log(rand())*r) + 1 # infection delay
    xi0 = (rand() < γp); #zero patient
    d = rand(dist) #parameter of the degree distribution
    oi = rand() < p_test_pla # oi = 1 if the particle is observed, oi = 0 if the particle is not observed 
    sympt = rand() < p_sympt_pla #the symptom appears 
    ci = rand() < fr #ci is the corruption bit. If ci==1 then the observation is corrupted (i.e. is flipped)
    ti_obs = rand(obs_range) 
    return xi0, sij, sji, d, oi, sympt, ci, ti_obs
end


function original_normalization(M,ν,sji)
    # first we sum on ti and tj
    tmp = sum(sum(ν,dims=1),dims=2)
    norm = zero(eltype(ν))
    T = M.T
    # for fixed taui and disorder sji it is possible to compute tauj for every sigma
    # so we sum over tauj at each step of the cycle
    for taui = 0:T+1
        norm += max(0,taui-sji) * tmp[0,0,taui,0] + (taui-sji >= 0) * tmp[0,0,taui,1] + (T+2 - max(taui-sji+1,0)) * tmp[0,0,taui,2]
    end
    return norm
end

function Γ(Σ,ν,ti,tj,τj,sji,sij)
    return Σ[ti,tj,min(τj+sji-1,T+1),2] - (τj-sij>=0)*Σ[ti,tj,max(τj-sij,0),2]+(τj+sji<=T+1)*ν[ti,tj,min(τj+sji,T+1),1] + Σ[ti,tj,T+1,0] - Σ[ti,tj,min(τj+sji,T+1),0]
end

function inferred_times_msg!(ms,pos,M,sij,sji,xi0)
    T = M.T
    ν = M.ν
    ms[:,:,:,:,pos] .= 0
    Σ = cumsum(ν,dims=3)
    #ti,tj,taui,v
    for ti = 0:T+1
        for tj = 0:T+1
            if xi0 == 1
                taui = 0
                f = Γ(Σ,ν,tj,ti,taui,sij,sji) # first term of the sum
                s = (taui-sji>=0) * ν[tj,ti,max(taui-sji,0),2] # second term of the sum
                t = (taui-sji-1>=0) * Σ[tj,ti,max(taui-sji-1,0),2] # third term of the sum
                ms[ti,tj,taui,0,pos] = f + s + t
            elseif xi0 == 0
                for taui = 1:T+1
                    f = Γ(Σ,ν,tj,ti,taui,sij,sji) # first term of the sum
                    s = (taui-sji>=0) * ν[tj,ti,max(taui-sji,0),2] # second term of the sum
                    ms[ti,tj,taui,1,pos] = f + s
                    ms[ti,tj,taui,2,pos] =  f
                end
            end
        end
    end
end

function update_μ!(M,l,sij,sji)
    @unpack T,Λ,μ,Paux,ν = M
    μ[:,:,:,:,l] .= 0
    # First we calculate and store the cumulated of ν with respect to 
    # planted time, i.e. the third argument. We call Σ this cumulated 
    Σ = cumsum(ν,dims=3)
    @inbounds for tj = 0:T+1
        for τj = 0:T+1
            #First of all we set to 0 the function we want to update
            #because later we want to sum over it
            Paux .= 0
            for ti = 0:T+1
                #we pre calculate the value of the summed part
                # so not to calculate it twice
                Γ = Σ[ti,tj,min(τj+sji-1,T+1),2] - (τj-sij>=0)*Σ[ti,tj,max(τj-sij,0),2]+(τj+sji<=T+1)*ν[ti,tj,min(τj+sji,T+1),1]+
                    Σ[ti,tj,T+1,0] - Σ[ti,tj,min(τj+sji,T+1),0]
                
                for c = 0:1
                    Paux[c,0] += Λ[tj-ti-c] * (τj-sij-1>=0) * Σ[ti,tj,max(τj-sij-1,0),2]
                    Paux[c,1] += Λ[tj-ti-c] * (τj-sij>=0) * ν[ti,tj,max(τj-sij,0),2]
                    Paux[c,2] += Λ[tj-ti-c] * Γ
                end
            end
            μ[tj,:,τj,:,l] = Paux
        end
    end
    S = sum(@view μ[:,:,:,:,l])
    if iszero(S)
        println("sum-zero μ  at $(M.λi), $(M.γi)")
        @show(ν)
        return
    end 
    if any(isnan,μ[:,:,:,:,l])
        println("NaN in μ")
        @show sum(ν)
    end
    μ[:,:,:,:,l] ./= S
end


function update_measure!(F∂i,F∂iold,pos,d,msg,M)
    T = M.T
    F∂iold .= F∂i
    F∂i .= 0
    for ti = 0:T+1
        for tk = 0:T+1
            θik = (ti - tk - 1 >= 0)
            tik1 = θik ? (ti - tk - 1) : 0
            for S1 = tik1 : d * T
                for S2 = θik : d
                    for v = 0:2
                        for taui = 0:T+1
                            F∂i[ti,S1,S2,v,taui] += msg[ti,tk,taui,v,pos] * F∂iold[ti,S1 - tik1, S2 - θik,v,taui]
                        end
                    end
                end
            end
        end
    end
end

function psi(M,ti,S1,S2)
    seed = (ti == 0 ? M.γi : 1 - M.γi)
    return seed * ((1 - M.λi) ^ S1) * (1 - (1 <= ti <= M.T) * (1 - M.λi) ^ S2)
end


function unif_initializ!(M)
    N = popsize(M)
    M.ν .= 1/(T+2)^2
    for l = 1:N
        xi0,sij,sji,d,oi,sympt,ci,ti_obs = rand_disorder(M,M.distribution)
        update_μ!(M,l,sij,sji)
    end
end
