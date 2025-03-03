function calculate_ν_stab!(M,val,neighbours,xi0,oi,sympt,ci,ti_obs)
    @unpack T,γi,Λ,μ,tmpν,mom1μ = M
    tmpν .= 0
    if xi0 == 0
        for τi = 1:T+1
            for ti = 0:T+1
                ξ = obs(M,ti,τi,oi,sympt,ci,ti_obs)
                if ξ == 0.0
                    continue  
                end
                seed = ti == 0 ? γi : 1 - γi
                phi = ti == 0 || ti == T + 1 ? 0 : 1
                m1, m2, m3, m4 = one(eltype(μ)),one(eltype(μ)),one(eltype(μ)),one(eltype(μ))
                for k in neighbours 
                    if k==val
                        m1 *= mom1μ[ti,1,τi,1,k] + mom1μ[ti,1,τi,2,k]
                        m2 *= mom1μ[ti,0,τi,1,k] + mom1μ[ti,0,τi,2,k]
                        m3 *= mom1μ[ti,1,τi,2,k]
                        m4 *= mom1μ[ti,0,τi,2,k]
                    else
                        m1 *= μ[ti,1,τi,1,k] + μ[ti,1,τi,2,k]
                        m2 *= μ[ti,0,τi,1,k] + μ[ti,0,τi,2,k]
                        m3 *= μ[ti,1,τi,2,k]
                        m4 *= μ[ti,0,τi,2,k]
                    end
                end
                for tj=0:T+1                
                    tmpν[ti,tj,τi,1] = ξ  * seed * (Λ[ti-tj-1] * m1 - phi * Λ[ti-tj] * m2)
                    tmpν[ti,tj,τi,2] = tmpν[ti,tj,τi,1] + ξ * (τi<T+1) * seed * (phi * Λ[ti-tj] * m4 - Λ[ti-tj-1] * m3)
                end
            end
        end
    else
        for tj = 0:T+1
            for ti = 0:T+1
                ξ = obs(M,ti,0,oi,sympt,ci,ti_obs)
                if ξ == 0.0  #if the observation is NOT satisfied
                    continue
                end
                seed = (ti==0 ? γi : (1-γi) )
                phi = (ti==0 || ti==T+1) ? 0 : 1
                m1, m2 = one(eltype(μ)), one(eltype(μ))
                for k in neighbours    
                    if k==val 
                        m1 *= mom1μ[ti,1,0,0,k] + mom1μ[ti,1,0,1,k] + mom1μ[ti,1,0,2,k]
                        m2 *= mom1μ[ti,0,0,0,k] + mom1μ[ti,0,0,1,k] + mom1μ[ti,0,0,2,k]
                    else
                        m1 *= μ[ti,1,0,0,k] + μ[ti,1,0,1,k] + μ[ti,1,0,2,k]
                        m2 *= μ[ti,0,0,0,k] + μ[ti,0,0,1,k] + μ[ti,0,0,2,k]
                    end           
                end
                tmpν[ti,tj,0,:] .= ξ * seed * (Λ[ti-1-tj] * m1 - phi * Λ[ti-tj] * m2)
            end
        end
    end
    if any(isnan.(tmpν))
        println("NaN tmpν at $(M.λi), $(M.γi), $(popsize(M)), $(M.fr)")
        exit(0)
    end 
end

function original_L2norm(M,ν,sji)
    # first we sum on ti and tj
    tmp = sum(sum(ν.*ν,dims=1),dims=2)
    norm = zero(eltype(ν))
    T = M.T
    # for fixed taui and disorder sji it is possible to compute tauj for every sigma
    # so we sum over tauj at each step of the cycle
    for taui = 0:T+1
        norm += max(0,taui-sji) * tmp[0,0,taui,0] + (taui-sji >= 0) * tmp[0,0,taui,1] + (T+2 - max(taui-sji+1,0)) * tmp[0,0,taui,2]
    end
    norm/=(T+2)^4
    return sqrt(norm)
end

function update_mom1μ!(M,l,sij,sji)
    @unpack T,Λ,mom1μ,Paux,mom1ν = M
    mom1μ[:,:,:,:,l] .= 0
    Σ = cumsum(mom1ν,dims=3)
    @inbounds for tj = 0:T+1
        for τj = 0:T+1
            Paux .= 0
            for ti = 0:T+1
                Γ = Σ[ti,tj,min(τj+sji-1,T+1),2] - (τj-sij>=0)*Σ[ti,tj,max(τj-sij,0),2]+(τj+sji<=T+1)*mom1ν[ti,tj,min(τj+sji,T+1),1]+
                    Σ[ti,tj,T+1,0] - Σ[ti,tj,min(τj+sji,T+1),0]
                
                for c = 0:1
                    Paux[c,0] += Λ[tj-ti-c] * (τj-sij-1>=0) * Σ[ti,tj,max(τj-sij-1,0),2]
                    Paux[c,1] += Λ[tj-ti-c] * (τj-sij>=0) * mom1ν[ti,tj,max(τj-sij,0),2]
                    Paux[c,2] += Λ[tj-ti-c] * Γ
                end
            end
            mom1μ[tj,:,τj,:,l] = Paux
        end
    end
    if any(isnan,mom1μ[:,:,:,:,l])
        println("NaN in mom1μ")
        exit(0)
    end
end
