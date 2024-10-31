using Random

#mutable struct ParametricModel_1graph{D,D2,M,M1,M2,O,Tλ}
mutable struct ParametricModel_1graph{M,M1,M2,O,Tλ}
    N::Int
    Nedges::Int
    T::Int
    γp::Float64
    λp::Float64
    γi::Float64
    λi::Tλ
    μ::M
    mom1μ::M
    belief::M2
    ν::M1
    tmpν::M1
    mom1ν::M1
    fr::Float64
    Λ::O
    Neigh::Vector{Vector{Int64}}
    Observations::BitVector
end

function ParametricModel_1graph(; N, T, γp, λp, γi=γp, λi=λp, fr=0.0, dilution=0.0, distribution, maxd)
    Λ = OffsetArray([t <= 0 ? 1.0 : (1-λi)^t for t = -T-2:T+1], -T-2:T+1)
    μ = fill(0.0, 0:T+1, 0:1, maxd, N)
    for i in 1:N
        for deg in 1:maxd
            for ti in 0:T+1
                for tj in 0:T+1
                    μ[tj,0,deg,i]+=Λ[tj-ti]/((T+2)^2)
                    μ[tj,1,deg,i]+=Λ[tj-ti-1]/((T+2)^2)
                end
            end
        end
    end
    mom1μ = OffsetArray(rand(-1:2:1, T+2, 2, maxd, N).*0.00001, 0:T+1, 0:1, 1:maxd, 1:N)
    belief = fill(0.0, 0:T+1, N)
    ν = fill(0.0, 0:T+1, 0:T+1)
    tmpν = fill(0.0, 0:T+1, 0:T+1)
    mom1ν = fill(0.0, 0:T+1, 0:T+1)

    G=makeGraph(N,distribution);
    Neigh = Vector{Int64}[]
    Nedges=0
    for i in vertices(G)
        neigh = Int64[j for j in outneighbors(G,i)]
        append!(Neigh,[neigh])
        Nedges += length(outneighbors(G,i))
    end    
    Nedges=Int64(Nedges/2)
    x=Bool.(zeros(N,T+1))
    sample!(x,G,λp,γp)
    Observations=x[:,T+1]
    ParametricModel_1graph(N, Nedges, T, γp, λp,γi, λi, μ, mom1μ, belief, ν, tmpν, mom1ν,fr, Λ,Neigh,Observations)

end

function obs_1graph(M::ParametricModel_1graph, ti::Int64, oi::Bool)
    @unpack T, fr = M
    xT=(ti<=T)
    res = xT==oi ? (1.0 - fr) : fr
    return res
end

function makeGraph(Ngraph,degree_dist::Dirac)
    return random_regular_graph(Ngraph,degree_dist.value) |> IndexedBiDiGraph 
end
function sample!(x, G, λi, γi)
    x .= 0
    N, T = size(x)
    for i=1:N
        x[i,1] = rand() < γi
    end
    for t = 2:T
        for i = 1:N
            if x[i,t-1] == 1
                x[i,t] = 1
                continue
            end
            r = 1
            for j in inedges(G,i)
                r *= 1 - λi * x[j.src,t-1] 
            end
            x[i,t] = rand() > r
        end
    end
    x
end

function calculate_ν_1graph!(M::ParametricModel_1graph, i::Int64, j::Int64, oi::Bool)
    @unpack ν, μ, Neigh, γi, T, Λ = M
    ind_ij=findall(Neigh[i].==j)[1]
    if (length(ind_ij)> 1) error("double edges") end
    neigh=[Neigh[i][1:ind_ij-1];Neigh[i][ind_ij+1:end]]
    ν .= 0.0
    for ti=0:T+1
        ξ=obs_1graph(M, ti,oi)
        if (ξ == 0.0) continue end
        seed= ti == 0 ? γi : 1.0-γi
        phi = ti == 0 || ti == T + 1 ? 0 : 1
        m0,m1 = one(eltype(μ)), one(eltype(μ))
        for k in neigh
            ind_ki=findall(Neigh[k].==i)[1]
            if (length(ind_ki)> 1) error("double edges") end
            m0 *= μ[ti,0,ind_ki,k]
            m1 *= μ[ti,1,ind_ki,k]
        end
        for tj = 0:T+1
            ν[ti,tj] = seed*ξ*(m1*Λ[ti-tj-1] - phi*m0*Λ[ti-tj])
        end
    end
    z_ij = sum(ν)
    if (sum(ν .< 0)>0) error("negative elements in ν") end
    if (any(isnan.(ν))) error("NaN in ν") end 
    if (z_ij <= zero(eltype(ν))) error("sum(ν)=", z_ij) end    
    return z_ij
end
function calculate_ν_1graph_stab!(M::ParametricModel_1graph, i::Int64, j::Int64, m::Int64, oi::Bool)
    @unpack μ, tmpν, mom1μ, Neigh, γi, T, Λ = M
    ind_ij=findall(Neigh[i].==j)[1]
    if (length(ind_ij)> 1) error("double edges") end
    neigh=[Neigh[i][1:ind_ij-1];Neigh[i][ind_ij+1:end]]
    tmpν .= 0.0
    for ti=0:T+1
        ξ=obs_1graph(M, ti,oi)
        if (ξ == 0.0) continue end
        seed= ti == 0 ? γi : 1.0-γi
        phi = ti == 0 || ti == T + 1 ? 0 : 1
        m0,m1 = one(eltype(μ)), one(eltype(μ))
        for k in neigh
            ind_ki=findall(Neigh[k].==i)[1]
            if (length(ind_ki)> 1) error("double edges") end
            if k==m
                m0 *= mom1μ[ti,0,ind_ki,k]
                m1 *= mom1μ[ti,1,ind_ki,k]                
            else
                m0 *= μ[ti,0,ind_ki,k]
                m1 *= μ[ti,1,ind_ki,k]
            end
        end
        for tj = 0:T+1
            tmpν[ti,tj] = seed*ξ*(m1*Λ[ti-tj-1] - phi*m0*Λ[ti-tj])
        end
    end
    z_ijm = sum(tmpν)
    if (any(isnan.(tmpν))) error("NaN in tmpν") end 
    return z_ijm
end

function calculate_belief_1graph!(M::ParametricModel_1graph, i::Int64, oi::Bool)
    @unpack belief, μ, Neigh, γi, T = M
    belief[:,i].=zero(eltype(belief))
    for ti=0:T+1
        ξ=obs_1graph(M, ti,oi)
        if (ξ == 0.0) continue end
        seed= ti == 0 ? γi : 1.0-γi
        phi = ti == 0 || ti == T + 1 ? 0 : 1
        m0,m1 = one(eltype(μ)), one(eltype(μ))
        for k in Neigh[i]
            ind_ki=findall(Neigh[k].==i)[1]
            if (length(ind_ki)> 1) error("double edges") end
            m0 *= μ[ti,0,ind_ki,k]
            m1 *= μ[ti,1,ind_ki,k]            
        end
        belief[ti,i]=seed*ξ*(m1-phi*m0)
    end
    z_i = sum(belief[:,i])
    if (sum(belief[:,i].<0)>0) error("negative elements in belief at i=", i) end
    if (any(isnan.(belief[:,i]))) error("NaN in belief") end 
    if (z_i <= zero(eltype(belief[:,i]))) error("sum belief at i=", i, " is ", z_i) end    
    return z_i    
end

function update_μ_1graph!(M::ParametricModel_1graph, i::Int64, j::Int64)
    @unpack ν, μ, Λ, Neigh, T = M
    ind_ij=findall(Neigh[i].==j)[1]
    if (length(ind_ij)> 1) error("double edges") end    
    μ[:,:,ind_ij, i].=zero(eltype(μ))
    for tj=0:T+1
        m0,m1 = zero(eltype(μ)), zero(eltype(μ))
        for ti=0:T+1
            m0+=ν[ti,tj]*Λ[tj-ti]
            m1+=ν[ti,tj]*Λ[tj-ti-1]
        end
        μ[tj,0,ind_ij,i]=m0
        μ[tj,1,ind_ij,i]=m1
    end
    S=sum(@view μ[:,:,ind_ij,i])
    if (sum(μ[:,:,ind_ij,i] .<0) >0) error("negative elements in μ at i=", i, "indj=", ind_ij) end
    if (any(isnan.(μ[:,:,ind_ij,i]))) error("Nan μ") end
    if (S<=0) error("sum-zero μ") end
    return
end
function update_μ_1graph_stab!(M::ParametricModel_1graph, i::Int64, j::Int64)
    @unpack mom1ν, mom1μ, Λ, Neigh, T = M
    ind_ij=findall(Neigh[i].==j)[1]
    if (length(ind_ij)> 1) error("double edges") end    
    mom1μ[:,:,ind_ij, i].=zero(eltype(mom1μ))
    for tj=0:T+1
        m0,m1 = zero(eltype(mom1μ)), zero(eltype(mom1μ))
        for ti=0:T+1
            m0+=mom1ν[ti,tj]*Λ[tj-ti]
            m1+=mom1ν[ti,tj]*Λ[tj-ti-1]
        end
        mom1μ[tj,0,ind_ij,i]=m0
        mom1μ[tj,1,ind_ij,i]=m1
    end
    if (any(isnan.(mom1μ[:,:,ind_ij,i]))) error("Nan in mom1μ") end
    return
end

function sweep_1graph!(M::ParametricModel_1graph)
    @unpack N, Neigh, ν, belief, Observations = M
    F_i=0.0
    F_ij=0.0    
    for i in shuffle(1:N)
        oi=Observations[i]
        di=length(Neigh[i])
        for j in shuffle(Neigh[i])
            zψij=calculate_ν_1graph!(M,i,j,oi)
            ν./=zψij
            update_μ_1graph!(M,i,j)
            F_ij += log(zψij)
        end
        z_i = calculate_belief_1graph!(M, i, oi)
        belief[:,i]./=z_i
        F_i += (0.5*di-1.0)*log(z_i)
    end
    F=F_i-0.5*F_ij
    return F/N
end
function sweep_1graph_stab!(M::ParametricModel_1graph, iter::Int64, maxiter::Int64; nbstab = round(maxiter/2))
    @unpack N, Nedges, T, Neigh, ν, tmpν, mom1ν, belief, Observations = M
    F_i=0.0
    F_ij=0.0    
    Δ=0.0
    for i in shuffle(1:N)
        oi=Observations[i]
        neighs_i=Neigh[i]
        di=length(neighs_i)
        for j in shuffle(neighs_i)
            mom1ν.=0.0
            zψij=calculate_ν_1graph!(M,i,j,oi)
            F_ij += log(zψij)
            if (iter > maxiter - nbstab)
                ind_ij=findall(neighs_i.==j)[1]
                neighs_ij=[neighs_i[1:ind_ij-1];neighs_i[ind_ij+1:end]]
                for m in neighs_ij
                    zψijm=calculate_ν_1graph_stab!(M,i,j,m,oi)
                    mom1ν .+= tmpν .- ν.*(zψijm/zψij)
                end
                mom1ν./=zψij
                update_μ_1graph_stab!(M,i,j)
                Δ+=sqrt(sum(mom1ν .* mom1ν))/(T+2)
            end
            ν./=zψij
            update_μ_1graph!(M,i,j)
        end
        z_i = calculate_belief_1graph!(M, i, oi)
        belief[:,i]./=z_i
        F_i += (0.5*di-1.0)*log(z_i)
    end
    F=F_i-0.5*F_ij
    if (iter > maxiter - nbstab)
        Δ/=Nedges
        return F/N, Δ
    else
        return F/N, -1
    end
end

function pop_dynamics_1graph(M::ParametricModel_1graph; tot_iter=5)
    println("#1.iter 2.conv_crit")
    flush(stdout)
    Fold=Inf
    for iterations=1:tot_iter
        F=sweep_1graph!(M)
        println(iterations, "\t", abs(F-Fold))
        Fold=F
        flush(stdout)
    end
end
function pop_dynamics_1graph_stab(M::ParametricModel_1graph; tot_iter=5)
    println("#1.iter 2.conv_crit 3.Δ")
    flush(stdout)
    Fold=Inf
    for iterations=1:tot_iter
        F, Δ = sweep_1graph_stab!(M, iterations, tot_iter)
        println(iterations, "\t", abs(F-Fold), "\t", Δ)
        Fold=F
        flush(stdout)
    end
end
