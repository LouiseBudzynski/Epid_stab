using Random

mutable struct ParametricModel_1graph{D,D2,M,M1,M2,O,Tλ}
    N::Int
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
    distribution::D
    residual::D2
    Λ::O
    Neigh::Vector{Vector{Int64}}
    Observations::BitVector
end

function ParametricModel_1graph(; N, T, γp, λp, γi=γp, λi=λp, fr=0.0, dilution=0.0, distribution, maxd)
    μ = fill(1.0 / ((T+2)^2), 0:T+1, 0:1, maxd, N)
    mom1μ = OffsetArray(rand(-1:2:1, T+2, 2, maxd, N).*0.00001, 0:T+1, 0:1, 1:maxd, 1:N)
    belief = fill(0.0, 0:T+1, N)
    ν = fill(0.0, 0:T+1, 0:T+1)
    tmpν = fill(0.0, 0:T+1, 0:T+1)
    mom1ν = fill(0.0, 0:T+1, 0:T+1)
    Λ = OffsetArray([t <= 0 ? 1.0 : (1-λi)^t for t = -T-2:T+1], -T-2:T+1)

    G=makeGraph(N,distribution);
    Neigh = Vector{Int64}[]
    for i in vertices(G)
        neigh = Int64[j for j in outneighbors(G,i)]
        append!(Neigh,[neigh])
    end    
    x=Bool.(zeros(N,T))
    sample!(x,G,λp,γp)
    Observations=x[:,T]
    ParametricModel_1graph(N, T, γp, λp,γi, λi, μ, mom1μ, belief, ν, tmpν, mom1ν,fr, distribution, residual(distribution), Λ,Neigh,Observations)

end

function obs_1graph(ti::Int64, oi::Bool)
    @unpack T, fr = M
    xT=(ti<=T)
    res = xT==oi ? (1.0 - fr) : fr
    return res
end

function calculate_ν_1graph!(M::ParametricModel_1graph, i::Int64, j::Int64, oi::Bool)
    @unpack ν, μ, Neigh, γi, Λ = M
    ind_ij=findall(Neigh[i].==j)[1]
    if (length(ind_ij)> 1) error("double edges") end
    neigh=[Neigh[i][1:ind_ij-1];Neigh[i][ind_ij+1:end]]
    ν .= 0.0
    for ti=0:T+1
        ξ=obs_1graph(ti,oi)
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

function calculate_ν_1graph_stab!()
end

function calculate_belief_1graph!(M::ParametricModel_1graph, i::Int64, oi::Bool)
    @unpack belief, μ, Neigh, γi= M
    belief[:,i].=zero(eltype(belief))
    for ti=0:T+1
        ξ=obs_1graph(ti,oi)
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
    @unpack ν, μ, Λ, Neigh = M
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

function update_momμ_1graph!()
end

function sweep_1graph!(M::ParametricModel_1graph)
    @unpack N, Neigh, ν, Observations = M
    F=0.0
    for i in 1:N #shuffle(1:N)
        oi=Observations[i]
        di=length(Neigh[i])
        for j in Neigh[i] #shuffle(Neigh[i])
            zψij=calculate_ν_1graph!(M,i,j,oi)
            ν./=zψij
            update_μ_1graph!(M,i,j)
            F += -0.5*log(zψij)
        end
        z_i = calculate_belief_1graph!(M, i, oi)
        F += (0.5*di-1.0)*log(z_i)
    end
    return F/N
end
function sweep_1graph_stab!(M::ParametricModel_1graph)

end

function pop_dynamics_1graph(M::ParametricModel_1graph; tot_iter=5)
    F=0.0
    F_window=zeros(10)
    converged=false
    err=-1
    println("#1.iter 2.F")
    for iterations=1:tot_iter
        F=sweep_1graph!(M)
        println(iterations, "\t", F)
    end
end
function pop_dynamics_1graph_stab()
end