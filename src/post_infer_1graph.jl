mutable struct ParametricModel_1graph{D,D2,M,M1,M2,O,Tλ}
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
    μ = fill(1.0 / (2*(T+2)), 0:T+1, 0:1, 1:maxd, 1:N)
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
    @show typeof(Observations)
    ParametricModel_1graph(T, γp, λp,γi, λi, μ, mom1μ, belief, ν, tmpν, mom1ν,fr, distribution, residual(distribution), Λ,Neigh,Observations)

end
