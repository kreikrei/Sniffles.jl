# =========================================================================
#    BASIC STRUCTURES
# =========================================================================

struct vtx
    #IDENTIFIERS
    name::String
    type::String

    #VALS
    x::Float64 #xcoor
    y::Float64 #ycoor
    MAX::Float64 #max inventory level
    MIN::Float64 #min inventory level
    START::Float64 #starting inventory level

    #COSTS
    h::Float64 #inv cost per unit
end

struct veh
    #IDENTIFIERS
    name::String
    type::String

    #CHAR
    cover::Vector{Int64}
    BP::Dict{Int64,Int64}
    Q::Int64

    #COSTS
    vx::Float64
    vl::Float64
    fp::Float64
    fd::Float64
end

struct dt
    V::Dict{Int64,vtx}
    dist::JuMP.Containers.DenseAxisArray
    K::Dict{Int64,veh}
    T::Vector{Int64}
    d::JuMP.Containers.DenseAxisArray
end

struct col
    #quantity
    u::JuMP.Containers.DenseAxisArray
    v::JuMP.Containers.DenseAxisArray
    l::JuMP.Containers.DenseAxisArray

    #decision
    y::JuMP.Containers.DenseAxisArray
    z::JuMP.Containers.DenseAxisArray
    x::JuMP.Containers.DenseAxisArray
end

struct dval
    #three master constraints
    λ::JuMP.Containers.DenseAxisArray
    δ::JuMP.Containers.SparseAxisArray
    ϵ::JuMP.Containers.DenseAxisArray

    #bounding constraints
    ρ::JuMP.Containers.DenseAxisArray
    σ::JuMP.Containers.DenseAxisArray
end

struct β{T<:Any,S<:Any}
    q::T
    i::S
    k::S
    t::S
    v::S
end

struct bound
    β::β #component bound set
    type::Symbol #≳ or ≲
    κ::Int64 #value of aggregation
end

struct stabilizer
    slCoeff::Float64
    suCoeff::Float64
    slLim::JuMP.Containers.DenseAxisArray
    suLim::JuMP.Containers.DenseAxisArray
end

struct node
    #IDENTIFIER
    parent::UUID
    self::UUID

    #Dynamic SET
    bounds::Vector{bound}
    columns::Vector{Dict}

    #SUPPORT
    stab::stabilizer

    #STATUS
    status::Vector{String}
end

function node()
    return node(
        uuid1(), uuid1(),
        Vector{bound}(),Vector{col}(),
        initStab(),
        ["UNVISITED"]
    )
end
