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
