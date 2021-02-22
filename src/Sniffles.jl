module Sniffles

using JuMP
using XLSX
using DataFrames
using Distances
using UUIDs
using StatsBase
using Combinatorics

include("struct.jl")
include("settings.jl")
include("base.jl")
include("core.jl")
include("branch.jl")

#struct.jl
export vtx
export veh
export dt
export col
export dval
export kt
export β
export bound
export stabilizer
export node

#base.jl
export b
export extract!
export initStab
export root

#settings.jl
export set_optimizer!
export get_optimizer
export reset_optimizer
export set_slack_coeff!
export set_surp_coeff!
export su_C
export sl_C
export silence!
export silent

#core.jl
export Q
export f
export s
export master
export getDuals
export callMx!
export maxq
export colStructure!
export callSubstruct
export sub
export getCols
export colGen
export origin

#branch.jl
export separate
export Btest
export integerCheck
export createBranch

end
