module Sniffles

using JuMP
using XLSX
using DataFrames
using Distances
using UUIDs

include("struct.jl")
include("settings.jl")
include("base.jl")

#struct.jl
export vtx
export veh
export dt
export col
export dval
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

end
