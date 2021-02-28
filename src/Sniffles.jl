module Sniffles

using JuMP
using XLSX
using DataFrames
using Distances
using UUIDs
using Combinatorics

include("struct.jl")
include("settings.jl")
include("base.jl")
include("core.jl")
include("branch.jl")

end
