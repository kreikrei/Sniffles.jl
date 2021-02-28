# =========================================================================
#    ENVIRONMENT SETTINGS
# =========================================================================

#OPTIMIZER SETTING
const default_optimizer = Ref{Any}(nothing)
set_optimizer!(O) = default_optimizer[] = O
get_optimizer() = default_optimizer[]
reset_optimizer() = default_optimizer[] = nothing

#SLACK SURP COEFF SETTING
const slack_coeff = Ref{Any}(nothing)
const surp_coeff = Ref{Any}(nothing)
set_slack_coeff!(O) = slack_coeff[] = O
set_surp_coeff!(O) = surp_coeff[] = O
su_C() = surp_coeff[]
sl_C() = slack_coeff[]

#SILENCE SETTING
const verbosity = Ref{Any}(nothing)
silence!(O) = verbosity[] = O
silent() = verbosity[]
