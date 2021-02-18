# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

function master(n::node)
    mp = buildMaster(n)
    optimize!(mp)

    return mp
end

function buildMaster(n::node)
    mp = Model(get_default_optimizer())
    if silent()
        set_silent(mp)
    end

    R = Dict(1:length(n.columns) .=> n.columns)
    B = Dict(1:length(n.bounds) .=> n.bounds)
    ≲ = filter(b -> last(b).type == "≲",B)
    ≳ = filter(b -> last(b).type == "≳",B)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    iter_k = collect(keys(b().K))
    iter_i = collect(keys(b().V))

    @variable(mp, θ[keys(R), iter_k, b().T] >= 0)
    @variable(mp, I[iter_i, vcat(first(b().T) - 1, b().T)])
    @variable(mp, 0 <= slack[i = iter_i, t = b().T] <= n.stab.slLim[i,t])
    @variable(mp, 0 <= surp[i = iter_i, t = b().T] <= n.stab.suLim[i,t])

    @objective(mp, Min,
        sum(
            θ[r,k,t] * (
                sum(
                    b().dist[i,j] * (
                        b().K[k].vx * R[r].x[i,j,k,t] +
                        b().K[k].vl * R[r].l[i,j,k,t]
                    )
                    for i in b().K[k].cover, j in b().K[k].cover
                ) +
                sum(
                    b().K[k].fd * R[r].u[i,k,t]
                    for i in b().K[k].cover
                ) +
                sum(
                    b().K[k].fp * R[r].z[i,k,t]
                    for i in b().K[k].cover
                )
            )
            for r in keys(R), k in keys(b().K), t in b().T
        ) + #column costs
        sum(
            b().V[i].h * I[i,t]
            for i in iter_i, t in b().T
        ) + #inventory costs
        sum(
            n.stab.slCoeff * slack[i,t]
            for i in iter_i, t in b().T
        ) - #stabilizer
        sum(
            n.stab.suCoeff * surp[i,t]
            for i in iter_i, t in b().T
        ) #stabilizer
    )

    @constraint(mp, λ[i = iter_i, t = b().T],
        I[i,t - 1] +
        sum(R[r].u[i,k,t] * θ[r,k,t] for r in keys(R), k in iter_k) +
        slack[i,t] - surp[i,t] ==
        sum(R[r].v[i,k,t] * θ[r,k,t] for r in keys(R), k in iter_k) +
        b().d[i,t] + I[i,t]
    )

    @constraint(mp, δ[k = iter_k, i = b().K[k].cover, t = b().T],
        sum(R[r].z[i,k,t] * θ[r,k,t] for r in keys(R)) <= b().K[k].BP[i]
    )

    @constraint(mp, ϵ[k = iter_k, t = b().T],
        sum(θ[r,k,t] for r in keys(R)) <= sum(b().K[k].BP[i] for i in b().K[k].cover)
    )

    @constraint(mp, [i = keys(b().V), t = b().T],
        b().V[i].MIN <= I[i,t] <= b().V[i].MAX
    )

    @constraint(mp, [i = keys(b().V)],
        I[i,first(b().T)-1] == b().V[i].START
    )

    return mp
end

function getDuals(mp::Model)
    λ = dual.(mp.obj_dict[:λ])
    δ = dual.(mp.obj_dict[:δ])
    ϵ = dual.(mp.obj_dict[:ϵ])
    ρ = dual.(mp.obj_dict[:ρ])
    σ = dual.(mp.obj_dict[:σ])

    return dval(λ,δ,ϵ,ρ,σ)
end
