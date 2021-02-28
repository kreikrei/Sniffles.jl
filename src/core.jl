# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]

function column(n::node)
    R = Dict{Tuple,Model}()

    for k in K(), t in T()
        sp = Model(get_optimizer())
        set_silent(sp)

        #VARIABLE DEFINITION
        @variable(sp, u[K(k).cover] >= 0, Int)
        @variable(sp, v[K(k).cover] >= 0, Int)
        @variable(sp, l[K(k).cover, K(k).cover] >= 0, Int)
        @variable(sp, y[K(k).cover], Bin)
        @variable(sp, z[K(k).cover], Bin)
        @variable(sp, x[K(k).cover, K(k).cover], Bin)

        #BASIC CONSTRAINTS
        @constraint(sp, sum(u[i] for i in K(k).cover) == sum(v[i] for i in K(k).cover))

        @constraint(sp, [i = K(k).cover], sum(x[j,i] for j in K(k).cover) == y[i] + z[i])
        @constraint(sp, [i = K(k).cover], sum(x[i,j] for j in K(k).cover) == y[i] + z[i])

        @constraint(sp, [i = K(k).cover],
            sum(l[j,i] for j in K(k).cover) - sum(l[i,j] for j in K(k).cover) == u[i] - v[i]
        ) #vehicle load balance

        @constraint(sp, [i = K(k).cover], u[i] <= K(k).Q * y[i]) #UY
        @constraint(sp, [i = K(k).cover], v[i] <= K(k).Q * z[i]) #VZ
        @constraint(sp, [i = K(k).cover, j = K(k).cover], l[i,j] <= K(k).Q * x[i,j]) #XL

        @constraint(sp, sum(z[i] for i in K(k).cover) <= 1) #one start

        #SUBPROBLEM MODIFICATIONS
        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).type == "<=" && F[j].S.k == k && F[j].S.t == t, F)
        lB = filter(f -> last(f).type == ">=" && F[j].S.k == k && F[j].S.t == t, F)

        @variable(sp, g[keys(uB)], Bin)
        @variable(sp, h[keys(lB)], Bin)

        q = col(u,v,l,y,z,x)

        for j in keys(uB)
            η = @variable(sp, [F[j].S.sequence], Bin)
            @constraint(sp, g[j] >= 1 - sum((1 - η[e]) for e in F[j].S.sequence))
            for e in F[j].S.sequence
                if e.sense == 1
                    @constraint(sp, (2 - e.v) * η[e] >= getproperty(q,e.q)[e.i] - e.v + 1)
                else #if e.sense == -1
                    @constraint(sp, e.v * η[e] >= e.v - getproperty(q,e.q)[e.i])
                end
            end
        end

        for j in keys(lB)
            η = @variable(sp, [F[j].S.sequence], Bin)
            @constraint(sp, [e = F[j].S.sequence], h[j] <= η[e])
            for e in F[j].S.sequence
                if e.sense == 1
                    @constraint(sp, e.v * η[e] <= getproperty(q,e.q)[e.i])
                else #if e.sense == -1
                    @constraint(sp, (2 - e.v) * η[e] <= 1 - getproperty(q,e.q)[e.i])
                end
            end
        end

        optimize!(sp)
        R[(k,t)] = sp
    end

    return R
end

function master(n::node)
    mp = Model(get_optimizer())
    if silent()
        set_silent(mp)
    else
        unset_silent(mp)
    end

    R = Dict(1:length(n.columns) .=> n.columns)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    @variable(mp, θ[keys(R), K(), T()] >= 0)
    @variable(mp, I[V(), vcat(first(T()) - 1, T())])
    @variable(mp, 0 <= slack[i = V(), t = T()] <= n.stab.slLim[i,t])
    @variable(mp, 0 <= surp[i = V(), t = T()] <= n.stab.suLim[i,t])

    @objective(mp, Min,
        sum(
            θ[r,k,t] * (
                sum(
                    dist(i,j) * (
                        K(k).vx * R[r][(k,t)].x[i,j] +
                        K(k).vl * R[r][(k,t)].l[i,j]
                    )
                    for i in K(k).cover, j in K(k).cover
                ) +
                sum(
                    K(k).fd * R[r][(k,t)].u[i]
                    for i in K(k).cover
                ) +
                sum(
                    K(k).fp * R[r][(k,t)].z[i]
                    for i in K(k).cover
                )
            )
            for r in keys(R), k in K(), t in T()
        ) + #column costs
        sum(
            V(i).h * I[i,t]
            for i in V(), t in T()
        ) + #inventory costs
        sum(
            n.stab.slCoeff * slack[i,t]
            for i in V(), t in T()
        ) - #stabilizer
        sum(
            n.stab.suCoeff * surp[i,t]
            for i in V(), t in T()
        ) #stabilizer
    )

    @constraint(mp, λ[i = V(), t = T()],
        I[i,t - 1] +
        sum(R[r][(k,t)].u[i] * θ[r,k,t] for r in keys(R), k in passes(i)) +
        slack[i,t] - surp[i,t] ==
        sum(R[r][(k,t)].v[i] * θ[r,k,t] for r in keys(R), k in passes(i)) +
        d(i,t) + I[i,t]
    )

    @constraint(mp, δ[k = K(), i = K(k).cover, t = T()],
        sum(R[r][(k,t)].z[i] * θ[r,k,t] for r in keys(R)) <= K(k).BP[i]
    )

    @constraint(mp, ϵ[k = K(), t = T()],
        sum(θ[r,k,t] for r in keys(R)) <= sum(K(k).BP[i] for i in K(k).cover)
    )

    @constraint(mp, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    )

    @constraint(mp, [i = V()],
        I[i,first(T())-1] == V(i).START
    )

    # ================================
    #    BOUND GENERATOR
    # ================================
    F = Dict(1:length(n.bounds) .=> n.bounds)
    uB = filter(f -> last(f).sense == "<=",F)
    lB = filter(f -> last(f).sense == ">=",F)

    @constraint(mp, ρ[j = keys(uB)], sum(θ[q.r,q.k,q.t] for q in Q(F[j].S,R)) <= F[j].κ)
    @constraint(mp, σ[j = keys(lB)], sum(θ[q.r,q.k,q.t] for q in Q(F[j].S,R)) >= F[j].κ)

    optimize!(mp)

    return mp
end

function sub(n::node,duals::dv)
    for k in K(), t in T()
        sp = callSub()[(k,t)]

        # ================================
        #    BOUND IDENTIFICATION
        # ================================
        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).type == "<=" && last(f).S.k == k && last(f).S.t == t,F)
        lB = filter(f -> last(f).type == ">=" && last(f).S.k == k && last(f).S.t == t,F)

        #ADD OBJECTIVE
        @objective(sp, Min,
            sum(
                dist(i,j) * (
                    K(k).vx * sp.obj_dict[:x][i,j] +
                    K(k).vl * sp.obj_dict[:l][i,j]
                )
                for i in K(k).cover, j in K(k).cover
            ) +
            sum(
                K(k).fd * sp.obj_dict[:u][i]
                for i in K(k).cover
            ) +
            sum(
                K(k).fp * sp.obj_dict[:z][i]
                for i in K(k).cover
            ) -
            sum(
                (sp.obj_dict[:u][i] - sp.obj_dict[:v][i]) * duals.λ[i,t]
                for i in K(k).cover
            ) -
            sum(
                sp.obj_dict[:z][i] * duals.δ[k,i,t]
                for i in K(k).cover
            ) -
            duals.ϵ[k,t] -
            sum(
                sp.obj_dict[:g][j] * duals.ρ[j]
                for j in keys(uB)
            ) -
            sum(
                sp.obj_dict[:h][j] * duals.σ[j]
                for j in keys(lB)
            )
        )

        optimize!(sp)
        if !silent()
            println("($k,$t): $(objective_value(sp))")
        end
    end

    return callSub()
end
