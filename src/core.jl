# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]

const column_structure = Ref{Any}(nothing)
callSub() = column_structure[]

function Q(key,R)
    q = Vector{NamedTuple}()

    for seq in key.sequence
        for r in keys(R)
            if seq.sense == 1
                if getproperty(R[r][key.k,key.t],seq.q)[seq.i] >= seq.v
                    push!(q,(r=r,k=k,t=t))
                end
            else #if seq.sense == -1
                if getproperty(R[r][key.k,key.t],seq.q)[seq.i] < seq.v
                    push!(q,(r=r,k=k,t=t))
                end
            end
        end
    end

    return unique!(q)
end

s(key,R,θ) = sum(θ[q.r,q.k,q.t] for q in Q(key,R))

function column!(n::node)
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

    return column_structure = R
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
            sum(K(k).fd * sp.obj_dict[:u][i] for i in K(k).cover) +
            sum(K(k).fp * sp.obj_dict[:z][i] for i in K(k).cover) -
            sum(sp.obj_dict[:u][i] - sp.obj_dict[:v][i] * duals.λ[i,t] for i in K(k).cover)-
            sum(sp.obj_dict[:z][i] * duals.δ[k,i,t] for i in K(k).cover) -
            duals.ϵ[k,t] -
            sum(sp.obj_dict[:g][j] * duals.ρ[j] for j in keys(uB)) -
            sum(sp.obj_dict[:h][j] * duals.σ[j] for j in keys(lB))
        )

        optimize!(sp)
        if !silent()
            println("($k,$t): $(objective_value(sp))")
        end
    end

    return callSub()
end

function getDuals(mp::Model)
    λ = dual.(mp.obj_dict[:λ])
    δ = dual.(mp.obj_dict[:δ])
    ϵ = dual.(mp.obj_dict[:ϵ])

    ρ = dual.(mp.obj_dict[:ρ])
    σ = dual.(mp.obj_dict[:σ])

    return dval(λ,δ,ϵ,ρ,σ)
end

function getCols(sp)
    if isa(sp,Model)
        u = value.(sp.obj_dict[:u])
        v = value.(sp.obj_dict[:v])
        l = value.(sp.obj_dict[:l])
        y = value.(sp.obj_dict[:y])
        z = value.(sp.obj_dict[:z])
        x = value.(sp.obj_dict[:x])

        return col(u,v,l,y,z,x)
    elseif isa(sp,Dict)
        new = Dict{Tuple,col}()
        for r in keys(sp)
            new[r] = getCols(sp[r])
        end

        return new
    end
end

function colvals()
    collection = 0
    for k in K(), t in T()
        collection += objective_value(callSub()[(k,t)])
    end

    return sum(collection)
end

function updateStab!(stab::stabilizer,param::Float64)
    for i in first(stab.slLim.axes),t in last(stab.slLim.axes)
        stab.slLim[i,t] = param * stab.slLim[i,t]
        if stab.slLim[i,t] < 1
            stab.slLim[i,t] = 0
        end
    end

    for i in first(stab.suLim.axes),t in last(stab.suLim.axes)
        stab.suLim[i,t] = param * stab.suLim[i,t]
        if stab.suLim[i,t] < 1
            stab.suLim[i,t] = 0
        end
    end

    return stab
end

checkStab(mp::Model) = sum(value.(mp.obj_dict[:slack])) + sum(value.(mp.obj_dict[:surp]))

function colGen(n::node;maxCG::Float64,track::Bool)
    terminate = false
    iter = 0
    column!(n)

    while !terminate
        if iter < maxCG
            mp = master(n)

            if has_values(mp) && has_duals(mp)
                if track #print master problem obj
                    println("obj: $(objective_value(mp))")
                end

                duals = getDuals(mp)
                sp = sub(n,duals)

                if track #print subproblem price
                    println("price: $(colvals())")
                end

                if isapprox(colvals(),0,atol = 1e-8) || colvals() > 0
                    if isapprox(checkStab(mp),0,atol = 1e-8)
                        terminate = true #action
                        push!(n.status,"EVALUATED") #report
                        if track
                            println("EVALUATED")
                        end
                    else
                        updateStab!(n.stab,0.5) #action
                        push!(n.status,"STABILIZED") #report
                        if track
                            println("STABILIZED")
                        end
                    end
                else
                    push!(n.columns,getCols(sp)) #action
                    push!(n.status,"ADD_COLUMN") #report
                    if track
                        println("ADD_COLUMN")
                    end
                end

                iter += 1 #iteration update
            else
                terminate = true #action
                push!(n.status,"NO_SOLUTION")
                if track
                    println("NO_SOLUTION")
                end
            end
        else
            terminate = true #action
            push!(n.status,"EVALUATED") #report
            if track
                println("EVALUATED")
            end
        end
    end

    if n.status[end] == "NO_SOLUTION"
        println("NODE $(n.self) FAILED.")
    else
        if integerCheck(n)
            push!(n.status,"INTEGER")
            println("NODE $(n.self) INTEGER")
        else
            println("NODE $(n.self) FINISHED.")
        end
    end

    return n
end

function origin(n::node)
    z = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    y = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    u = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    v = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),K(),T())
    x = JuMP.Containers.DenseAxisArray{Float64}(
        undef,collect(V()),collect(V()),collect(K()),T()
    )
    l = JuMP.Containers.DenseAxisArray{Float64}(
        undef,collect(V()),collect(V()),collect(K()),T()
    )

    z .= 0
    y .= 0
    u .= 0
    v .= 0
    l .= 0
    x .= 0

    R = Dict(1:length(n.columns) .=> n.columns)
    mp = master(n)
    θ = mp.obj_dict[:θ]

    for k in K(), i in K(k).cover,t in T()
        z[i,k,t] = value(sum(R[r][(k,t)].z[i] * θ[r,k,t] for r in keys(R)))
        y[i,k,t] = value(sum(R[r][(k,t)].y[i] * θ[r,k,t] for r in keys(R)))
        u[i,k,t] = value(sum(R[r][(k,t)].u[i] * θ[r,k,t] for r in keys(R)))
        v[i,k,t] = value(sum(R[r][(k,t)].v[i] * θ[r,k,t] for r in keys(R)))
    end

    for k in K(), i in K(k).cover,t in T(), j in K(k).cover
        x[i,j,k,t] = value(sum(R[r][(k,t)].x[i,j] * θ[r,k,t] for r in keys(R)))
        l[i,j,k,t] = value(sum(R[r][(k,t)].l[i,j] * θ[r,k,t] for r in keys(R)))
    end

    return col(
        u,v,l,
        y,z,x
    )
end
