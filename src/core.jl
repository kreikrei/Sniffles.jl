# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

function Q(key,R)
    if isa(key,β)
        q = Vector{NamedTuple}()

        if key.q == :k
            for r in keys(R), k in keys(b().K), t in b().T
                if k == key.i #CARI K yg sama
                    push!(q,(r=r,k=k,t=t))
                end
            end
        elseif key.q == :t
            for r in keys(R), k in keys(b().K), t in b().T
                if t == key.i #CARI T yang sama
                    push!(q,(r=r,k=k,t=t))
                end
            end
        else
            for r in keys(R), k in keys(b().K), t in b().T
                if getproperty(R[r],key.q)[key.i,k,t] >= key.v
                    push!(q,(r=r,k=k,t=t))
                end
            end
        end

        return q
    elseif isa(key,Vector{β})
        q = Vector{Vector{NamedTuple}}()

        for b in key
            push!(q,Q(b,R))
        end

        if !isempty(q)
            return reduce(intersect,q)
        else
            return q
        end
    elseif isempty(key)
        q = Vector{NamedTuple}()

        for r in keys(R), k in keys(b().K), t in b().T
            push!(q,(r=r,k=k,t=t))
        end

        return q
    end
end

function f(key,R,θ)
    if !isempty(Q(key,R))
        return sum(θ[q.r,q.k,q.t] - floor(θ[q.r,q.k,q.t]) for q in Q(key,R))
    else
        return 0
    end
end

function s(key,R,θ)
    if !isempty(Q(key,R))
        return sum(θ[q.r,q.k,q.t] for q in Q(key,R))
    else
        return 0
    end
end

function master(n::node)
    mp = buildMaster(n)
    optimize!(mp)

    return mp
end

function buildMaster(n::node)
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
            for r in keys(R), k in iter_k, t in b().T
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

    # ================================
    #    BOUND GENERATOR
    # ================================
    F = Dict(1:length(n.bounds) .=> n.bounds)
    uB = filter(f -> last(f).type == :≲,F)
    lB = filter(f -> last(f).type == :≳,F)

    @constraint(mp, ρ[j = keys(uB)], sum(θ[q.r,q.k,q.t] for q in Q(F[j].B,R)) <= F[j].κ)
    @constraint(mp, σ[j = keys(lB)], sum(θ[q.r,q.k,q.t] for q in Q(F[j].B,R)) >= F[j].κ)

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

const max_component = Ref{Any}(nothing)
maxq(q::Symbol,i::Int64,k::Int64,t::Int64) = getproperty(max_component[],q)[i,k,t]

function computeMax!()
    max_val = col(
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T), #u
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T), #v
        JuMP.Containers.DenseAxisArray{Float64}(
            undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K)),b().T
        ), #l
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T), #y
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T), #z
        JuMP.Containers.DenseAxisArray{Float64}(
            undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K)),b().T
        ) #x
    )

    max_val.y .= 1
    max_val.z .= 1
    max_val.x .= 1

    for k in keys(b().K)
        max_val.u .= b().K[k].Q
        max_val.v .= b().K[k].Q
        max_val.l .= b().K[k].Q
    end

    return max_component[] = max_val
end

const column_structure = Ref{Any}(nothing)
callSubstruct() = column_structure[] #called everytime a node is processed

function colStructure!(n::node)
    sp = Model(get_optimizer())
    set_silent(sp)
    if solver_name(sp) == "Gurobi"
        set_optimizer_attribute(sp,"Cuts",2)
        set_optimizer_attribute(sp,"MIPFocus",2)
        set_optimizer_attribute(sp,"Presolve",2)
        set_optimizer_attribute(sp,"NodefileStart",0.5)
    end

    # ==========================================
    #    MODEL CONSTRUCTION (BASIC COSNTRAINTS)
    # ==========================================
    iter_k = collect(keys(b().K))
    iter_i = collect(keys(b().V))

    @variable(sp, u[iter_i, iter_k, b().T] >= 0, Int)
    @variable(sp, v[iter_i, iter_k, b().T] >= 0, Int)
    @variable(sp, l[iter_i, iter_i, iter_k, b().T] >= 0, Int)
    @variable(sp, y[iter_i, iter_k, b().T], Bin)
    @variable(sp, z[iter_i, iter_k, b().T], Bin)
    @variable(sp, x[iter_i, iter_i, iter_k, b().T], Bin)

    @constraint(sp, [k = iter_k, t = b().T],
        sum(u[i,k,t] for i in b().K[k].cover) ==
        sum(v[i,k,t] for i in b().K[k].cover) #all pickup delivered
    )

    @constraint(sp, [k = iter_k, i = b().K[k].cover, t = b().T],
        sum(x[j,i,k,t] for j in b().K[k].cover) == y[i,k,t] + z[i,k,t] #traverse in
    )

    @constraint(sp, [k = iter_k, i = b().K[k].cover, t = b().T],
        sum(x[i,j,k,t] for j in b().K[k].cover) == y[i,k,t] + z[i,k,t] #traverse out
    )

    @constraint(sp, [k = iter_k, i = b().K[k].cover, t = b().T],
        sum(l[j,i,k,t] for j in b().K[k].cover) -
        sum(l[i,j,k,t] for j in b().K[k].cover) == u[i,k,t] - v[i,k,t] #load balance
    )

    @constraint(sp, [i = iter_i, k = iter_k, t = b().T],
        u[i,k,t] <= b().K[k].Q * y[i,k,t] #u-y corr
    )

    @constraint(sp, [i = iter_i, k = iter_k, t = b().T],
        v[i,k,t] <= b().K[k].Q * z[i,k,t] #u-y corr
    )

    @constraint(sp, [i = iter_i, j = iter_i, k = iter_k, t = b().T],
        l[i,j,k,t] <= b().K[k].Q * x[i,j,k,t] #l-x corr
    )

    @constraint(sp, [k = iter_k, t = b().T],
        sum(z[i,k,t] for i in b().K[k].cover) <= 1 #only one starting point
    )

    for k in iter_k, t in b().T
        o = [i for i in iter_i if !(i in b().K[k].cover)] #sets not in cover

        if !isempty(o)
            @constraint(sp, [forbidden = o], z[forbidden,k,t] == 0)
            @constraint(sp, [forbidden = o], y[forbidden,k,t] == 0)
        end
    end

    # ================================
    #    BOUND IDENTIFICATION
    # ================================
    F = Dict(1:length(n.bounds) .=> n.bounds)
    uB = filter(f -> last(f).type == :≲,F)
    lB = filter(f -> last(f).type == :≳,F)

    @variable(sp, g[keys(uB), iter_k, b().T], Bin)
    @variable(sp, h[keys(lB), iter_k, b().T], Bin)

    q = col(u,v,l,y,z,x)

    for j in keys(uB)
        η = @variable(sp, [F[j].B, iter_k, b().T], Bin)

        @constraint(sp, [k = iter_k, t = b().T],
            g[j,k,t] >= 1 - sum(1 - η[e,k,t] for e in F[j].B)
        )
        @constraint(sp, [e = F[j].B, k = iter_k, t = b().T],
            (maxq(e.q,e.i,k,t) - e.v + 1) * η[e,k,t] >=
            (getproperty(q,e.q)[e.i,k,t] - e.v + 1)
        )
    end

    for j in keys(lB)
        η = @variable(sp, [F[j].B, iter_k, b().T], Bin)

        @constraint(sp, [e = F[j].B, k = iter_k, t = b().T],
            h[j,k,t] <= η[e,k,t]
        )
        @constraint(sp, [e = F[j].B, k = iter_k, t = b().T],
            e.v * η[e,k,t] <= getproperty(q,e.q)[e.i,k,t]
        )
    end

    optimize!(sp) #first call biar model kebuild

    return column_structure[] = sp
end

function sub(n::node,duals::dval)
    sp = buildSub(n,duals)
    optimize!(sp)

    return sp
end

function buildSub(n::node,duals::dval)
    sp = callSubstruct()
    if silent()
        set_silent(sp)
    else
        unset_silent(sp)
    end

    # ================================
    #    BOUND IDENTIFICATION
    # ================================
    F = Dict(1:length(n.bounds) .=> n.bounds)
    uB = filter(f -> last(f).type == :≲,F)
    lB = filter(f -> last(f).type == :≳,F)

    #ADD OBJECTIVE
    @objective(sp, Min,
        sum(
            sum(
                b().dist[i,j] * (
                    b().K[k].vx * sp.obj_dict[:x][i,j,k,t] +
                    b().K[k].vl * sp.obj_dict[:l][i,j,k,t]
                )
                for i in b().K[k].cover, j in b().K[k].cover
            ) +
            sum(
                b().K[k].fd * sp.obj_dict[:u][i,k,t]
                for i in b().K[k].cover
            ) +
            sum(
                b().K[k].fp * sp.obj_dict[:z][i,k,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                (sp.obj_dict[:u][i,k,t] - sp.obj_dict[:v][i,k,t]) * duals.λ[i,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                sp.obj_dict[:z][i,k,t] * duals.δ[k,i,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            duals.ϵ[k,t]
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                sp.obj_dict[:g][j,k,t] * duals.ρ[j]
                for j in keys(uB)
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                sp.obj_dict[:h][j,k,t] * duals.σ[j]
                for j in keys(lB)
            )
            for k in keys(b().K), t in b().T
        )
    )

    return sp
end

function getCols(sp::Model)
    u = value.(sp.obj_dict[:u])
    v = value.(sp.obj_dict[:v])
    l = value.(sp.obj_dict[:l])
    y = value.(sp.obj_dict[:y])
    z = value.(sp.obj_dict[:z])
    x = value.(sp.obj_dict[:x])

    return col(u,v,l,y,z,x)
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

function checkStab(mp::Model)
    s = (sum(value.(mp.obj_dict[:slack])) + sum(value.(mp.obj_dict[:surp])))

    return s
end

function colGen(n::node;maxCG::Float64,track::Bool)
    terminate = false
    iter = 0
    colStructure!(n)

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
                    println("price: $(objective_value(sp))")
                end

                if isapprox(objective_value(sp),0,atol = 1e-8) || objective_value(sp) > 0
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
        println("NODE $(n.self) FINISHED.")
    end

    return n
end

function origin(n::node)
    z = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    y = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    u = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    v = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    x = JuMP.Containers.DenseAxisArray{Float64}(
        undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K)),b().T
    )
    l = JuMP.Containers.DenseAxisArray{Float64}(
        undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K)),b().T
    )

    R = Dict(1:length(n.columns) .=> n.columns)
    mp = master(n)
    θ = mp.obj_dict[:θ]

    for i in keys(b().V), k in keys(b().K), t in b().T
        z[i,k,t] = value(sum(R[r].z[i,k,t] * θ[r,k,t] for r in keys(R)))
        y[i,k,t] = value(sum(R[r].y[i,k,t] * θ[r,k,t] for r in keys(R)))
        u[i,k,t] = value(sum(R[r].u[i,k,t] * θ[r,k,t] for r in keys(R)))
        v[i,k,t] = value(sum(R[r].v[i,k,t] * θ[r,k,t] for r in keys(R)))

        for j in keys(b().V)
            x[i,j,k,t] = value(sum(R[r].x[i,j,k,t] * θ[r,k,t] for r in keys(R)))
            l[i,j,k,t] = value(sum(R[r].l[i,j,k,t] * θ[r,k,t] for r in keys(R)))
        end
    end

    return col(
        u,v,l,
        y,z,x
    )
end
