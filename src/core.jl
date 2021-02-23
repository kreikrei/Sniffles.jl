# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

function Q(key,R)
    if isa(key,β)
        q = Vector{NamedTuple}()

        if key.q == :k
            for r in keys(R), k in keys(b().K), t in b().T
                if k == key.i #CARI VEHICLE
                    push!(q,(r=r,k=k,t=t))
                end
            end
        elseif key.q == :t
            for r in keys(R), k in keys(b().K), t in b().T
                if t == key.i #CARI PERIOD
                    push!(q,(r=r,k=k,t=t))
                end
            end
        else
            for r in keys(R), k in keys(b().K), t in b().T
                if getproperty(R[r][(k,t)],key.q)[key.i] >= key.v
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
    else
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
                        b().K[k].vx * R[r][(k,t)].x[i,j] +
                        b().K[k].vl * R[r][(k,t)].l[i,j]
                    )
                    for i in b().K[k].cover, j in b().K[k].cover
                ) +
                sum(
                    b().K[k].fd * R[r][(k,t)].u[i]
                    for i in b().K[k].cover
                ) +
                sum(
                    b().K[k].fp * R[r][(k,t)].z[i]
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
        sum(R[r][(k,t)].u[i] * θ[r,k,t] for r in keys(R), k in iter_k) +
        slack[i,t] - surp[i,t] ==
        sum(R[r][(k,t)].v[i] * θ[r,k,t] for r in keys(R), k in iter_k) +
        b().d[i,t] + I[i,t]
    )

    @constraint(mp, δ[k = iter_k, i = b().K[k].cover, t = b().T],
        sum(R[r][(k,t)].z[i] * θ[r,k,t] for r in keys(R)) <= b().K[k].BP[i]
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

    optimize!(mp)

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
maxq(q::Symbol,i::Int64,k::Int64) = getproperty(max_component[],q)[i,k]

function callMx!()
    max_val = col(
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K)), #u
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K)), #v
        JuMP.Containers.DenseAxisArray{Float64}(
            undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K))
        ), #l
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K)), #y
        JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K)), #z
        JuMP.Containers.DenseAxisArray{Float64}(
            undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K))
        ) #x
    )

    max_val.y .= 0
    max_val.z .= 0
    max_val.x .= 0
    max_val.u .= 0
    max_val.v .= 0
    max_val.l .= 0

    for k in keys(b().K), i in b().K[k].cover
        max_val.y[i,k] = 1
        max_val.z[i,k] = 1
        max_val.u[i,k] = b().K[k].Q
        max_val.v[i,k] = b().K[k].Q
        for j in b().K[k].cover
            max_val.x[i,j,k] = 1
            max_val.l[i,j,k] = b().K[k].Q
        end
    end

    return max_component[] = max_val
end

const column_structure = Ref{Any}(nothing)
callSubstruct() = column_structure[] #called everytime a node is processed

function colStructure!(n::node)
    # ==========================================
    #    MODEL CONSTRUCTION (BASIC COSNTRAINTS)
    # ==========================================
    R = Dict{Tuple,Model}()

    iter_k = collect(keys(b().K))
    iter_i = collect(keys(b().V))

    @inbounds for k in iter_k, t in b().T
        sp = Model(get_optimizer())
        set_silent(sp)
        if solver_name(sp) == "Gurobi"
            set_optimizer_attribute(sp,"Cuts",2)
            set_optimizer_attribute(sp,"MIPFocus",2)
            set_optimizer_attribute(sp,"Threads",1)
            set_optimizer_attribute(sp,"NodefileStart",0.5)
            set_optimizer_attribute(sp, "NumericFocus",3)
        end

        @variable(sp, u[iter_i] >= 0, Int)
        @variable(sp, v[iter_i] >= 0, Int)
        @variable(sp, l[iter_i, iter_i] >= 0, Int)
        @variable(sp, y[iter_i], Bin)
        @variable(sp, z[iter_i], Bin)
        @variable(sp, x[iter_i, iter_i], Bin)

        @constraint(sp,
            sum(u[i] for i in b().K[k].cover) ==
            sum(v[i] for i in b().K[k].cover) #all pickup delivered
        )

        @constraint(sp, [i = b().K[k].cover],
            sum(x[j,i] for j in b().K[k].cover) == y[i] + z[i] #traverse in
        )

        @constraint(sp, [i = b().K[k].cover],
            sum(x[i,j] for j in b().K[k].cover) == y[i] + z[i] #traverse out
        )

        @constraint(sp, [i = b().K[k].cover],
            sum(l[j,i] for j in b().K[k].cover) -
            sum(l[i,j] for j in b().K[k].cover) == u[i] - v[i] #load balance
        )

        @constraint(sp, [i = iter_i],
            u[i] <= b().K[k].Q * y[i] #u-y corr
        )

        @constraint(sp, [i = iter_i],
            v[i] <= b().K[k].Q * z[i] #u-y corr
        )

        @constraint(sp, [i = iter_i, j = iter_i],
            l[i,j] <= b().K[k].Q * x[i,j] #l-x corr
        )

        @constraint(sp,
            sum(z[i] for i in b().K[k].cover) <= 1 #only one starting point
        )

        o = [i for i in iter_i if !(i in b().K[k].cover)] #sets not in cover
        if !isempty(o)
            @constraint(sp, [forbidden = o], z[forbidden] == 0)
            @constraint(sp, [forbidden = o], y[forbidden] == 0)
        end

        # ================================
        #    BOUND IDENTIFICATION
        # ================================
        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).type == :≲,F)
        lB = filter(f -> last(f).type == :≳,F)

        @variable(sp, g[keys(uB)], Bin)
        @variable(sp, h[keys(lB)], Bin)

        q = col(u,v,l,y,z,x)

        for j in keys(uB)
            inclusion = filter(f -> f.q in [:k,:t], F[j].B)
            included = true
            if !isempty(inclusion)
                for f in inclusion
                    if (f.q == :k && f.i != k) || (f.q == :t && f.i != t)
                        included = false
                    end
                end
            end

            #IF K AND T CHECKS OUT EXECUTE BOUND
            if included
                println("upper bound $j included on (k=$k,t=$t).")
                unit = filter(f -> !(f.q in [:k,:t]), F[j].B)

                η = @variable(sp, [unit], Bin)
                @constraint(sp, g[j] >= 1 - sum(1 - η[e] for e in unit))
                for e in unit
                    @constraint(sp,
                        (maxq(e.q,e.i,k) - e.v + 1) * η[e] >=
                        (getproperty(q,e.q)[e.i] - e.v + 1)
                    )
                end
            else
                @constraint(sp, g[j] == 0)
            end
        end

        for j in keys(lB)
            inclusion = filter(f -> f.q in [:k,:t], F[j].B)
            included = true
            if !isempty(inclusion)
                for f in inclusion
                    if (f.q == :k && f.i != k) || (f.q == :t && f.i != t)
                        included = false
                    end
                end
            end

            #IF K AND T CHECKS OUT EXECUTE BOUND
            if included
                println("lower bound $j included on (k=$k,t=$t).")
                unit = filter(f -> !(f.q in [:k,:t]), F[j].B)

                η = @variable(sp, [unit], Bin)
                @constraint(sp, [e = unit], h[j] <= η[e])
                for e in unit
                    @constraint(
                        sp, e.v * η[e] <= getproperty(q,e.q)[e.i]
                    )
                end
            else
                @constraint(sp, h[j] == 0)
            end
        end

        optimize!(sp) #first call biar model kebuild

        R[(k,t)] = sp
    end

    return column_structure[] = R
end

function sub(n::node,duals::dval)
    # ==========================================
    #    ADD OBJECTIVE AND SOLVE (for each kt)
    # ==========================================
    iter_k = collect(keys(b().K))
    iter_i = collect(keys(b().V))

    @inbounds for k in iter_k, t in b().T
        sp = callSubstruct()[(k,t)]
        set_silent(sp)

        # ================================
        #    BOUND IDENTIFICATION
        # ================================
        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).type == :≲,F)
        lB = filter(f -> last(f).type == :≳,F)

        #ADD OBJECTIVE
        @objective(sp, Min,
            sum(
                b().dist[i,j] * (
                    b().K[k].vx * sp.obj_dict[:x][i,j] +
                    b().K[k].vl * sp.obj_dict[:l][i,j]
                )
                for i in b().K[k].cover, j in b().K[k].cover
            ) +
            sum(
                b().K[k].fd * sp.obj_dict[:u][i]
                for i in b().K[k].cover
            ) +
            sum(
                b().K[k].fp * sp.obj_dict[:z][i]
                for i in b().K[k].cover
            ) -
            sum(
                (sp.obj_dict[:u][i] - sp.obj_dict[:v][i]) * duals.λ[i,t]
                for i in b().K[k].cover
            ) -
            sum(
                sp.obj_dict[:z][i] * duals.δ[k,i,t]
                for i in b().K[k].cover
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

    return callSubstruct()
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
        @inbounds for r in keys(sp)
            new[r] = getCols(sp[r])
        end

        return new
    end
end

function colvals()
    collection = 0
    @inbounds for k in keys(b().K), t in b().T
        collection += objective_value(callSubstruct()[(k,t)])
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

function checkStab(mp::Model)
    s = (sum(value.(mp.obj_dict[:slack])) + sum(value.(mp.obj_dict[:surp])))

    return s
end

function colGen(n::node;maxCG::Float64,track::Bool)
    terminate = false
    iter = 0
    colStructure!(n) #defines the structure of the node

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
            println("NODE $(n.id) INTEGER")
        else
            println("NODE $(n.self) FINISHED.")
        end
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
