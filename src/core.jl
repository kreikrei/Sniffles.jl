# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

const max_i = Ref{Any}(nothing)
rmax() = max_i[] #buat manggil nilai max di i

function callMax!()
    val = col(
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(b().V)), #u
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(b().V)), #v
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(b().V), keys(b().V)), #empty l
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(b().V)), #y
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(b().V)), #z
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(b().V), keys(b().V)) #empty x
    )

    for i in keys(b().V)
        z = [b().K[k].BP[i] for k in keys(b().K)
            if i in b().K[k].cover
        ]
        val.z[i] = findmax(z)[1] #dapet nilai max z tiap titik

        y = [sum(b().K[k].BP[p] for p in b().K[k].cover) for k in keys(b().K)
            if i in b().K[k].cover
        ]
        val.y[i] = findmax(y)[1] #dapet nilai max z tiap titik

        v = [b().K[k].BP[i] * b().K[k].Q for k in keys(b().K)
            if i in b().K[k].cover
        ]
        val.v[i] = findmax(v)[1]

        u = [sum(b().K[k].BP[p] for p in b().K[k].cover) * b().K[k].Q
            for k in keys(b().K) if i in b().K[k].cover
        ]
        val.u[i] = findmax(u)[1]
    end

    return max_i[] = val
end

function Q(key,R::Dict)
    if isa(key,β)
        q = Vector{NamedTuple}()

        for r in keys(R), k in keys(b().K), t in b().T
            if getproperty(R[r],key.q)[key.i,k,t] >= key.v
                push!(q,(r=r,k=k,t=t))
            end
        end

        return q
    elseif isa(key,Vector{β})
        if !isempty(key)
            q = Vector{Vector{NamedTuple}}()

            for b in key
                push!(q,Q(b,R))
            end

            if !isempty(q)
                #ngembaliin semua rkt
                return reduce(intersect,q)
            else
                #ngembaliin vector kosong
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
end

function f(B,R,θ)
    if !isempty(Q(B,R))
        return sum(θ[q.r,q.k,q.t] - floor(θ[q.r,q.k,q.t]) for q in Q(B,R))
    else
        return 0
    end
end

function tot(B,R,θ)
    if !isempty(Q(B,R))
        return sum(θ[q.r,q.k,q.t] for q in Q(B,R))
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

    # ================================
    #    BOUND GENERATOR
    # ================================
    @constraint(mp, ρ[b = keys(≲)],
        sum(θ[q.r,q.k,q.t] for q in Q(B[b].B,R)) <= B[b].κ
    )

    @constraint(mp, σ[b = keys(≳)],
        sum(θ[q.r,q.k,q.t] for q in Q(B[b].B,R)) >= B[b].κ
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

function sub(n::node,duals::dval)
    sp = buildSub(n,duals)
    optimize!(sp)

    return sp
end

const subproblem_model = Ref{Any}(nothing)
colStructure() = subproblem_model[] #pemanggil model subproblem

function callSub!()
    sp = Model(get_optimizer())
    if solver_name(sp) == "Gurobi"
        set_optimizer_attribute(sp,"Cuts",2)
        set_optimizer_attribute(sp,"MIPFocus",2)
        set_optimizer_attribute(sp,"Presolve",2)
        set_optimizer_attribute(sp,"Cutoff",0)
        set_optimizer_attribute(sp,"NodefileStart",0.5)
    end

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    iter_k = collect(keys(b().K))
    iter_i = collect(keys(b().V))

    q = col(
        @variable(sp, u[iter_i, iter_k, b().T] >= 0, Int),
        @variable(sp, v[iter_i, iter_k, b().T] >= 0, Int),
        @variable(sp, l[iter_i, iter_i, iter_k, b().T] >= 0, Int),
        @variable(sp, y[iter_i, iter_k, b().T], Bin),
        @variable(sp, z[iter_i, iter_k, b().T], Bin),
        @variable(sp, x[iter_i, iter_i, iter_k, b().T], Bin)
    )

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

    for k in keys(b().K), t in b().T
        n = [i for i in iter_i if !(i in b().K[k].cover)] #sets not in cover

        if !isempty(n)
            @constraint(sp, [forbidden = n], z[forbidden,k,t] == 0)
            @constraint(sp, [forbidden = n], y[forbidden,k,t] == 0)
        end
    end

    optimize!(sp) #first call biar model kebuild

    return subproblem_model[] = sp
end

function buildSub(n::node,duals::dval)
    sp = colStructure()

    if silent()
        set_silent(sp)
    end

    #BRANCHING PREP
    B = Dict(1:length(n.bounds) .=> n.bounds)
    ≲ = filter(b -> last(b).type == "≲",B)
    ≳ = filter(b -> last(b).type == "≳",B)

    g = @variable(sp, [keys(≲), keys(b().K), b().T], Bin) #anonymous variable addition
    h = @variable(sp, [keys(≳), keys(b().K), b().T], Bin) #anonymous variable addition

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
                g[j,k,t] * duals.ρ[j]
                for j in keys(≲)
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                h[j,k,t] * duals.σ[j]
                for j in keys(≳)
            )
            for k in keys(b().K), t in b().T
        )
    )

    #BRANCHING CONSTRAINTS FOR SUBPROBLEM
    for j in keys(≲), k in keys(b().K), t in b().T
        η = @variable(sp, [B[j].B], Bin)

        @constraint(sp, -g[j,k,t] >= 1 - sum(1 - η[e] for e in B[j].B))
        for e in B[j].B
            @constraint(sp, (getproperty(rmax(),e.q)[e.i] - e.v + 1) * η[e] >=
                getproperty(q,e.q)[e.i,k,t] - e.v + 1
            )
        end
    end

    for j in keys(≳), k in keys(b().K), t in b().T
        η = @variable(sp, [B[j].B], Bin)

        @constraint(sp, [e = B[j].B], h[j,k,t] <= η[e])
        for e in B[j].B
            @constraint(sp, e.v * η[e] <= getproperty(q,e.q)[e.i,k,t])
        end
    end

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
    mem = 0

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
