# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

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

        @constraint(sp, [i = K(k).cover], u[i] <= b().K[k].Q * y[i]) #UY
        @constraint(sp, [i = K(k).cover], v[i] <= b().K[k].Q * z[i]) #VZ
        @constraint(sp, [i = K(k).cover, j = K(k).cover], l[i,j] <= b().K[k].Q * x[i,j]) #XL

        @constraint(sp, sum(z[i] for i in K(k).cover) <= 1) #one start

        #SUBPROBLEM MODIFICATIONS
        F = Dict(1:length(n.bounds) .=> n.bounds)
        uB = filter(f -> last(f).type == ">=" && F[j].S.k == k && F[j].S.t == t, F)
        lB = filter(f -> last(f).type == "=<" && F[j].S.k == k && F[j].S.t == t, F)

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

    return mp
end

function sub(n::node,duals::dv)
    
    return sp
end
