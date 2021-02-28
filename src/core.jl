# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

function master(n::node)

    return mp
end

function sub(n::node,duals::dv)

    return sp
end

function column(n::node)
    G = Dict{Tuple,Model}()

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
        uB = filter(f -> last(f).type == ">=", F)
        lB = filter(f -> last(f).type == "=<", F)

        @variable(sp, g[keys(uB)], Bin)
        @variable(sp, h[keys(lB)], Bin)

        q = col(u,v,l,y,z,x)

        for j in keys(uB)

        end

        for j in keys(lB)

        end

        optimize!(sp)
        G[(k,t)] = sp
    end

    return G
end
