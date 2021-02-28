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

        @variable(sp, u[K(k).cover] >= 0, Int)
        @variable(sp, v[K(k).cover] >= 0, Int)
        @variable(sp, l[K(k).cover, K(k).cover] >= 0, Int)
        @variable(sp, y[K(k).cover], Bin)
        @variable(sp, z[K(k).cover], Bin)
        @variable(sp, x[K(k).cover, K(k).cover], Bin)

        @constraint(sp, sum(u[i] for i in K(k).cover) == sum(v[i] for i in K(k).cover))

        @constraint(sp, [i = K(k).cover], sum(x[j,i] for j in K(k).cover) == y[i] + z[i])
        @constraint(sp, [i = K(k).cover], sum(x[i,j] for j in K(k).cover) == y[i] + z[i])

        @constraint(sp, [i = K(k).cover],
            sum(l[j,i] for j in K(k).cover) - sum(l[i,j] for j in K(k).cover) == u[i] - v[i]
        )

        @constraint(sp, [i = K(k).cover], u[i] <= b().K[k].Q * y[i])
        @constraint(sp, [i = K(k).cover], v[i] <= b().K[k].Q * z[i])
        @constraint(sp, [i = K(k).cover, j = K(k).cover], l[i,j] <= b().K[k].Q * x[i,j])

        @constraint(sp, sum(z[i] for i in K(k).cover) <= 1)

        optimize!(sp)
        G[(k,t)] = sp
    end

    return G
end
