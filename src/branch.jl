# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function separate(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    res = Vector{NamedTuple}() #SEMUA RKT yang fractional
    for r in θ.axes[1], k in θ.axes[2], t in θ.axes[3]
        if !isinteger(θ[r,k,t])
            push!(res,(r=r,k=k,t=t))
        end
    end

    group = Dict{Tuple,Vector{Int64}}()
    for k in [p.k for p in res], t in [p.t for p in res if p.k == k]
        group[(k,t)] = [p.r for p in res if p.k == k && p.t == t]
    end

    return group
end

function Btest(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    fract = separate(n)
    for id in keys(fract)
        stack = Vector{Vector{β}}()
        key = β[β(:k,id[1]),β(:t,id[2])]

        #GENERATE STACK OF q and i for (k,t)
        for q in [:u,:v,:y,:z], i in b().K[id[1]].cover
            store = Vector{Int64}()
            for r in fract[id]
                val = Int64(getproperty(R[r][id],q)[i])
                push!(store,val)
            end
            unique!(store) #distinct values of qi

            test = Vector{Int64}()
            for w in 1:length(store)-1
                comp = ceil((store[w] + store[w+1]) / 2)
                push!(test,comp)
            end #summarize the values of qi

            if !isempty(test)
                for v in test
                    push!(stack,vcat(key,[β(q,i,v)]))
                end
            end
        end

        #combination and testing of stack
        for pair in [[:y,:z],[:u,:v]]
            cardinality = 1
            while cardinality <= floor(log2(f([],R,θ))) + 1
                q_combo = collect(combinations(pair,cardinality))
                i_combo = collect(combinations(b().K[id[1]].cover,cardinality))

                for q in q_combo, i in i_combo
                    raw = filter(p -> last(p).i in i && last(p).q in q,stack)
                    if !isempty(raw)
                        pure = reduce(union,raw)
                        val = s(pure,R,θ)
                        if !Sniffles.issinteger(val,1e-10)
                            return pure
                        end
                    end
                end
                cardinality += 1
            end
        end
        #println("k: $(s[1]),t: $(s[2])")
        #println(stack)
    end
end

issinteger(val,tol) = abs(round(val) - val) < tol

function integerCheck(n::node)
    integer = true
    θ = value.(master(n).obj_dict[:θ])

    for val in θ
        if !issinteger(val,1e-8)
            integer = false
            break
        end
    end

    return integer
end

function createBranch(n::node)
    branches = Vector{node}()
    seeds = Btest(n)
    if isempty(seeds)
        println("no B to be found")
    end
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    for br in [:≳,:≲]
        push!(branches,
            node(
                n.self, #parent
                uuid1(), #self

                vcat(n.bounds, #bounds
                    bound(
                        seeds,br,
                        if br == :≲
                            floor(s(seeds,R,θ))
                        else
                            ceil(s(seeds,R,θ))
                        end
                    )
                ),
                deepcopy(n.columns), #columns
                initStab(), #stabilizer
                ["UNVISITED"]
            )
        )
    end

    return branches
end

function leaf(n::node,upperBound::Float64,maxiter::Float64)
    println("leaf search initiated.")
    #INITIALIZE DFS
    terminate = false
    stack = Vector{node}()
    visited = Vector{node}()
    iter = 0

    #PUSH ROOT TO STACK
    push!(stack,n)

    #ITERATION
    while !(terminate || isempty(stack) || iter > maxiter)
        #POP FROM STACK
        u = pop!(stack)

        #CHECK AND PROCESS
        if u.status[end] == "UNVISITED"
            #PROCESS THE NODE
            colGen(u;track=false,maxCG=Inf)

            #NODE STATUSES
            if u.status[end] == "INTEGER"
                obj = objective_value(master(u))

                #CHANGE UPPERBOUND
                if obj < upperBound
                    upperBound = obj
                    push!(u.status,"NEW_BOUND")
                end

                #ALWAYS TERMINATE IF INTEGER NODE IS FOUND
                terminate = true
            elseif u.status[end] == "EVALUATED"
                obj = objective_value(master(u))

                #BRANCH IF NOT ABOVE uB
                if obj <= upperBound
                    append!(stack,createBranch(u))
                end
            end

            #SEND TO VISITED
            push!(visited,u)
            iter += 1
        end
    end

    return (stack=stack,visited=visited,upperBound=upperBound)
end


function traverse(max::Float64)
    #VALUE STORE
    stack = Vector{node}()
    visited = Vector{node}()
    upperBound = Inf

    #FLOW CONTROL
    terminate = false

    #CREATE ROOT
    newroot = root()

    #PUSH ROOT TO STACK
    push!(stack,newroot)

    #ITERATION
    while !(terminate || isempty(stack))
        #POP FIRST ELEMENT OF STACK
        n = pop!(stack)

        #FIND LEAVE (INTEGER NODE) FROM NODE
        results = leaf(n,upperBound,Inf)

        #USE THE RESULTS
        append!(stack,results.stack)
        append!(visited,results.visited)
        if results.upperBound < upperBound
            upperBound = results.upperBound
        end

        #PRINT CURRENT CONDITIONS
        println("size of stack: $(length(stack))")
        println("size of visited: $(length(visited))")
        println("current upperBound: $upperBound")

        #MAX ITER REACHED
        if length(visited) >= max
            terminate = true
        end
    end

    return stack,visited,upperBound
end
