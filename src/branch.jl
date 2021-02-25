# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function separate(bounds,R,θ)
    #FIND x TO BRANCH
    for k in keys(b().K), t in b().T
        for i in b().K[k].cover, j in b().K[k].cover
            key = β[β(:x,i,j,k,t,1)]
            val = s(key,R,θ)
            if !issinteger(val,1e-8)
                return key
            end
        end
    end

    #FIND l TO BRANCH
    for k in keys(b().K), t in b().T
        for i in b().K[k].cover, j in b().K[k].cover
            for v in LV(i,j,k,t,R,θ)
                key = β[β(:l,i,j,k,t,v)]
                val = s(key,R,θ)
                if !issinteger(val,1e-8)
                    return key
                end
            end
        end
    end
end

function LV(i::Int64,j::Int64,k::Int64,t::Int64,R,θ)
    res = Vector{Int64}()
    for r in keys(R)
        if !isinteger(θ[r,k,t])
            push!(res,r)
        end
    end

    store = Vector{Int64}() #collect all values contained in rkt
    for r in res
        val = getproperty(R[r][(k,t)],:l)[i,j]
        push!(store,val)
    end

    unique!(store)

    test = Vector{Int64}() #summarize the values, we average 1-by-1
    for w in 1:length(store)-1
        comp = ceil((store[w] + store[w+1]) / 2)
        push!(test,comp)
    end

    return test
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
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    seeds = separate(n.bounds,R,θ)
    println("branch on $seeds: $(s(seeds,R,θ))")

    for br in [:≲,:≳]
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
                println("found integer.")
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
