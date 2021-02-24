# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function separate(R,θ)
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

function Btest(R,θ)
    fract = separate(R,θ)

    for denom in 1:1000
        weight = [p/denom for p in 1:denom-1]
        println(weight)
        for p in keys(fract)
            for q in [:y,:z,:u,:v]
                for i in b().K[p[1]].cover
                    for v in ceil.(weight .* maxq(q,i,p[1]))
                        key = β(q,i,p[1],p[2],v)
                        if !issinteger(s(key,R,θ),1e-8)
                            return key
                        end
                    end
                end
            end
        end
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
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])
    seeds = Btest(R,θ)

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
