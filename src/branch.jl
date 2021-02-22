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

    #println(res)

    minBoundstack = Vector{Vector{β}}() #HITUNG V tiap qi di rkt
    maxBoundstack = Vector{Vector{β}}()
    for q in [:u,:v,:y,:z]
        for i in keys(b().V), k in keys(b().K), t in b().T
            store = Vector{Int64}() #collect all values contained in rkt
            for raw in res
                if k == raw.k && t == raw.t
                    val = getproperty(R[raw.r],q)[i,k,t]
                    push!(store,val)
                end
            end
            unique!(store)

            test = Vector{Int64}() #summarize the values, we average 1-by-1
            for w in 1:length(store)-1
                comp = ceil((store[w] + store[w+1]) / 2)
                push!(test,comp)
            end

            if !isempty(test)
                vmax = findmax(test)[1]
                vmin = findmin(test)[1]

                newmax = Vector{β}()
                push!(newmax,β(:k,k))
                push!(newmax,β(:t,t))
                push!(newmax,β(q,i,vmax))
                push!(maxBoundstack,newmax)

                newmin = Vector{β}()
                push!(newmin,β(:k,k))
                push!(newmin,β(:t,t))
                push!(newmin,β(q,i,vmin))
                push!(minBoundstack,newmin)
            end
        end
    end

    return maxBoundstack,minBoundstack
end

function Btest(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    newstack = separate(n)

<<<<<<< HEAD
    for stacks in newstack
        for k in keys(b().K), t in b().T
            for q in [:y,:z,:u,:v]
                totest = filter(p -> p[1].i == k && p[2].i == t && p[3].q == q,stacks)
                candidate = [last(p).i for p in totest]
                cardinality = 1

                while cardinality <= floor(log2(f([],R,θ))) + 1
                    combo = collect(combinations(candidate,cardinality))

                    for i in combo
                        tes = filter(p -> last(p).i in i,totest)
                        if !isempty(tes)
                            hehe = reduce(union,tes)
                            nilaii = s(hehe,R,θ)
                            if !Sniffles.issinteger(nilaii,1e-8)
                                return hehe
                            end
=======
    for q in [:y,:z,:u,:v]
        for k in keys(b().K), t in b().T
            totest = filter(p -> p[1].i == k && p[2].i == t && p[3].q == q,newstack)
            candidate = [last(p).i for p in totest]
            cardinality = 1

            while cardinality <= floor(log2(f([],R,θ))) + 1
                combo = collect(combinations(candidate,cardinality))

                for i in combo
                    tes = filter(p -> last(p).i in i,totest)
                    if !isempty(tes)
                        hehe = reduce(union,tes)
                        nilaii = s(hehe,R,θ)
                        if !Sniffles.issinteger(nilaii,1e-8)
                            return hehe
>>>>>>> e52615bd26b22b53ff45c202b381aad5e98800a4
                        end
                    end

                    #fractional not found increase cardinality
                    cardinality += 1
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
    seeds = Btest(n)
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
    #INITIALIZE DFS
    terminate = false
    stack = Vector{node}()
    visited = Vector{node}()
    iter = 0
    integer = node()

    #PUSH ROOT TO STACK
    push!(stack,n)

    #ITERATION
    while !(terminate || isempty(stack))
        if iter < maxiter
            #POP FROM STACK
            u = pop!(stack)

            #CHECK AND PROCESS
            if u.status[end] == "UNVISITED"
                #PROCESS THE NODE
                colGen(u;track=false,maxCG=Inf)

                #NODE STATUSES
                if u.status[end] == "EVALUATED" && integerCheck(u)
                    obj = objective_value(master(u))

                    #CHANGE UPPERBOUND
                    if obj < upperBound
                        upperBound = obj
                        push!(u.status,"NEW_BOUND")
                    end

                    #ALWAYS TERMINATE IF INTEGER NODE IS FOUND
                    terminate = true
                    println("we have integer solution.")
                    integer=u
                elseif u.status[end] == "EVALUATED"
                    obj = objective_value(master(u))

                    #BRANCH IF NOT ABOVE uB
                    if obj <= upperBound
                        append!(stack,createBranch(u))
                    end
                end

                #SEND TO VISITED
                push!(visited,u)
            end
            iter += 1
        else
            terminate = true
        end
    end

    return (stack=stack,visited=visited,upperBound=upperBound,integer=integer)
end
