# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================
function separate(R,θ,idx,Seq,record)
    F = fract(Seq,R,θ)

    if isempty(F)
        return record
    else
        xagg = Dict()
        for i in idx
            xagg[i] = sum(
                R[f.r][(f.k,f.t)].u[i] * θ[f.r,f.k,f.t] for f in F
            )
        end

        found = false
        for i in idx
            if abs(xagg[i] - round(xagg[i])) > 1e-8
                println("$edge : $(xagg[edge])")
                Seq0 = S(Seq.k,Seq.t,vcat(Seq.sequence,β(:u,i,1,floor(xagg[i]))))
                push!(record, Seq0)
                #println(record)
                found = true
            end
            if found
                #println(record)
                return record
            end
        end
        #=
        J = filter(p -> 0 < xagg[p] < sF(Seq,R,θ),idx)
        if !isempty(J)
            star = pop!(J)

            #CREATE NEW INPUT
            Seq1 = S(Seq.k,Seq.t,vcat(Seq.sequence,β(:x,star,1,1)))
            Seq2 = S(Seq.k,Seq.t,vcat(Seq.sequence,β(:x,star,-1,1)))
            record = separate(R,θ,J,Seq1,record)
            record = separate(R,θ,J,Seq2,record)
        end

        println(record)
        return record
        =#
    end
end

function assess(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    res = Vector{S}()
    for k in K(), t in T()
        Seq = S(k,t,β[])
        idx = K(k).cover
        record = Vector{S}()

        final = separate(R,θ,idx,Seq,record)
        println(final)
        if !isempty(final)
            res = final
            return res
        end
    end

    return res
end


issinteger(val,tol) = abs(round(val) - val) < tol

function integerCheck(n::node)
    integer = true
    ori = origin(n)

    for k in K(), t in T(), i in K(k).cover, j in K(k).cover
        val = ori.x[i,j,k,t]
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

    candidate = assess(n)
    seeds = candidate[1]

    println(seeds)

    val = sQ(seeds,R,θ)
    println("branch on $seeds: $val")

    for br in ["<=",">="]
        push!(branches,
            node(
                n.self, #parent
                uuid1(), #self

                vcat(n.bounds, #bounds
                    bound(
                        seeds,br,
                        if br == "<="
                            floor(sQ(seeds,R,θ))
                        else
                            ceil(sQ(seeds,R,θ))
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
