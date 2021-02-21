# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function separate(n::node)
    #PROTOTYPE SEPARATION
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    F = Vector{NamedTuple}()
    for r in keys(R), k in keys(b().K), t in b().T
        if θ[r,k,t] - floor(θ[r,k,t]) > 0
            push!(F,(r=r,k=k,t=t))
        end
    end #semua q yang nilai θ-nya masih fractional

    qF = DataFrame(q = Symbol[], i = Int64[], val = Int64[])
    for f in F, q in [:u,:v,:y,:z], i in keys(b().V)
        append!(qF,
            DataFrame(
                q = q,
                i = i,
                val = getproperty(R[f.r],q)[i,f.k,f.t]
            )
        )
    end #collect existing values of the variable (q,i)

    return qF
end

function vtest(q::Symbol,i::Int64,qF::DataFrame)
    test_v = Vector{Int64}()

    distinct = filter(p -> p.i == i && p.q == q,qF).val
    push!(distinct,0)
    sort!(distinct)

    for w in 1:length(distinct)-1
        comp = ceil((distinct[w] + distinct[w+1]) / 2)
        if comp > 0
            push!(test_v, comp)
        end
    end

    res = Vector{β}()

    for v in unique!(test_v)
        push!(res,β(q,i,v))
    end

    return reverse(res)
end

function qstack()
    res = Vector{β}()

    for q in [:u,:v,:y,:z], i in keys(b().V)
        push!(res,β(q,i,nothing))
    end

    return res
end

issinteger(val,tol) = abs(round(val) - val) < tol

function Btest(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    test_stack = qstack()
    test_B = Vector{β}()
    lastf1 = deepcopy(test_B)
    terminate = false
    sep = separate(n)

    while !terminate
        if isempty(test_stack)
            append!(test_B,lastf1)
            test_stack = qstack()
        else
            test_β = pop!(test_stack)
            #println(test_β)
            #println("stack sisa: $(length(test_stack))")

            if isnothing(test_β.v)
                append!(test_stack,vtest(test_β.q,test_β.i,sep))
                #println("new stack added")
            else
                push!(test_B,test_β)

                fract = f(test_B,R,θ)
                if fract >= 1
                    lastf1 = deepcopy(test_B)
                end

                if !issinteger(s(test_B,R,θ),1e-8)
                    break
                else
                    pop!(test_B)
                end
            end
        end
    end

    return test_B
end

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

    for br in ["≳","≲"]
        push!(branches,
            node(
                n.self, #parent
                uuid1(), #self

                vcat(n.bounds, #bounds
                    bound(
                        seeds,br,
                        if br == "≲"
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

    return (stack=stack,visited=visited,upperBound=upperBound)
end
