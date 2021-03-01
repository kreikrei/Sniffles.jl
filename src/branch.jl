# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================
function separate(R,θ,idx,Seq,record)
    F = fract(Seq,R,θ)

    if isempty(F)
        return record
    else
        xagg = Dict()
        for edge in idx
            xagg[edge] = sum(
                R[f.r][(f.k,f.t)].x[first(edge),last(edge)] * θ[f.r,f.k,f.t] for f in F
            )
        end

        found = false
        for edge in idx
            if abs(xagg[edge] - round(xagg[edge])) > 1e-8
                #println("$edge : $(xagg[edge])")
                Seq0 = S(Seq.k,Seq.t,vcat(Seq.sequence,β(:x,edge,1,1)))
                push!(record, Seq0)
                #println(record)
                found = true
            end
            if found
                #println(record)
                return record
            end
        end

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
    end
end

function assess(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    res = Vector{S}()
    for k in K(), t in T()
        Seq = S(k,t,β[])
        idx = edges(k)
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
