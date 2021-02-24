# =========================================================================
#    BASIC DATA GENERATION
# =========================================================================

const template_data = Ref{Any}(nothing)
b() = template_data[] #buat manggil dt default

function extract!(path::String;distformula::String) #extract from excel
    xf = XLSX.readxlsx(path) #READ WORKSHEET
    data = Dict{Symbol,DataFrame}() #DATAFRAME DICT

    for sheets in XLSX.sheetnames(xf) #TURN SHEETS INTO DATAFRAME
        df = DataFrame(XLSX.gettable(xf[sheets])...) #TRANSFORM TABLE INTO DATAFRAME
        data[Symbol(sheets)] = df #DEFINE THE NAME FROM THE WORKSHEET
    end

    V = Dict{Int64,vtx}() #INITIATE VERTICES
    for v in eachrow(data[:vertices]) #ITERATE OVER DATA
        V[v.id] = vtx(
            v.name, v.type,
            v.x, v.y, v.MAX, v.MIN, v.START,
            v.h
        )
    end

    dist = JuMP.Containers.DenseAxisArray{Float64}(undef, keys(V), keys(V))
    for i in keys(V), j in keys(V)
        if distformula == "haversine"
            if i != j
                dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],6378.137)
            else
                dist[i,j] = 999999999
            end
        elseif distformula == "euclidean"
            if i != j
                dist[i,j] = round(euclidean([V[i].x,V[i].y],[V[j].x,V[j].y]))
            else
                dist[i,j] = 999999999
            end
        end
    end

    K = Dict{Int64,veh}() #INITIATE VEHICLES
    for k in eachrow(data[:vehicles]) #ITERATE OVER DATA
        K[k.id] = veh(
            k.name, k.type,
            parse.(Int64,split(k.cover)),
            Dict(parse.(Int64,split(k.cover)) .=>  parse.(Int64,split(k.BP))), k.Q,
            k.vx, k.vl, k.fp, k.fd
        )
    end

    T = collect( #range from starting month for duration
        range(
            last(data[:periods].start),
            length = last(data[:periods].T),
            step = 1
        )
    )

    d = JuMP.Containers.DenseAxisArray(
        Array{Float64}(data[:demands][:,string.(T)]), #dataset
        Array{Int64}(data[:demands].point), #dims 1
        T #dims 2
    )

    res = dt(V,dist,K,T,d) #WRAPPING RESULT TO TYPE dt

    return template_data[] = res
end

function initStab(res = b())
    slackCoeff = sl_C()
    surpCoeff = su_C()
    slackLim = abs.(res.d)
    surpLim = abs.(res.d)

    return stabilizer(slackCoeff,surpCoeff,slackLim,surpLim)
end

function root(res = b())
    id = uuid1()
    root = node(
        id, id,
        Vector{bound}(),Vector{col}(),
        initStab(),
        ["UNVISITED"]
    )

    return root
end
