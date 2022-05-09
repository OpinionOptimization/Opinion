struct Graph
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int, 1}
    v :: Array{Int, 1} # uv is an edge
#    w :: Array{Float64, 1} # weight of each edge
    nbr :: Array{Array{Int, 1}, 1}
end



function get_graph_direct(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()
    #E = Set{Tuple{Int32, Int32, Float32}}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1

    fname = string("data/",ffname)
    fin = open(fname, "r")


    str = readline(fin)
    str = split(str)
    #n   = parse(Int, str[1])
    m   = parse(Int, str[3])
    u = Int[]
    v = Int[]

    tot = 0
    for i = 1 : m
        str = readline(fin)
        str = split(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        if x!=y
            u1 = getID(x)
            v1 = getID(y)
            Origin[u1] = x
            Origin[v1] = y
            push!(u, u1)
            push!(v, v1)
            tot += 1
        end
    end
    nbr=[ [ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr[u1],v1);
        #push!(nbr[v1],u1);
    end

    close(fin)
    return Graph(n, tot, u, v,nbr)
end

function lap_direct(G :: Graph)
    F = zeros(G.n, G.n);
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        #F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        #F[G.v[i], G.v[i]] += 1
    end
    return F
end
