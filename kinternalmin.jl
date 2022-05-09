include("graph.jl")
using LinearAlgebra
using Random
using Statistics
using Laplacians
using SparseArrays
using DataFrames


function influence_centrality(G, tau)
    d = [length(G.nbr[j]) for j = 1:G.n];
    InForest = [false for j = 1:G.n];
    Next = [-1 for j = 1:G.n];
    rho = zeros(G.n);
    rootindex = zeros(Int, G.n);
    for p = 1:tau
        InForest .= false;
        Next .= -1;
        rootindex .= 0;
        for i = 1:G.n
            u = i;
            while !InForest[u]
                seed = rand();
                if seed <= 1 / (1 + d[u])
                    InForest[u] = true;
                    Next[u] = -1;
                    rootindex[u] = u;
                    rho[u] += 1;
                else
                    k = floor(Int, seed * (1 + d[u]));
                    Next[u] = G.nbr[u][k];
                    u = Next[u];
                end
            end
            rootnow = rootindex[u];
            u = i;
            while !InForest[u]
                InForest[u] = true;
                rootindex[u] = rootnow;
                rho[rootnow] += 1;
                u = Next[u];
            end
        end
    end
    return rho ./ tau;
end

function influence_centrality_real(G)
    rho_real = zeros(G.n);
    L = lap_direct(G);
    for i = 1:G.n
        L[i, i] += 1;
    end
    L = inv(L);
    for i = 1:G.n
        for j = 1:G.n
            rho_real[i] += L[j, i];
        end
    end
    return rho_real
end

function k_internal(G,s,tau,k)
    ansnode = [];
    ansvalue = 0;
    rho = influence_centrality(G,tau);
    rho = rho.*s;
    for i = 1:k
        p = findmax(rho)[2][1];
        rho[p] = 0;
        push!(ansnode,p);
    end
    ansvalue = sum(rho)/G.n;
    return ansnode
end

function k_internal_real(G,s,k,rho)
    ansnode = [];
    ansvalue = 0;
    rho = rho.*s;
    for i = 1:k
        p = findmax(rho)[2][1];
        rho[p] = 0;
        push!(ansnode,p);
    end
    ansvalue = sum(rho)/G.n;
    return ansnode
end

function k_internal_random(G,s,k,rho)
    ansnode = rand(1:G.n,k);
    ansvalue = 0;
    rho = rho.*s;
    for i = 1:k
        p = ansnode[i];
        rho[p] = 0;
    end
    ansvalue = sum(rho)/G.n;
    return ansvalue
end



function runtime(filename) #effectiveness
    l = length(filename);
    n = zeros(l);
    m = zeros(l);
    rf500value = zeros(l);
    rf500time = zeros(l);
    rf1000value = zeros(l);
    rf1000time = zeros(l);
    rf2000value = zeros(l);
    rf2000time = zeros(l);
    rf500err = zeros(l);
    rf1000err = zeros(l);
    rf2000err = zeros(l);
    exacttime = zeros(l);
    exactvalue = zeros(l);
    for i=1:l
        G = get_graph_direct(filename[i]);
        n[i] = G.n;
        m[i] = G.m;
        s = rand(G.n,1);
        k=50;
        tau = 500;
        t1=time();
        ansrf = k_internal(G,s,tau,k);
        t2=time();
        rf500time[i] = t2-t1;
        rf500value[i] = ansrf;
        tau = 1000;
        t1=time();
        ansrf = k_internal(G,s,tau,k);
        t2=time();
        rf1000time[i] = t2-t1;
        rf1000value[i] = ansrf;
        tau = 2000;
        t1=time();
        ansrf = k_internal(G,s,tau,k);
        t2=time();
        rf2000time[i] = t2-t1;
        rf2000value[i] = ansrf;
        t2=time();
        rho = influence_centrality_real(G);
        ansreal = k_internal_real(G,s,k,rho);
        t3=time();
        exacttime[i] = t3-t2;
        exactvalue[i] = ansreal;
        rf500err[i] = abs(rf500value[i]-exactvalue[i])/exactvalue[i];
        rf1000err[i] = abs(rf1000value[i]-exactvalue[i])/exactvalue[i];
        rf2000err[i] = abs(rf2000value[i]-exactvalue[i])/exactvalue[i];
    end
    df = DataFrame(name = filename, n=n,m=m,rf500value = rf500value,
    rf500time = rf500time,
    rf1000value = rf1000value,
    rf1000time = rf1000time,
    rf2000value = rf2000value,
    rf2000time = rf2000time,
    exacttime = exacttime,
    exactvalue = exactvalue,
    rf500err = rf500err,
    rf1000err = rf1000err,
    rf2000err = rf2000err)
    return df
end



function runtime2(filename)
    l = length(filename);
    n = zeros(l);
    m = zeros(l);
    rf1000value = zeros(l);
    rf1000time = zeros(l);
    rf500value = zeros(l);
    rf500time = zeros(l);
    rf2000value = zeros(l);
    rf2000time = zeros(l);
    for i=1:l
        G = get_graph_direct(filename[i]);
        n[i] = G.n;
        m[i] = G.m;
        s = rand(G.n,1);
        k=50;
        tau = 500;
        t1=time();
        ansrf = k_internal(G,s,tau,k);
        t2=time();
        rf500time[i] = t2-t1;
        rf500value[i] = ansrf;
        tau = 1000;
        t1=time();
        ansrf = k_internal(G,s,tau,k);
        t2=time();
        rf1000time[i] = t2-t1;
        rf1000value[i] = ansrf;
        tau = 2000;
        t1=time();
        ansrf = k_internal(G,s,tau,k);
        t2=time();
        rf2000time[i] = t2-t1;
        rf2000value[i] = ansrf;
    end
    df = DataFrame(name=filename,n=n,m=m,rf500value = rf500value,
    rf500time = rf500time,
    rf1000value = rf1000value,rf1000time = rf1000time,rf2000value = rf2000value,rf2000time = rf2000time)
    return df
end




function runeffect(G,s,tau,rho) #k = 10,20,30,40,50
    #exact, rf, s-order, d-order，random,z-order
    k=50;
    effect = zeros(5,6);
    U = k_internal_real(G,s,k,rho);
    ss = deepcopy(s);
    for i = 1:5
        for j = 1:10
            ss[U[10*(i-1)+j]] =0;
        end
        effect[i,1] = sum(rho.*ss)/G.n;
    end

    U = k_internal(G,s,tau,k);
    ss = deepcopy(s);
    for i = 1:5
        for j = 1:10
            ss[U[10*(i-1)+j]] =0;
        end
        effect[i,2] = sum(rho.*ss)/G.n;
    end

    ss = deepcopy(s);
    U = [];
    for i = 1:k
        index = findmax(ss)[2][1];
        ss[index] = 0;
        push!(U,index);
    end
    ss = deepcopy(s);
    for i = 1:5
        for j = 1:10
            ss[U[10*(i-1)+j]] =0;
        end
        effect[i,3] = sum(rho.*ss)/G.n;
    end


    din = zeros(G.n);
    for i = 1:G.m
        din[G.v[i]] += 1;
    end
    #d = [length(G.nbr[j]) for j = 1:G.n];
    U = [];
    for i = 1:k
        index = findmax(din)[2];
        push!(U,index);
        din[index] = 0;
    end
    ss = deepcopy(s);
    for i = 1:5
        for j = 1:10
            ss[U[10*(i-1)+j]] =0;
        end
        effect[i,4] = sum(rho.*ss)/G.n;
    end

    U= [];
    while length(U)<50
        a = rand(1:G.n,1);
        if !(a in U)
            push!(U,a[1])
        end
    end
    ss = deepcopy(s);
    for i = 1:5
        for j = 1:10
            ss[U[10*(i-1)+j]] =0;
        end
        effect[i,5] = sum(rho.*ss)/G.n;
    end


    L = lap_direct(G);
    for i = 1:G.n
        L[i, i] += 1;
    end
    L = inv(L);
    z = L*s;
    U = [];
    for i = 1:k
        index = findmax(z)[2][1];
        push!(U,index);
        z[index] = 0;
    end
    ss = deepcopy(s);
    for i = 1:5
        for j = 1:10
            ss[U[10*(i-1)+j]] =0;
        end
        effect[i,6] = sum(rho.*ss)/G.n;
    end
    return effect
end
function runeffectuniform(filename)
    #filename = "Humanproteins.txt";
    G = get_graph_direct(filename);
    rho = influence_centrality_real(G);
    s = rand(G.n,1);
    tau = 500;
    eff = runeffect(G,s,tau,rho);
end



function runeffectgauss(filename)
    #filename = "P2p-Gnutella08.txt";
    G = get_graph_direct(filename);
    rho = influence_centrality_real(G);
    s = rand(G.n,1);
    for i = 1:G.n
        s[i] = randn(1)[1];
    end
    a1 = findmax(s)[1];
    a2 = findmin(s)[1];
    k = 1/(a1-a2);
    b = 1-a1/(a1-a2);
    for i = 1:G.n
        s[i] = k*s[i]+b;
    end
    tau = 500;
    eff = runeffect(G,s,tau,rho);
end
