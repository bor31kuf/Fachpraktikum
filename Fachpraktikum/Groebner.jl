using Oscar 
include("Version2.jl")

function Groebner(G, ordering::MonomialOrdering = default_ordering(parent(G[1])))
    LG =length(G)
    gMax = [leading_term(G[i],ordering=ordering) for i=1:LG]
    gMaxVec= [exponent_vector(gMax[i],1) for i=1:LG]
    for i=1:LG-1
        for j=i+1:LG-1
            println(j)
            Sij = SPoly(G[i],G[j],ordering)
            SijG = NewReduce(Sij,G,gMax,gMaxVec,ordering)
            if SijG != 0
                G = RReduce(SijG,G,ordering,gMax,gMaxVec)
            end
        end
    end
    return G
end

function RReduce(S,G,ordering,gMax,gMaxVec)
    push!(G,S)
    s = leading_term(S,ordering=ordering)
    push!(gMax,s)
    push!(gMaxVec,exponent_vector(s,1))
    for i=1:length(G)-1
        Sij = SPoly(G[i],S)
        SijG = NewReduce(Sij,G,gMax,gMaxVec,ordering)
        if SijG != 0
            G = RReduce(SijG,G,ordering,gMax,gMaxVec)
        end
    end
    return G
end



function SPoly(g,h, ordering::MonomialOrdering = default_ordering(parent(g)))
    ltg = leading_term(g,ordering =ordering)
    lth = leading_term(h,ordering=ordering)
    lmg = monomial(ltg,1)
    lmh = monomial(lth,1)
    mgh = ltg/gcd(lmg,lmh)
    mhg = lth/gcd(lmg,lmh)
    Sgh = mhg*g-mgh*h
    return Sgh
end

function my_isgb(L, ordering::MonomialOrdering  = default_ordering(parent(L[1])))
    t = length(L)
    for i = 1:t-1
        for j=i+1:t
            if DivRestFam(SPoly(L[i],L[j]),L,ordering) != 0
                return false
            end
        end
    end
    return true
end

function reduce_groebner(G,ordering::MonomialOrdering = default_ordering(parent(G[1])))
    x = true

    while x
        x =false
        i=1
        while i<=length(G) 
            G2 = copy(G)
            G2 = deleteat!(G2,i)
            A = [leading_term(G2[i],ordering=ordering) for i=1:length(G2)]
            A3 = [exponent_vector(A[i],1) for i=1:length(G2)]
            alpha = NewReduce(G[i],G2,A,A3,ordering)
            if alpha != G[i]
                if alpha == 0
                    deleteat!(G,i)
                    i-=1
                else
                    G[i] = alpha+G[i]*0
                x = true
                end
            end
            i+=1
        end
    end
    return G
end

function Groebner2(G,ordering::MonomialOrdering=default_ordering(parent(G[1])))
    H= copy(G)
    Groebner(H,ordering)
    """sort!(H,
          by=g->leading_monomial(g,ordering=ordering),
          order=ordering)"""
    H = reduce_groebner(H,ordering)
    sort!(H,
          by=g->leading_monomial(g,ordering=ordering),
          order=ordering)
    H = Oscar.IdealGens(H)
    return H
end