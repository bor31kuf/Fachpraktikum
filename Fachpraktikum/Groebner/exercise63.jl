## Blatt 6, Aufgabe 3
## Joris Wirsen, Helena Belzl

"""using Oscar
include("/home/joris/Dokumente/Fachpraktikum/Version2.jl")

include("/home/joris/Dokumente/Fachpraktikum/Groebner/exercise62.jl")
include("exercise61.jl")

function my_gb(L)
    i = length(L)
    i2 = 2
    L2 = deepcopy(L)
    while true
        for k=1:i
            for j=min(k+1,i2):i
                if gcd(leading_monomial(L2[j]),leading_monomial(L2[k])) != 1
                    S = my_spoly(L2[j],L2[k])
                    S2 = (L2,S)
                    if S2 != 0
                        L2 = push!(L2,S2)
                    end
                end
            end
        end
        if i == length(L2)
            break
        end
        i2 = i
        i = length(L2)
    end
    return L2
end

R, (x,y,z) = QQ[:x,:y,:z]
G = [x^2*y^2-x^2*y+x,x^2*y+y^2]
println(my_isgb(my_gb(G)))"""

#test

"function test_ord(G1)
    for i=1:length(G1)
        my_gb(G1[i])
    end
    return
end
function my_isgb(L)
    t = length(L)
    for i = 1:t-1
        for j=i+1:t
            if div_mit_rest(L, my_spoly(L[i],L[j])) != 0
                return false
            end
        end
    end
    return true
end"


"function test()
    G = []
    for i = 1:4
        G1 = []
        for k = 1:1
            G1 = push!(G1,gens(katsura(i)))
        end
        G = push!(G,G1)
    end
    R1, (x1,x2) = QQ[:x1,:x2]
    R2, (x1,x2,x3) = QQ[:x1,:x2,:x3]
    R3, (x1,x2,x3,x4) = QQ[:x1,:x2,:x3,:x4]
    R4, (x1,x2,x3,x4,x5) = QQ[:x1,:x2,:x3,:x4,:x5]
    R = [R1,R2,R3,R4]
    for k=1:4
        F = G[k]
        with_ordering(R[k],lex(R[k])) do
            @time test_ord(F)
        end
    end
    for k=1:4
        F = G[k]
        with_ordering(R[k],deginvlex(R[k])) do
            @time test_ord(F)
        end
    end
    K = []
    for i = 1:4
        K1 = []
        for k = 1:10
            K1 = push!(K1,gens(katsura(i)))
        end
        K = push!(K,K1)
    end
    for k=1:4
        for i=1:length(K[k])
            if my_isgb(my_gb(K[k][i])) != true
                return false
            end
        end
    end
    return
end"




