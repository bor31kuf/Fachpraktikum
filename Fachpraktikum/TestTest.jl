using Oscar
using BenchmarkTools
include("Groebner.jl")

PolAlg, (x1,x2,x3,x4) =  polynomial_ring(QQ,[:x1;:x2;:x3;:x4])
G = [x1 + 2*x2 + 2*x3 - 1,x1^2 - x1 + 2*x2^2 + 2*x3^2,2*x1*x2 + 2*x2*x3 - x2]

#B = gens(Groebner2(G,lex(PolAlg)))
#println(B)

#A= gens(groebner_basis(ideal(G),ordering=lex(PolAlg),complete_reduction=true))
#println(A)

function Test()
    G = []
    G = [x1 + 2*x2 + 2*x3 - 1,x1^2 - x1 + 2*x2^2 + 2*x3^2,2*x1*x2 + 2*x2*x3 - x2]
    B = gens(Groebner2(G,lex(PolAlg))) 
end

function Test2()
    G = []
    G = [x1 + 2*x2 + 2*x3 - 1,x1^2 - x1 + 2*x2^2 + 2*x3^2,2*x1*x2 + 2*x2*x3 - x2]
    B = gens(groebner_basis(ideal(G),ordering=lex(PolAlg),complete_reduction=true))
end