using SIMD,Oscar,Revise,BenchmarkTools,DataStructures

include("CircularSIMDNormal.jl")



using .SIMD4

n=3


vars = ["x$i" for i in 1:n]
PolAlg, gens = polynomial_ring(QQ,vars,internal_ordering=:degrevlex)
ord = degrevlex(PolAlg)

println(typeof(ord))
A = [[(rand(PolAlg,-1:2,0:4,-10:6)) for j=1:8] for i=1:30]
B = []
while length(B) != 30
    m= rand(PolAlg,-1:8,3:11,1:30)^3
    if m!= []
        push!(B,m)
    end
end




c = Vec{n+1,Int64}(ntuple(i->1,n+1))
Ane = [[SIMD4.PolNeu(A[i][j],ord=ord) for j=1:8] for i = 1:30]
Bne = [SIMD4.PolNeu(B[i],ord=ord) for i= 1:30]
Anee = [[SIMD9.PolNeu(A[i][j],ord=ord) for j=1:8] for i = 1:30]
Bnee = [SIMD9.PolNeu(B[i],ord=ord) for i= 1:30]


for i=1:30
    j = 1
    while j<= length(Ane[i])
        while length(Ane[i][j].Monome)==0
            A[i][j] = rand(PolAlg,-1:2,0:4,-10:6)
            Ane[i][j] = SIMD4.PolNeu(A[i][j],ord=ord)
            Anee[i][j] = SIMD9.PolNeu(A[i][j],ord=ord)
        end
        j+=1
    end
end


function Test1()
    for i=1:length(Ane)
        SIMD4.DIV2(Bne[i],Ane[i])
    end
    return "fanta"
end 


function Test2()
    for i=1:length(Ane)
        divrem(B[i],A[i])
    end
end


println(@btime Test1())
println(@btime Test2())
