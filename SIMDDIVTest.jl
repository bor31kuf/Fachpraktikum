using SIMD,Oscar,Revise,BenchmarkTools,DataStructures

include("SIMDCircular.jl")
include("NotCIrc.jl")



using .SIMD4
using .SIMD9

n=5


vars = ["x$i" for i in 1:n]
PolAlg, gens = polynomial_ring(QQ,vars,internal_ordering=:degrevlex)
ord = degrevlex(PolAlg)
#PolAlg2, gens2 = polynomial_ring(QQ,vars,internal_ordering=ord)
println(typeof(ord))
A = [[(rand(PolAlg,-1:2,0:4,-10:6)) for j=1:8] for i=1:30]
B = []
while length(B) != 30
    m= rand(PolAlg,-1:8,3:11,1:30)^3
    if m!= []
        push!(B,m)
    end
end

#groebner_basis(ideal(A[1]),ordering=ord)


c = Vec{n+1,Int64}(ntuple(i->1,n+1))
Ane = [[SIMD4.PolNeu(A[i][j],ord=ord) for j=1:8] for i = 1:30]
Bne = [SIMD4.PolNeu(B[i],ord=ord) for i= 1:30]
Anee = [[SIMD9.PolNeu(A[i][j],ord=ord) for j=1:8] for i = 1:30]
Bnee = [SIMD9.PolNeu(B[i],ord=ord) for i= 1:30]

#println(Bne)  
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
    #Bne = [SIMD4.PolNeu(B[i],c,ord=ord) for i= 1:30]
    for i=1:length(Ane)
        SIMD4.DIV2(Bne[i],Ane[i])
        #println("hi")
        #println(length(Bne[i].Koeffizienten))

        #println(collect(exponents(SIMD4.NeuPol(SIMD4.DIV2(Bne[i],Ane[i]),PolAlg))) == collect(exponents(divrem(B[i],A[i])[2])))
        
        #if SIMD4.NeuPol(SIMD4.DIV2(Bne[i],Ane[i]),PolAlg,ord=ord) != divrem(B[i],A[i])[2]
         #  println("warum")
          #  break
        #end

    end
    return "fanta"
end 

function Test3()
    for i =1:length(Ane)
        SIMD9.DIV2(Bnee[i],Anee[i])
        #if SIMD9.NeuPol(SIMD9.DIV2(Bnee[i],Anee[i]),PolAlg) != divrem(B[i],A[i])[2]   
         #   println("warum")
          #  break
        #end
    end
    return "fanta"
end

function Test2()
    for i=1:length(Ane)
        divrem(B[i],A[i])
    end
end


#println(@benchmark(Test3()))
#println(B[1])

#println(B[1])

#Test3()
#println(Bne[1])
#Test1()
println("Was")

#println(Bne[1]) 

#println("hi")
#println(Test1())
println(@btime Test1())
println(@btime Test3())
#println(Bne[1])
println(@btime Test2())
#if SIMD4.NeuPol(SIMD4.DIV2(Bne[i],Ane[i]),PolAlg) != divrem(B[i],A[i])[2]
#println("warum")
#break
#end
