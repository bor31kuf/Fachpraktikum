using Oscar,DataStructures,SIMD,Revise,BenchmarkTools



struct ReamainderOrdering <: Base.Order.Ordering
end

function Ordnung(a,b)
    for i in 1:4
        if a[1][i] < b[1][i]
            return false
        elseif a[1][i] > b[1][i]
            return true
        end
    end
    return false
end
import Base.Order.lt
lt(o::ReamainderOrdering,a,b) = Ordnung(a,b)



struct PolyNom
    PolyNOM::BinaryHeap{}
end
 

PolAlg, (x,y,z) = polynomial_ring(QQ,[:x,:y,:z],internal_ordering=:lex)


function PolNeu(f,c::Vec{4,Int64})
    #besser Speichern
    o= lex([x,y,z])
    A = collect(coefficients(f,ordering=o))
    B = collect(exponents(f,ordering=o))
    X = []
    for i=1:length(A)
        push!(X,[Vec{4,Int64}((sum([B[i][j]*c[j+1] for j=1:length(c)-1]),B[i][1],B[i][2],B[i][3])),A[i]])
    end
    #println(X)
    F = PolyNom(BinaryHeap(ReamainderOrdering(),X))
    return F 
end  

#println(PolNeu(x*y+z,Vec{4,Int64}((0,1,1,1))))


function Vgl(a::Vec{4,Int64},b::Vec{4,Int64})
    for i in 1:4
        if a[1][i] < b[1][i]
            return 0
        elseif a[1][i] > b[1][i]
            return 1
        end
    end
    return 2
end

function NeuPol(f::PolyNom) 
    a = x*0
    for i=1:length(f.PolyNOM)
        m = pop!(f.PolyNOM)
        a += monomial(PolAlg,collect(Tuple([m[1][j] for j= 2:4])))*m[2]
    end
    return a
end


function DIV2(f::PolyNom,G::Array{PolyNom})
    #println((f.PolyNOM).valtree)
    #f2 = PolyNom(deepcopy(f.PolyNOM))
    f2 = PolyNom(BinaryHeap(ReamainderOrdering(),copy(f.PolyNOM.valtree)))
    #collect(f2.PolyNOM)
    r = PolyNom(BinaryHeap(ReamainderOrdering(),[[Vec{4,Int64}((0,0,0,0)),QQFieldElem(0)]]))
    L = length(G)
    o=1
    A= [(G[i].PolyNOM).valtree for i=1:length(G)]
    while !isempty(f2.PolyNOM)
        w = false
        for i=o:L
            if sum(first(f2.PolyNOM)[1]>=first(G[i].PolyNOM)[1])==4
                DIV1 = first(f2.PolyNOM)[1]-A[i][1][1] 
                DIV2 = -first(f2.PolyNOM)[2]/A[i][1][2]
                L2 = length(A[i])
                j=2
                pop!(f2.PolyNOM)
                for j =2:L2
                    push!(f2.PolyNOM,[A[i][j][1]+DIV1,A[i][j][2]*DIV2])
                end
                w =true
                break
            end 
        end
        if w == false
            push!(r.PolyNOM,first(f2.PolyNOM))
            pop!(f2.PolyNOM)
        end
    end
    return r
end


function Mul1(f::PolyNom,g::Vec{4,Int64},l)
    h =deepcopy(f)
    for i=1:length(f.Monome)
        h.Monome[i] += g
        h.Koeffizienten[i] *= l
    end
    return h
end




f1 =x*y
f2 =x^3+z
g = x^12+x*z^3*y+x^4*z+z^13+z^7*y
c = Vec{4,Int64}((0,0,0,0))
A = PolNeu(g,c)
B = PolyNom[PolNeu(f2,c),PolNeu(f1,c)]
function MinTest1()
    return NeuPol(DIV2(A,B))
    
end
function MinTest2()
    return divrem(g,[f1,f2])
end

println(MinTest1())
println(MinTest2())

"""A = [[rand(PolAlg,-1:2,0:4,-10:6) for j=1:12] for i=1:30]
B = []
while length(B) != 30
    m= rand(PolAlg,-1:8,3:7,1:10)
    if m!= []
        push!(B,m)
    end
end
Ane = [[PolNeu(A[i][j],Vec{4,Int64}((0,0,0,0))) for j=1:12] for i = 1:30]
Bne = [PolNeu(B[i],Vec{4,Int64}((0,0,0,0))) for i= 1:30]
#println(Bne)
for i=1:30
    j=1
    while j<= length(Ane[i])
        #println(Ane)
        if isempty(Ane[i][j].PolyNOM)
            deleteat!(Ane[i],j)
            deleteat!(A[i],j)
        else
            j+=1
        end
    end
end


function Test1()
    for i=1:30
        #println(Bne[i])
        #println(Ane[i])
        DIV2(Bne[i],Ane[i])
        #println(divrem(B[i],A[i])[2])
    end
end

function Test2()
    for i=1:30
        divrem(B[i],A[i])
    end
end

Test1()

println(@benchmark Test1())
println(@benchmark Test2())"""