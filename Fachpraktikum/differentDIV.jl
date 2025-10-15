module SIMD6


export PolNeu, DIV2, NeuPol, Sub1




using Oscar,SIMD,BenchmarkTools,DataStructures


mutable struct PolyNom{W}
    Monome::CircularDeque{Vec{W,Int64}}
    Koeffizienten::CircularDeque{FieldElem}
end

PolAlg, (x1,x2,x3,x4,x5,x6,x7) = polynomial_ring(QQ,[:x1,:x2,:x3,:x4,:x5,:x6,:x7],internal_ordering=:lex)
ord=lex(PolAlg)

function PolNeu(f,c::Vec{W,Int64};ord::MonomialOrdering=default_ordering(parent(f))) where {W}
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    D = PolyNom(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
    for i=1:length(A)
        push!(D.Monome,Vec{W,Int64}((sum([B[i][j]*c[j+1] for j=1:length(c)-1]),B[i]...)))
        push!(D.Koeffizienten,A[i])
    end
    return D

end  

function Vgl(a::Vec{W,Int64},b::Vec{W,Int64}) where{W}
    eq = a== b
    gt = a>b
    
    for i in 1:W
        if !eq[i]
            return Int(gt[i])
        end
    end
    return 2
end

function NeuPol(f::PolyNom,PolAlg) 
    a=zero(PolAlg)
    k = length(f.Monome)
    for i=1:k
        a += monomial(PolAlg,collect(Tuple(f.Monome[i]))[2:end])*f.Koeffizienten[i]
    end
    return a
end

struct ReamainderOrdering <: Base.Order.Ordering
end
import Base.Order.lt
lt(o::ReamainderOrdering,a,b) = lex_less(a,b)

function DIV2(f::PolyNom{W},G::Vector{PolyNom{W}}) where W
    #f2 = PolyNom(deepcopy(f.Monome),deepcopy(f.Koeffizienten))
    L = length(f.Koeffizienten)
    if L==0
        return f
    end
    f2 = geobucketpol([PolyNom(CircularDeque{Vec{W,Int64}}(1),CircularDeque{FieldElem}(1))])
    A = CircularDeque(copy(f.Monome.buffer),f.Monome.capacity,f.Monome.n,f.Monome.first,f.Monome.last)
    B = CircularDeque(copy(f.Koeffizienten.buffer),f.Koeffizienten.capacity,f.Koeffizienten.n,f.Koeffizienten.first,f.Koeffizienten.last)

    f2 =addgeobucket(f2,PolyNom(A,B))

    LTf2 = Leitterm(f2)
    a = 1
    b=1
    r = PolyNom(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
    while true
        w = false
        for i=1:D
              
            if sum((first(LTf2.Monome)>=first(G[i].Monome)))==W
                
                
                DIV2 = first(LTf2.Koeffizienten)/first(G[i].Koeffizienten)
                

                #Binomische Formel implementieren
                n = sum(first(LTf2.Monome))
                for t=2:W
                    if first(G[i].Monome)[t] != 0 && first(LTf2.Monome)[t] != 0
                        n2 = first(LTf2.Monome)[t]/first(G[i].Monome)[t]
                        if n2 < n
                            n = Int(floor(n2))
                        end
                    end
                end
                m = length(G[i].Monome)
                Z =collect(weak_compositions(n,m-1))
                nc = length(Z)
                A= CircularDeque{Vec{W,Int64}}(nc)
                B = CircularDeque{FieldElem}(nc)
                DIV1 = first(LTf2.Monome)-n*first(G[i].Monome)
                for t in 1:nc
                    Koeff = factorial(big(n))
                    #println(n)
                
                    for j = 1:m-1
                        Koeff //= factorial(big(Z[t][j]))
                    end
                    Mono = Vec{W,Int64}(ntuple(i->0,W))
                    coeff = DIV2*Koeff
                    for j=1:m-1
                        Mono +=G[i].Monome[j+1]*Z[t][j]
                        coeff *= G[i].Koeffizienten[j+1]^Z[t][j]
                    end
                    coeff *=(-1)^n
                    Mono += DIV1
                    push!(A,Mono)
                    push!(B,coeff)
                end

                A2 = A.buffer
                B2 = B.buffer
                M = sortperm(eachindex(A2);lt = (x,y) -> begin
                    n = length(x)
                    @inbounds for i in 1:n
                        if x[i] != y[i]
                            return x[i] < y[i]
                        end
                    end
                    return false 
                end)
                A2 =A2[M]
                B2 = B2[M]


                w = true
                if length(A)!=0
                    g = PolyNom(A,B)
                    #println(DIV1)
                    #println(n)
                    #println(NeuPol(g,PolAlg))
                    
                    f2= addgeobucket(f2,g)
                end
                a+=1
                b+=1

                LTf2 = Leitterm(f2)
                if isempty(LTf2.Koeffizienten)
                    #println(a," ",b)
                    return r
                end
                break
            
                emtpy!(A)
                empty!(B)
            end
        end
        if w == false
            a+=1
            r = pushing(r,LTf2)
            LTf2 = Leitterm(f2)
            if isempty(LTf2.Koeffizienten)
                #println(a, " ",b)
                return r
            end
        end
    end
end
 
function lex_less(a,b)
    n = length(a[1])
    @inbounds for i in 1:n
        if a[1][i] != b[1][i]
            return a[1][i] < b[1][i]
        end
    end
    return false
end

function pushing(r::PolyNom{W},f::PolyNom{W}) where{W}
    if capacity(r.Koeffizienten) > length(r.Koeffizienten)
        push!(r.Monome,first(f.Monome))
        push!(r.Koeffizienten,first(f.Koeffizienten))
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.Monome))
        r22 = CircularDeque{FieldElem}(2*capacity(r.Monome))
        for i=1:length(r.Koeffizienten)
            push!(r21,r.Monome[i])
            push!(r22,r.Koeffizienten[i]) 
        end
        push!(r21,first(f.Monome))
        push!(r22,first(f.Koeffizienten))
    end
    #empty!(r.Monome)
    #empty!(r.Koeffizienten)
    return PolyNom(r21,r22)
end

function Sub1(f::PolyNom{W},g::PolyNom{W},mf,mg,kf,kg) where{W}
    j=1
    k=1
    
    lg = length(g.Koeffizienten)
    lf = length(f.Koeffizienten)
    A = CircularDeque{Vec{W,Int64}}(lg+lf)
    C = CircularDeque{FieldElem}(lg+lf)
    

    while k <=lf && j <= lg
      
        x = Vgl(f.Monome[k]+mf,g.Monome[j]+mg)
        #potentiell aufpassen
        if x == 0
            push!(A,g.Monome[j]+mg)
            push!(C,g.Koeffizienten[j]*kg)
            j+=1
        elseif x==2
            if g.Koeffizienten[j]*kg+f.Koeffizienten[k]*kf != 0
                push!(C,f.Koeffizienten[k]*kf + g.Koeffizienten[j]*kg)
                push!(A,f.Monome[k]+mf)
            end
            k+=1
            j+=1
        else
            push!(A,f.Monome[k]+mf)
            push!(C,f.Koeffizienten[k]*kf)
            k+= 1
        end    
    end
    while j <=lg
        push!(A,g.Monome[j]+mg)
        push!(C,g.Koeffizienten[j]*kg)
        j+=1
    end
    while k <=lf
        push!(A,f.Monome[k]+mf)
        push!(C,f.Koeffizienten[k]*kf)
        k+=1
    end


    f2 = PolyNom(A,C)

    return f2
end





struct geobucketpol{W}
    Bucket::Vector{PolyNom{W}}
end


function addgeobucket(B::geobucketpol{W},f::PolyNom) where{W}
    #nochmal hinschauen
    log = cld(64-leading_zeros(length(f.Koeffizienten)),2)
    i=max(1,log)
    m = length(B.Bucket)
    if i <= m
        if isempty(B.Bucket[i].Koeffizienten)==false 
            f =add(f,B.Bucket[i])
        end
        while i <=m && length(f.Koeffizienten) > 4^i
            if i!=m
                if isempty(B.Bucket[i+1].Koeffizienten) == false
                    f=add(f,B.Bucket[i+1])
                end
            else
                push!(B.Bucket,f)
            end
            empty!(B.Bucket[i].Monome)
            empty!(B.Bucket[i].Koeffizienten)
            i+=1
        end
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNom(CircularDeque{Vec{W,Int64}}(1),CircularDeque{FieldElem}(1)))
    end
    B.Bucket[i] = f
    return B
end


function add(f::PolyNom{W},g::PolyNom{W})where{W}
    lf = length(f.Koeffizienten)
    lg = length(g.Koeffizienten)
    k= 1
    j= 1

    A = CircularDeque{Vec{W,Int64}}(lf+lg)
    C = CircularDeque{FieldElem}(lf+lg)
    while k <=lf && j <= lg
        
        x = Vgl(f.Monome[k],g.Monome[j])

        #potentiell aufpassen
        if x == 0
            push!(A,g.Monome[j])
            push!(C,g.Koeffizienten[j])
            j+=1
        elseif x==2
            if f.Koeffizienten[k]+g.Koeffizienten[j] != 0
                push!(C,f.Koeffizienten[k]+ g.Koeffizienten[j])
                push!(A,f.Monome[k])
            end
            k+=1
            j+=1
        else
            push!(A,f.Monome[k])
            push!(C,f.Koeffizienten[k])
            k+=1
        end   
    end
    while j <=lg
        push!(A,g.Monome[j])
        push!(C,g.Koeffizienten[j])
        j+=1
    end

    while k <=lf
        push!(A,f.Monome[k])
        push!(C,f.Koeffizienten[k])
        k+=1
    end
    h = PolyNom(A,C)
    return h
end

function Leitterm(B::geobucketpol{W}) where{W}
    m= length(B.Bucket)
    j= 0
    while true
        j= 0
        w = true
        for i=1:m
            if isempty(B.Bucket[i].Koeffizienten) == false
                if j == 0
                    j=i
                else
                    wt = Vgl(first(B.Bucket[i].Monome),first(B.Bucket[j].Monome))
                    if wt==1
                        j=i
                    elseif wt==2
                        if first(B.Bucket[i].Koeffizienten) + first(B.Bucket[j].Koeffizienten)!=0
                            B.Bucket[j].Koeffizienten.buffer[B.Bucket[j].Koeffizienten.first]+=first(B.Bucket[i].Koeffizienten)
                            popfirst!(B.Bucket[i].Koeffizienten)
                            popfirst!(B.Bucket[i].Monome)
                        else
                            popfirst!(B.Bucket[i].Koeffizienten)
                            popfirst!(B.Bucket[i].Monome)
                            popfirst!(B.Bucket[j].Koeffizienten)
                            popfirst!(B.Bucket[j].Monome)
                            w = false
                            break
                        end
                    end
                end
            end
        end
        if j==0 || w== true
            break
        end
    end
    if j== 0
        return PolyNom(CircularDeque{Vec{W,Int64}}(1),CircularDeque{FieldElem}(1)) 
    end
    #return
    h = PolyNom(CircularDeque{Vec{W,Int64}}(1),CircularDeque{FieldElem}(1))

    push!(h.Monome,popfirst!(B.Bucket[j].Monome))
    push!(h.Koeffizienten,popfirst!(B.Bucket[j].Koeffizienten))
    return h
end





function Testa()
    F2 = PolNeu(g,c,ord=ord)
    return DIV2(F2,K)
end
c = Vec{8,Int64}((0,0,0,0,0,0,0,0))
#G= [x1+x2+x3+x4+x5,x1*x2+x2*x3+x4*x5,x1*x2*x3+x2*x3*x4+x3*x4*x5,x1*x2*x3*x4+x2*x3*x4*x5,x1*x2*x3*x4*x5-1,x2^2+x2*x4] 
G= [x1 + x2 + x3 + x4 + x5, x1*x2 + x2*x3 + x3*x4 + x4*x5, x1*x2*x3 + x2*x3*x4 + x3*x4*x5, x1*x2*x3*x4 + x2*x3*x4*x5, x1*x2*x3*x4*x5 - 1, x2^2 + x2*x4 + x2*x5 - x3*x4 - x4*x5, x2*x3^2 - x2*x3*x4 + x3^2*x4, x2*x3*x4^2 - x2*x3*x4*x5 + x3*x4^2*x5, x2*x3*x4*x5^2 + 1, -x2*x3*x4*x5 - x3^3*x4 + x3^2*x4^2 - 2*x3^2*x4*x5 + 2*x3*x4^2*x5, -x3^3*x4^2 - x3^2*x4^2*x5 + 2*x3*x4^2*x5^2 + 2, -x2 - x3^2*x4^2*x5^2 - 2*x3*x4^2*x5^3 - 2*x5]
#G = [x1+x2+x3]
K= Vector{PolyNom{8}}()
for i = 1:length(G)
    push!(K,PolNeu(G[i],c,ord=ord))
end
f = x1^3 
F=  PolNeu(f,c,ord=ord)

#println(NeuPol(DIV2(F,K),PolAlg))
#println(divrem(x4*(x1*x2)^4,x1*x2+x2^2+x3)[2])
#println(divrem(f,G)[2])
#println(K)
#println(@btime DIV2(F,K))
#println(@btime divrem(f,G)[2])

function MinTest1()
    return DIV2(F,K)
end
function MinTest2()
    return divrem(g,[f1,f2])
end



#println()
#println(@btime MinTest1())
#println(@btime MinTest2())
#println(NeuPol(F,PolAlg))"""


end
