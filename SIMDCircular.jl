module SIMD4


export PolNeu, DIV2, NeuPol, Sub1




using Oscar,SIMD,BenchmarkTools,DataStructures


mutable struct PolyNom{W}
    Monome::CircularDeque{Vec{W,Int64}}
    Koeffizienten::CircularDeque{FieldElem}
end

PolAlg, (x1,x2,x3,x4,x5,x6,x7) = polynomial_ring(QQ,[:x1,:x2,:x3,:x4,:x5,:x6,:x7],internal_ordering=:lex)
ord=lex(PolAlg)

function PolNeu(f;ord::MonomialOrdering=default_ordering(parent(f)))
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    L = length(B)
    W= length(gens(parent(f)))+1
    D = PolyNom(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L)) 
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:lex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((0,B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdeglex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:deglex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),B[i]...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    if typeof(ord.o) ==Oscar.Orderings.WSymbOrdering{:wdegrevlex}
        c = ord.o.weights
        for i=1:length(A)
            push!(D.Monome,Vec{W,Int64}((sum(c[j]*B[i][j] for j=1:W-1),reverse(B[i])...)))
            push!(D.Koeffizienten,A[i])
        end
        return D
    end
    #return D

end  

function Vgl(a::Vec{W,Int64},b::Vec{W,Int64}) where{W}
    @inbounds for i in 1:W
        if a[i] != b[i]
            return a[i] > b[i]
        end
    end
    return 2
end

function NeuPol(f,PolAlg;ord=default_ordering(PolAlg))
    a=zero(PolAlg)
    k = length(f.Monome)
    if typeof(ord.o) == Oscar.Orderings.WSymbOrdering{:wdegrevlex} ||  typeof(ord.o) == Oscar.Orderings.SymbOrdering{:degrevlex}
        for i=1:k
            a += monomial(PolAlg,reverse(collect(Tuple(f.Monome[i]))[2:end]))*f.Koeffizienten[i]
        end
    else 
        for i=1:k
            a += monomial(PolAlg,collect(Tuple(f.Monome[i]))[2:end])*f.Koeffizienten[i]
        end
    end
    return a
end


function DIV2(f::PolyNom{W},G::Vector{PolyNom{W}}) where W
    #f2 = PolyNom(deepcopy(f.Monome),deepcopy(f.Koeffizienten))
    L = length(f.Koeffizienten)
    if L==0
        return f
    end
    f2 = geobucketpol([PolyNom(CircularDeque{Vec{W,Int64}}(8),CircularDeque{FieldElem}(8))])
    f2 =addgeobucket(f2,f)
    LTf2M= first(f.Monome)
    LTf2K =first(f.Koeffizienten)
    r = PolyNom(CircularDeque{Vec{W,Int64}}(L),CircularDeque{FieldElem}(L))
    D = length(G)
    while true
        w = false
        for i=1:D
            if sum(LTf2M>=first(G[i].Monome))==W
                
                DIV1 =LTf2M-first(G[i].Monome)
                DIV2 = -LTf2K/first(G[i].Koeffizienten)
                L2 = length(G[i].Koeffizienten)
                w = true
                if L2!=1
                    f2= addgeobucket(f2,G[i],DIV1,DIV2)
                end
            
                LTf2M,LTf2K = Leitterm(f2)
                if LTf2K == 0
                    return r
                end
                break
                
               
                
            end
        end
        if w == false
            r= pushing(r,LTf2M,LTf2K)
            LTf2M,LTf2K = Leitterm(f2)
            if LTf2K == 0
                if length(r.Koeffizienten) == 0
                    return r
                end
                if first(r.Koeffizienten) == 1
                    return r
                else 
                    tt = first(r.Koeffizienten)
                    for i =1:length(r.Koeffizienten)
                        r.Koeffizienten.buffer[i] /= tt
                    end
                    return r
                end
            end
        end
    end
end


function pushing(r::PolyNom{W},LTf2M,LTf2K) where{W}
    if capacity(r.Koeffizienten) > length(r.Koeffizienten)
        push!(r.Monome,LTf2M)
        push!(r.Koeffizienten,LTf2K)
        return r
    else
        r21 = CircularDeque{Vec{W,Int64}}(2*capacity(r.Monome))
        r22 = CircularDeque{FieldElem}(2*capacity(r.Monome))
        for i=1:length(r.Koeffizienten)
            push!(r21,r.Monome[i])
            push!(r22,r.Koeffizienten[i]) 
        end
        push!(r21,LTf2M)
        push!(r22,LTf2K)
    end

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


function addgeobucket(B::geobucketpol{W},f::PolyNom,DIV1 = Vec{W,Int64}(ntuple(i-> 0,W)),DIV2 =QQFieldElem(1)) where{W}
    #nochmal hinschauen
    log = cld(64-leading_zeros(length(f.Koeffizienten)),2)
    i=max(1,log)
    m = length(B.Bucket)
    if i <= m
        add2(B.Bucket[i],f,DIV1,DIV2)
        while i <=m && length(B.Bucket[i].Koeffizienten) > 4^i
            if i!=m
                add(B.Bucket[i+1],B.Bucket[i])
                empty!(B.Bucket[i].Koeffizienten)
                empty!(B.Bucket[i].Monome)
            else
                push!(B.Bucket,PolyNom(CircularDeque{Vec{W,Int64}}(2*4^(m+1)),CircularDeque{FieldElem}(2*4^(m+1))))
                add(B.Bucket[m+1],B.Bucket[m])
                empty!(B.Bucket[m].Koeffizienten)
                empty!(B.Bucket[m].Monome)
            end
            i+=1
        end
        return B
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNom(CircularDeque{Vec{W,Int64}}(2*4^(t+1)),CircularDeque{FieldElem}(2*4^(t+1))))
    end
    add2(B.Bucket[i],f,DIV1,DIV2) 
    return B
end

function add2(f::PolyNom{W},g::PolyNom{W},DIV1,DIV2)where{W}
    lf = length(f.Koeffizienten)
    lg = length(g.Koeffizienten)
    k= 1
    j= 2
    A =f.Koeffizienten.last
    t = 0
    while k <=lf && j <= lg
        t+=1
        x = Vgl(f.Monome[k],g.Monome[j]+DIV1)
        #potentiell aufpassen
        if x == 0
            push!(f.Monome,g.Monome[j]+DIV1)
            push!(f.Koeffizienten,g.Koeffizienten[j]*DIV2)
            j+=1
        elseif x==2
            if f.Koeffizienten[k]+g.Koeffizienten[j]*DIV2 != 0
                push!(f.Koeffizienten,f.Koeffizienten[k]+ g.Koeffizienten[j]*DIV2)
                push!(f.Monome,f.Monome[k])
            else
                f.Koeffizienten.n +=1
                f.Monome.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.Monome,f.Monome[k])
            push!(f.Koeffizienten,f.Koeffizienten[k])
            k+=1
        end
        f.Koeffizienten.n -=1
        f.Monome.n -=1   
    end
    while j <=lg
        push!(f.Monome,g.Monome[j]+DIV1)
        push!(f.Koeffizienten,g.Koeffizienten[j]*DIV2)
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.Monome,f.Monome[k])
        push!(f.Koeffizienten,f.Koeffizienten[k])
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.Koeffizienten.capacity
        f.Monome.first = A+1
        f.Koeffizienten.first = A+1 
    else
        f.Monome.first = 1
        f.Koeffizienten.first = 1
    end
    f.Monome.n  = t
    f.Koeffizienten.n = t
    
    return 
    
end

function add(f::PolyNom{W},g::PolyNom{W})where{W}
    lf = length(f.Koeffizienten)
    lg = length(g.Koeffizienten)
    k= 1
    j=1
    A =f.Koeffizienten.last
    t = 0
    while k <=lf && j <= lg
        t+=1
        x = Vgl(f.Monome[k],g.Monome[j])
 
        if x == 0
            push!(f.Monome,g.Monome[j])
            push!(f.Koeffizienten,g.Koeffizienten[j])
            j+=1
        elseif x==2
            if f.Koeffizienten[k]+g.Koeffizienten[j] != 0
                push!(f.Koeffizienten,f.Koeffizienten[k]+ g.Koeffizienten[j])
                push!(f.Monome,f.Monome[k])
            else
                f.Koeffizienten.n +=1
                f.Monome.n +=1
                t-=1
            end
            k+=1
            j+=1
        else
            push!(f.Monome,f.Monome[k])
            push!(f.Koeffizienten,f.Koeffizienten[k])
            k+=1
        end
        f.Koeffizienten.n -=1
        f.Monome.n -=1   
    end
    while j <=lg
        push!(f.Monome,g.Monome[j])
        push!(f.Koeffizienten,g.Koeffizienten[j])
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        j+=1
    end

    while k <=lf
        push!(f.Monome,f.Monome[k])
        push!(f.Koeffizienten,f.Koeffizienten[k])
        f.Koeffizienten.n -=1
        f.Monome.n -=1
        t+=1
        k+=1
    end
    if A+1 <= f.Koeffizienten.capacity
        f.Monome.first = A+1
        f.Koeffizienten.first = A+1 
    else
        f.Monome.first = 1
        f.Koeffizienten.first = 1
    end
    f.Monome.n  = t
    f.Koeffizienten.n = t
    return 
    
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
        return 0,0
    end
    return popfirst!(B.Bucket[j].Monome),popfirst!(B.Bucket[j].Koeffizienten) 
end





function Testa()
    F2 = PolNeu(g,c,ord=ord)
    return DIV2(F2,K)
end
c = Vec{8,Int64}((0,0,0,0,0,0,0,0))
G = [5//2*x1^2*x2*x3 - 4*x1*x3^2, 5//3*x1*x2*x3^3, 4//3*x1^2*x2^2*x3^3 - 1//3*x1*x2, 2//5*x1*x2^2*x3^3, -3//5*x2^2*x3, -1//2*x1^2*x3^2, 1//3*x1*x2*x3^3, -2*x1*x2^4*x3]
K= Vector{PolyNom{8}}()
#for i = 1:length(G)
 #   push!(K,PolNeu(G[i],c,ord=ord))
#end
#f = 1//5*x1^11*x2^8*x3^7 + 3//4*x1^9*x2^8*x3^3 + 23//5*x1^9*x2^5*x3^8 + 9//29*x1^8*x2^10*x3^8 + 8//9*x1^8*x2^9*x3^5 + 9//19*x1^6*x2^11*x3^11 + 13//12*x1^5*x2^9*x3^6 + 1//2*x1^3*x2^7*x3^10
#F=  PolNeu(f,c,ord=ord)


#println(DIV2(F,K),PolAlg)
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
