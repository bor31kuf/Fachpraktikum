module SIMD3


export PolNeu, DIV2, NeuPol, Sub1



using Oscar,SIMD,BenchmarkTools,DataStructures


mutable struct PolyNom{W}
    Monome::Deque{Vec{W,Int64}}
    Koeffizienten::Deque{FieldElem}
end

PolAlg, (x1,x2,x3,x4,x5,x6,x7) = polynomial_ring(QQ,[:x1,:x2,:x3,:x4,:x5,:x6,:x7],internal_ordering=:lex)
ord=lex(PolAlg)

function PolNeu(f,c::Vec{W,Int64};ord::MonomialOrdering=default_ordering(parent(f))) where {W}
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))



    D = PolyNom(Deque{Vec{W,Int64}}(1),Deque{FieldElem}(1)) 
    for i=1:length(A)
        push!(D.Monome,Vec{W,Int64}((sum([B[i][j]*c[j+1] for j=1:length(c)-1]),B[i]...)))
        push!(D.Koeffizienten,A[i])
    end
    return D

end  

function Vgl(a::Vec{W,Int64},b::Vec{W,Int64}) where{W}
    for i in 1:W
        if a[i] < b[i]
            return 0
        elseif a[i] > b[i]
            return 1
        end
    end
    return 2
end

function NeuPol(f::PolyNom,PolAlg) 
    a=zero(PolAlg)
    k = length(f.Monome)
    for i=1:k
        t1 = popfirst!(f.Monome)
        t2 = popfirst!(f.Koeffizienten)
        push!(f.Monome,t1)
        push!(f.Koeffizienten,t2)
        a += monomial(PolAlg,collect(Tuple(t1))[2:end])*t2
    end
    return a
end


function DIV2(f::PolyNom{W},G::Vector{PolyNom{W}}) where W
    #f2 = PolyNom(deepcopy(f.Monome),deepcopy(f.Koeffizienten))

    f2 = geobucketpol([PolyNom(Deque{Vec{W,Int64}}(1),Deque{FieldElem}(1))])
    if length(f.Monome)==0
        return f
    end
    f2 = addgeobucket(f2,f)
    A = Deque{Vec{W,Int64}}(1)
    B = Deque{FieldElem}(1)
    r = PolyNom(Deque{Vec{W,Int64}}(1),Deque{FieldElem}(1))
    D = length(G)
    LTf2 = Leitterm(f2)
    a = 0
    while length(LTf2.Monome) != 0
        w = false
        #println(length(f2.Bucket))
        """a+=1
        if a <40
            println(LTf2)
        else
            return "hi"
        end"""
        for i=1:D
            if sum((first(LTf2.Monome)>=first(G[i].Monome)))==W
                #println("nein")
                DIV1 = first(LTf2.Monome)-first(G[i].Monome)
                DIV2 = first(LTf2.Koeffizienten)/first(G[i].Koeffizienten)
                empty!(A)
                empty!(B)
                for t=1:length(G[i].Monome)
                    x = popfirst!(G[i].Monome)
                    push!(G[i].Monome,x)  
                    y= popfirst!(G[i].Koeffizienten)
                    push!(G[i].Koeffizienten,y)
                    if t!= 1
                        push!(A,x+DIV1)
                        push!(B,y*(-1*DIV2))
                    end
                end
            
                w = true
                
                if length(A)!=0
                    #println("na")
                    g = PolyNom(A,B)
                   
                    f2= addgeobucket(f2,g)
                end
                LTf2 = Leitterm(f2)
                break
            end
        end
        if w == false
            push!(r.Monome,first(LTf2.Monome))
            push!(r.Koeffizienten,first(LTf2.Koeffizienten))
            LTf2 = Leitterm(f2)
            #println(t)
        end
    end
    return r
end

function Mul1(f::PolyNom,g::Vec{W,Int64},l) where{W}
    h =PolyNom(copy(f.Monome),copy(f.Koeffizienten))
    for i=1:length(f.Monome)
        h.Monome[i] += gf2
        h.Koeffizienten[i] *= l
    end
    return h
end

function Sub1(f::PolyNom{W},g::PolyNom{W},mf,mg,kf,kg) where{W}
    j=1
    k=1
    #f2 = PolyNom(Vector{Vec{4,Int64}}(),Vector{FieldElem}())
    lg = length(g.Monome)
    lf = length(f.Monome)
    A = Deque{Vec{W,Int64}}()
    C = Deque{FieldElem}()



    while k <=lf && j <= lg
    
        fM = first(f.Monome)
        gM = first(g.Monome)
        gK= first(g.Koeffizienten)
        fK = first(f.Koeffizienten)
        x = Vgl(fM+mf,gM+mg)
        #potentiell aufpassen
        if x == 0
            push!(A,gM+mg)
            push!(C,gK*kg)
            t1 = popfirst!(g.Monome)
            t2 = popfirst!(g.Koeffizienten)
            push!(g.Monome,t1)
            push!(g.Koeffizienten,t2)
            j+=1
        elseif x==2
            if gK*kg+fK*kf != 0
                push!(C,fK*kf + gK*kg)
                push!(A,fM+mf)
            end
            t1= popfirst!(g.Monome)
            t2 = popfirst!(g.Koeffizienten)
            push!(g.Monome,t1)
            push!(g.Koeffizienten,t2)
            t1= popfirst!(f.Monome)
            t2 = popfirst!(f.Koeffizienten)
            push!(f.Monome,t1)
            push!(f.Koeffizienten,t2)
            k+=1
            j+=1
        else
            push!(A,fM+mf)
            push!(C,fK*kf)
            t1= popfirst!(f.Monome)
            t2 = popfirst!(f.Koeffizienten)
            push!(f.Monome,t1)
            push!(f.Koeffizienten,t2)
            k+= 1
        end    
    end
    while j <=lg
        gM =popfirst!(g.Monome)
        gK = popfirst!(g.Koeffizienten)
        push!(A,gM+mg)
        push!(C,gK*kg)
        push!(g.Monome,gM)
        push!(g.Koeffizienten,gK)
        j+=1
    end
    while k <=lf
        fM =popfirst!(f.Monome)
        fK = popfirst!(f.Koeffizienten)
        push!(A,fM+mf)
        push!(C,fK*kf)
        push!(f.Monome,fM)
        push!(f.Koeffizienten,fK)
        k+=1
    end


    f2 = PolyNom(A,C)

    return f2
end






struct geobucketpol{W}
    Bucket::Vector{PolyNom{W}}
end


function addgeobucket(B::geobucketpol{W},f::PolyNom) where{W}
    i=max(1,ceil(Int,log(4,length(f.Monome))))
    m = length(B.Bucket)
    if i <= m
        f =add(f,B.Bucket[i])
        while i <m && length(f.Monome) > 4^i
            f=add(f,B.Bucket[i+1])
            empty!(B.Bucket[i].Monome)
            empty!(B.Bucket[i].Koeffizienten)
            i+=1
        end
    end
    for t=m:max(m,i)-1
        push!(B.Bucket, PolyNom(Deque{Vec{W,Int64}}(1),Deque{FieldElem}(1)))
    end
    m = max(m,i)
    B.Bucket[i] = f
    return B
end


function add(f::PolyNom{W},g::PolyNom{W})where{W}
    lf = length(f.Monome)
    lg = length(g.Monome)
    k= 1
    j= 1
    #f2 = deepcopy(f)
    M= f.Monome
    N = f.Koeffizienten

    A = Deque{Vec{W,Int64}}(1)
    C = Deque{FieldElem}(1)

    while k <=lf && j <= lg
    
        fM = first(M)
        gM = first(g.Monome)
        gK= first(g.Koeffizienten)
        fK = first(N)
        x = Vgl(fM,gM)
        if x == 0
            push!(A,gM)
            push!(C,gK)
            popfirst!(g.Monome)
            popfirst!(g.Koeffizienten)
            j+=1
        elseif x==2
            if gK+fK != 0
                push!(C,fK + gK)
                push!(A,fM)
            end
            popfirst!(g.Monome)
            
            
            popfirst!(g.Koeffizienten)
            popfirst!(N)
            popfirst!(M)
            k+=1
            j+=1
        else
            push!(A,fM)
            push!(C,fK)
            popfirst!(M)
            popfirst!(N)
            k+= 1
        end    
    end
    while j <=lg
        gM =popfirst!(g.Monome)
        gK = popfirst!(g.Koeffizienten)
        push!(A,gM)
        push!(C,gK)
        j+=1
    end
    while k <=lf
        fM =popfirst!(M)
        fK = popfirst!(N)
        push!(A,fM)
        push!(C,fK)
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
            if first(isempty(B.Bucket[i].Monome)) == false
                if j == 0
                    j=i
                else
                    wt = Vgl(first(B.Bucket[i].Monome),first(B.Bucket[j].Monome))
                    if wt==1
                        j=i
                    elseif wt==2
                        if first(B.Bucket[i].Koeffizienten) + first(B.Bucket[j].Koeffizienten)!=0
                            a= first(B.Bucket[j].Koeffizienten)+ first(B.Bucket[i].Koeffizienten)
                            popfirst!(B.Bucket[j].Koeffizienten)
                            pushfirst!(B.Bucket[j].Koeffizienten,a)
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
        return PolyNom(Deque{Vec{W,Int64}}(1),Deque{FieldElem}(1)) 
    end
    #return
    h = PolyNom(Deque{Vec{W,Int64}}(1),Deque{FieldElem}(1))

    push!(h.Monome,popfirst!(B.Bucket[j].Monome))
    push!(h.Koeffizienten,popfirst!(B.Bucket[j].Koeffizienten))
    return h
end






c = Vec{8,Int64}((0,0,0,0,0,0,0,0))
G= [x1+x2+x3+x4+x5,x1*x2+x2*x3+x4*x5,x1*x2*x3+x2*x3*x4+x3*x4*x5,x1*x2*x3*x4+x2*x3*x4*x5,x1*x2*x3*x4*x5-1,x2^2+x2*x4] 
K= Vector{PolyNom{8}}()
for i = 1:length(G)
    push!(K,PolNeu(G[i],c,ord=ord))
end
f = x1^2*x2*x3+x1*x2^2*x3+x1*x2*x3^2-x2*x3*x4*x5
F=  PolNeu(f,c,ord=ord)
#println(K)
#println(F)
#println(DIV2(F,K))
#println(divrem(f,G)[2])
function Testa()
    F2 = PolNeu(f,c,ord=ord)
    return DIV2(F2,K)
end

#println(Testa())
#Testa()
#Testa()
#println(Testa())
#println(Testa())
#println(@btime Leitterm(x))
#println(x)
#PolNeu(g,c,ord=ord)
#println(x)
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
