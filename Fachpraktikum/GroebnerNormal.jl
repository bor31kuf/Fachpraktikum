using SIMD,Oscar,Revise,BenchmarkTools


PolAlg, (x1,x2,x3,x4,x5,x6,x7,x8) = polynomial_ring(QQ,[:x1,:x2,:x3,:x4,:x5,:x6,:x7,:x8],internal_ordering=:lex)
ord = lex(PolAlg)

function SPoly(f,g,ord)
    f1 = leading_term(f,ordering=ord)
    g1 = leading_term(g,ordering=ord)
    kgv = lcm(f1,g1)
    #println(cmp(f1,g1))     
    return kgv/f1*f-kgv/g1*g
end

function Buchberger2(G,ord) 
    L = length(G)
    Queue = pairs(L)
    Bits = trues(Int(L*(L-1)/2))
    k = 1
    while k <= length(Bits)
        if Bits[k]
            
            a = Queue[k][1]
            b = Queue[k][2]
            Sij = SPoly(G[a],G[b],ord)
            #println(SIMD4.NeuPol(Sij,PolAlg))
            S = zero(parent(G[1]))
            with_ordering(parent(G[1]),ord) do
                S = divrem(Sij,G)[2]
            end
           
            if S!=0 
                push!(G,S)
                Queue, Bits = QUEUE(G,Queue,Bits,k,ord)
            end
            #break
        end
        k+=1
    end
    return G
end


function pairs(n::Int)::Vector{Tuple{Int,Int}}
    t = Vector{Tuple{Int,Int}}()
    for i=1:n
        for j=i+1:n
            push!(t,(i,j))
        end
    end
    return t
end

function my_isgb(G,ord)
    t = length(G)
    for i = 1:t-1
        for j=1
            with_ordering(parent(G[1]),ord) do
                if divrem(SPoly(G[i],G[j],ord),G) != 0
                    return false
                end
            end
        end
    end
    return true
end

function reduce_groebner(G,ord)
    i = 1
    w =false
    L  = length(G)
    m=0
    while m <10
        w = true
        m+=1
        while i <= length(G)
            G2 = deepcopy(G)
            deleteat!(G2,i)
            a = zero(parent(f))
            with_ordering(parent(G[1]),ord) do
                a = divrem(G[i],G2)[2]
                if a==0
                    deleteat!(G,i)
                    w = false 
                #elseif leading_monomial(G[i],ordering=ord) != leading_monomial(a,ordering=ord)
                #   w = false
                #  G[i]=a
                # i+=1
                else

                    G[i] = a
                    i+=1
                end
            end
        end
    end
    return G
end

function Groebner(G,ord=default_ordering(parent(G[1])))
    L = length(G)
    G = Buchberger2(G,ord)
    #G = G[1:L]
    #println(G)
    #H = reduce_groebner(G,ord)
    return G
end

function QUEUE(G,Pairs,Bits,k,ord)
    #nach caramba.inria.fr/sem-slides/201409111030
    #EDER,Faugere,Martani,Perry,Roune
    #Seminar of the CARAMEL Team in Nancy, France
    #11.9.2014
    PolAlg
    h = G[length(G)]
    c = length(Bits)
    for i=k+1:length(Bits)
        if Bits[i]
            f=G[Pairs[i][1]]
            g=G[Pairs[i][2]]
            r = lcm(leading_monomial(f,ordering=ord),leading_monomial(g,ordering=ord))
            w1 = divides(r,leading_monomial(h,ordering=ord))
            #println(w1)
            if w1[1] == true && cmp(leading_monomial(h,ordering=ord),leading_monomial(f,ordering=ord)) != 0 && cmp(leading_monomial(h,ordering=ord),leading_monomial(g,ordering=ord)) != 0
                #println("hi")
                Bits[i] == 0
            end
        end
    end

    for i=1:length(G)-1
        push!(Pairs,(i,length(G)))
        w = cmp(lcm(leading_monomial(G[i],ordering=ord),leading_monomial(G[length(G)],ordering=ord)),leading_monomial(G[i],ordering=ord)*leading_monomial(G[length(G)],ordering=ord))
        if w == 0
            #println("Jo")
            push!(Bits,false)
        else
            push!(Bits,true)
        end
    end

    for i=1:length(G)-1
        if Bits[c+i]
            for j=i+1:length(G)-1
                if Bits[c+j]
                    r1 = lcm(leading_monomial(G[length(G)],ordering=ord),leading_monomial(G[i],ordering=ord))
                    r2 = lcm(leading_monomial(G[length(G)],ordering=ord),leading_monomial(G[j],ordering=ord))

                    w1 = divides(r1,r2)
                    w2 = divides(r2,r1)

                    if w1[1] ==true
                        Bits[c+i] = false
                        break
                    elseif w2[1] ==true
                        Bits[c+j] = false
                    end
                end
            end
        end
    end
    #println(sum(Bits))
    return Pairs, Bits

end
G= [x1+x2+x3+x4+x5,x1*x2+x2*x3+x3*x4+x4*x5,x1*x2*x3+x2*x3*x4+x5*x3*x4,x1*x2*x3*x4+x5*x2*x3*x4,x1*x2*x3*x4*x5-1]

function Test()
    G= [x1+x2+x3+x4+x5,x1*x2+x2*x3+x3*x4+x4*x5,x1*x2*x3+x2*x3*x4+x5*x3*x4,x1*x2*x3*x4+x5*x2*x3*x4,x1*x2*x3*x4*x5-1]
    return Groebner(G,ord)
end    
#println(G)
println(@btime Test())
#println(@btime collect(groebner_basis(ideal(G),ordering=ord,complete_reduction=true)))
#println(SPoly(f,g))