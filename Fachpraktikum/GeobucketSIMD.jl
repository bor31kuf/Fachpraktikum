module SIMD2

export PolNeu, DIV2, Sub1, NeuPol


using Oscar,SIMD,BenchmarkTools

struct PolyNom
    Monome::Array{Vec{4,Int64}}
    Koeffizienten::Array{QQFieldElem}
end



PolAlg, (x,y,z) = polynomial_ring(QQ,[:x,:y,:z],internal_ordering=:lex)

function PolNeu(f,c::Vec{4,Int64})
    o= lex([x,y,z])
    A = collect(coefficients(f,ordering=o))
    B = collect(exponents(f,ordering=o))
    D = PolyNom([],[])
    for i=1:length(A)
        push!(D.Monome,Vec{4,Int64}((sum([B[i][j]*c[j+1] for j=1:length(c)-1]),B[i][1],B[i][2],B[i][3])))
        push!(D.Koeffizienten,A[i])
    end
    return D
end  

function Vgl(a::Vec{4,Int64},b::Vec{4,Int64})
    for i =1:4
        if a[i] < b[i]
            return 0
        elseif a[i] > b[i]
            return 1
        end
    end
    return 2
end

function NeuPol(f::PolyNom) 
    a = x*0
    for i=1:length(f.Monome)
        a += monomial(PolAlg,collect(Tuple(f.Monome[i]))[2:end])*f.Koeffizienten[i]
    end
    return a
end

function first(f,G::Array{PolyNom},A,B,C)
    a = 0
    x = Vec{4,Int64}((0,0,0,0))
    if A[1] != 0
        a = 1
        x = f.Monome[A[1]]
    end
    for t=1:length(B)
        if A[t+1]!=0
            w = Vgl(G[B[t]].Monome[A[t+1]]+C[t],x) 
            if w == 1
                a= t+1
                x = G[B[t]].Monome[A[t+1]]
            end
        end
    end
    return a
end


function DIV2(f::PolyNom,G::Array{PolyNom})
    #f2 = PolyNom(copy(f.Monome),copy(f.Koeffizienten))
    ABBA = [1]#wo in den PolyNomen wir sind
    BAAB = []#welche PolyNome wir haben
    SAAB = []#Was wir drauf addieren
    DAAB = []#multiplizieren
    r = PolyNom(Vec{4,Int64}[],QQFieldElem[])
    if length(f.Monome) == 0
        return r
    end
    firste = first(f,G,ABBA,BAAB,SAAB)
    M = f.Monome[1]
    K = f.Koeffizienten[1]
    Z = length(G)
    P = Z
    while true
        w = false
        L =P
        for i=1:L
            if sum(M>=G[i].Monome[1])==4
                #println(ABBA)
                #println(BAAB)
                DIV1 = M-G[i].Monome[1]
                DIV2 = -K/G[i].Koeffizienten[1]
                push!(BAAB,i)
                push!(SAAB,DIV1)
                push!(DAAB,DIV2)
                if length(G[i].Monome)>1
                    push!(ABBA,2)
                else
                    push!(ABBA,0)
                end
                w = true

                if firste == 1
                    if length(f.Monome) > ABBA[firste]
                        ABBA[firste] +=1
                    else
                        ABBA[firste] = 0
                    end
                else
                    if length(G[BAAB[firste-1]].Monome) >ABBA[firste]
                        ABBA[firste] +=1
                    else
                        ABBA[firste] = 0
                    end    
                end
                firste = first(f,G,ABBA,BAAB,SAAB)
                if firste == 0
                    return r
                end
                if firste == 1
                    M  = f.Monome[ABBA[1]]
                    K = f.Koeffizienten[ABBA[1]]
                    
                else
                    #println(ABBA[firste])
                    M = G[BAAB[firste-1]].Monome[ABBA[firste]]+SAAB[firste-1]
                    K = G[BAAB[firste-1]].Koeffizienten[ABBA[firste]]*DAAB[firste-1]
                end
                P =i
                #break
            end
        end
        if w == false
            push!(r.Monome,M)
            push!(r.Koeffizienten,K)
            P = Z
            if firste == 1
                if length(f.Monome) > ABBA[firste]
                    ABBA[firste] +=1
                else
                    ABBA[firste] = 0
                end
            else
                if length(G[BAAB[firste-1]].Monome) >ABBA[firste]
                    ABBA[firste] +=1
                else
                    ABBA[firste] = 0
                end    
            end
            firste = first(f,G,ABBA,BAAB,SAAB)
            if firste == 0
                return r
            end
            if firste == 1
                M  = f.Monome[ABBA[1]]
                K = f.Koeffizienten[ABBA[1]]
            else
                M = G[BAAB[firste-1]].Monome[ABBA[firste]]+SAAB[firste-1]
                K = G[BAAB[firste-1]].Koeffizienten[ABBA[firste]]*DAAB[firste-1]
            end
        end

    end
end

function Mul1(f::PolyNom,g::Vec{4,Int64},l)
    h =PolyNom(copy(f.Monome),copy(f.Koeffizienten))
    for i=1:length(f.Monome)
        h.Monome[i] += g
        h.Koeffizienten[i] *= l
    end
    return h
end

function Sub1(f::PolyNom,g::PolyNom,mf,mg,kf,kg)
    j=2
    k=2
    f2 = PolyNom(Vec{4,Int64}[],QQFieldElem[])
    #A = Vec{4,Int64}[]
    #B = QQFieldElem[]
    a = length(g.Monome)
    b = length(f.Monome)
    while j <=a&& k <= b
        x = Vgl(f.Monome[k]+mf,g.Monome[j]+mg)
        if x == 0
            push!(f2.Monome,g.Monome[j]+mg)
            
            push!(f2.Koeffizienten,-g.Koeffizienten[j]*kg)
            j+=1
        elseif x==2
            if g.Koeffizienten[j]*kg+f.Koeffizienten[k]*kf != 0
                push!(f2.Koeffizienten,f.Koeffizienten[k]*kf +g.Koeffizienten[j]*kg)
                push!(f2.Monome,f.Monome[k]+mf)
            end
            j+=1
        else
            push!(f2.Monome,f.Monome[k]+mf)
            push!(f2.Koeffizienten,f.Koeffizienten[k]*kf)
            k+= 1
        end    
    end
    while j<=a
        push!(f2.Monome,g.Monome[j]+mg)
        push!(f2.Koeffizienten,-g.Koeffizienten[j]*kg)
        j+=1
    end
    while k <= b
        push!(f2.Monome,f.Monome[k]+mf)
        push!(f2.Koeffizienten,f.Koeffizienten[k]*kf)
        k+=1
    end
    #f2 = PolyNom(A,B)

    return f2
end


f1 =x*y
f2 =x^3+z
f3 = y*z^2+x
g = x^12+x*z^3*y+x^4*z+z^13+z^7*y+y^2*x
c = Vec{4,Int64}((0,0,0,0))
F = PolNeu(g,c)
K = PolyNom[PolNeu(f2,c),PolNeu(f1,c),PolNeu(f3,c)]
function MinTest1()
    return NeuPol(DIV2(F,K))
end
function MinTest2()
    return divrem(g,[f1,f2,f3])
end

#DIV2(F,K)
end
