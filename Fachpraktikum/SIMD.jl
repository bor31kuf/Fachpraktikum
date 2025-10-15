module SIMD1


export PolNeu, DIV2, NeuPol, Sub1



using Oscar,SIMD,BenchmarkTools


struct PolyNom{W}
    Monome::Vector{Vec{W,Int64}}
    Koeffizienten::Vector{FieldElem}
end

PolAlg, (v,x,y,z) = polynomial_ring(QQ,[:v,:x,:y,:z],internal_ordering=:lex)
ord= lex(PolAlg)

function PolNeu(f,c::Vec{W,Int64};ord::MonomialOrdering=default_ordering(parent(f))) where {W}
    A = collect(coefficients(f,ordering=ord))
    B = collect(exponents(f,ordering=ord))
    """if ord==:lex
        println("naja")
    elseif ord==deglex(gens(parent(f)))
        println("cool")
    elseif ord==degrevlex(gens(parent(f)))
        println("Hi")
        println("Nicht definiert")
    end"""
    D = PolyNom(Vector{Vec{W,Int64}}(),Vector{FieldElem}()) 
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
    for i=1:length(f.Monome)
        a += monomial(PolAlg,collect(Tuple(f.Monome[i]))[2:end])*f.Koeffizienten[i]
    end
    return a
end


function DIV2(f::PolyNom{W},G::Vector{PolyNom{W}}) where W
    f2 = PolyNom(copy(f.Monome),copy(f.Koeffizienten))
    
    #f2 = f
    #M = length(f.Monome[1])
    r = PolyNom(Vector{Vec{W,Int64}}(),Vector{FieldElem}())
    D = length(G)
    a=0
    while length(f2.Monome) != 0
        w = false
        for i=1:D
            if sum((f2.Monome[1]>=G[i].Monome[1]))==W
                DIV1 = f2.Monome[1]-G[i].Monome[1]
                DIV2 = -f2.Koeffizienten[1]/G[i].Koeffizienten[1]
                L2 = length(G[i].Monome)
                L3 = length(f2.Monome)
                j=2
                k=2
                A = Vector{Vec{W,Int64}}()
                B = Vector{FieldElem}()
                while j <=L2 && k <= L3
                    x = Vgl(f2.Monome[k],G[i].Monome[j]+DIV1)
                    if x == 0
                        push!(A,G[i].Monome[j]+DIV1)
                        push!(B,G[i].Koeffizienten[j]*DIV2)
                        j+=1
                    elseif x==2
                        if G[i].Koeffizienten[j]*DIV2+f2.Koeffizienten[k] != 0
                            push!(B,f2.Koeffizienten[k] + G[i].Koeffizienten[j]*DIV2)
                            push!(A,f2.Monome[k])
                        end
                        k+=1
                        j+=1
                    else
                        push!(A,f2.Monome[k])
                        push!(B,f2.Koeffizienten[k])
                        k+= 1
                    end    
                end
                while j<= L2
                    push!(A,G[i].Monome[j]+DIV1)
                    push!(B,G[i].Koeffizienten[j]*DIV2)
                    j+=1
                end
                while k <= L3
                    push!(A, f2.Monome[k])
                    push!(B,f2.Koeffizienten[k])
                    k+=1
                end
                f2 = PolyNom(A,B)
                w =true
                if length(A)==0
                    return r
                    #macht es halb richtig
                end
                break
            end
        end
        if w == false
            push!(r.Monome,f2.Monome[1])
            push!(r.Koeffizienten,f2.Koeffizienten[1])
            deleteat!(f2.Monome,1)
            deleteat!(f2.Koeffizienten,1)
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
    j=2
    k=2
    #f2 = PolyNom(Vector{Vec{4,Int64}}(),Vector{FieldElem}())
    a = length(g.Monome)
    b = length(f.Monome)
    A = Vector{Vec{W,Int64}}()
    B = Vector{FieldElem}()

    while j <=a&& k <= b
        x = Vgl(f.Monome[k]+mf,g.Monome[j]+mg)
        if x == 0
            push!(A,g.Monome[j]+mg)
            
            push!(B,-g.Koeffizienten[j]*kg)
            j+=1
        elseif x==2
            if -g.Koeffizienten[j]*kg+f.Koeffizienten[k]*kf != 0
                push!(B,f.Koeffizienten[k]*kf -g.Koeffizienten[j]*kg)
                push!(A,f.Monome[k]+mf)
            end
            k+=1
            j+=1
        else
            push!(A,f.Monome[k]+mf)
            push!(B,f.Koeffizienten[k]*kf)
            k+= 1
        end    
    end
    while j<=a
        push!(A,g.Monome[j]+mg)
        push!(B,-g.Koeffizienten[j]*kg)
        j+=1
    end
    while k <= b
        push!(A,f.Monome[k]+mf)
        push!(B,f.Koeffizienten[k]*kf)
        k+=1
    end
    f2 = PolyNom(A,B)

    return f2
end










"""f1 =x*y
f2 =x^3+z
g = x^12+x*z^3*y+x^4*z+z^13+z^7*y+y^2*x+v
c = Vec{5,Int64}((0,0,0,0,0))
F = PolNeu(g,c,ord=ord)
K = Vector{PolyNom{5}}()
push!(K,PolNeu(f1,c,ord=ord))
push!(K,PolNeu(f2,c,ord=ord))

println(F)
function MinTest1()
    return DIV2(F,K)
end
function MinTest2()
    return divrem(g,[f1,f2])
end

println(NeuPol(MinTest1(),PolAlg))
println(NeuPol(F,PolAlg))
println(K)
println(MinTest2())

println(@btime MinTest1())
println(@btime MinTest2())"""
end
