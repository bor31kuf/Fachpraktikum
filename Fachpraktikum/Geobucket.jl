using Revise,Oscar,BenchmarkTools

R, (x,y,z) = polynomial_ring(QQ,[:x,:y,:z])

#geobucket
function NewReduce(f,G,gMax=[],gMaxVec=[],ordering=default_ordering(parent(f)))
    if f == 0
        return f
    end
    l=1
    while true
        Ah = false
        H = collect(terms(f,ordering=ordering))
        F = AbstractAlgebra.Generic.geobucket(parent(f))
        Base.push!(F,f)
        for i = 1: length(G)
            k = exponent_vector(H[l],1)
            if G[i] != 0
                if Vektorvergleich2(k,gMaxVec[i])
                    ge = -H[l]/gMax[i]*G[i]
                    Base.push!(F,ge)
                    l+=1
                    Ah = true
                    if l>length(H)
                        break
                    end
                end
            end
        end
        f = finish(F)
        if !Ah && l>= length(H)
            return f*leading_monomial(f,ordering=ordering)/leading_term(f,ordering=ordering)
        elseif !Ah 
            l+=1
        else 
            l=1
        end
        if f == 0
            return f
        end
    end
end

function Vektorvergleich2(a::Array,b::Array)
    for i= 1:length(a)
        if a[i] < b[i] 
            return false
        end
    end
    return true
end

f = x^2+y^3
g = y^4
h = x*y
o = 3z^7*x^3*y+13y^5*x^3+143x+x^8
G1 = [f,g,h]
gMaxVec = [[0,3,0],[0,4,0],[1,1,0]]
gMax = [y^3,y^4,x*y]

function Test()
    return NewReduce(o,G1,gMax,gMaxVec)
end

println(Test())
