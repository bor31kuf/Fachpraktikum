using Revise,Oscar,BenchmarkTools
include("Version2.jl")
include("Geobucket.jl")

PolAlg, (x,y,z) = polynomial_ring(QQ,[:x,:y,:z],internal_ordering=:lex) 

#Polynom speichern durch Zahlen
function PolNeu(f,c)
    #besser Speichern
    A = collect(coefficients(f))
    B = collect(exponents(f))
    D = QQFieldElem[]
    Dim = dim(parent(f))
    for i=1:length(A)
        push!(D,sum([B[i][j]*c[j] for j=1:length(c)])*10^(2*Dim)+ sum([B[i][l]*10^(2l-2) for l=1:Dim]))
        push!(D,A[i])
        if D[2i] == 0
            deleteat!(D,i)
        end
    end
    """sort!(D,rev=true)
    E = QQFieldElem[]
    for i=1:length(D)
        push!(E,D[i][1])
        push!(E,D[i][2])
    end"""
    return D
end  


f = x^2*z+y^3*z
g = y^4*z
h = x*y
o = 3z^7*x^3*y+13y^5*x^3+143x+x^8
G1 = [f,g,h]
l = PolNeu(f,[1,1,1])
m = PolNeu(g,[1,1,1])
p = PolNeu(h,[1,1,1])
q = PolNeu(o,[1,1,1])
G= [l,m,p]

function Test()
    return DIV2(q,G,3) 
end

function DIV2(f,G,Dim)
    
    f2 = deepcopy(f)
    r = QQFieldElem[]
    L = length(G)
    o=1
    while f2 != []
        while o!= L
            if G[o][1]>f2[1]
                o+=1
            else
                break
            end
        end
        if o==L
            for i =1:length(f2)
                push!(r,f2[i])
            end
            return r
        end
        A = [mod(f2[1],10^(2*(Dim-i+1))) for i=1:Dim]
        w = false
        for i=o:L
            ww = true
            B = mod(G[i][1],10^(2*Dim))
            #maybe speichern
            for j=1:Dim
                if A[j]<B
                    ww=false
                    break
                else
                    B = mod(B,10^(2*(Dim-j)))
                end
            end
            if ww == true
                DIV1 = f2[1]-G[i][1]
                DIV2 = -f2[2]/G[i][2]
                L2 = length(G[i])
                j=1
                k=1
                while j <=L2 && k <= length(f2)
                    if f2[k] < G[i][j]+DIV1
                        insert!(f2,k,G[i][j]+DIV1)
                        insert!(f2,k+1,G[i][j+1]*DIV2)
                        j+=2
                        k+=2
                    elseif f2[k] == G[i][j]+DIV1
                        if G[i][j+1]*DIV2+f2[k+1] != 0
                            f2[k] += G[i][j+1]*DIV2
                            k+=2
                        else
                            deleteat!(f2,k)
                            deleteat!(f2,k)
                        end
                        j+=2
                    else
                        k+= 2
                    end    
                end
                while j<=length(G[i])
                    push!(f2,G[i][j]+DIV1)
                    push!(f2,G[i][j+1]*DIV2)
                    j+=2
                end
                w =true
                break
            end 
        end
        if w == false
            push!(r,f2[1])
            push!(r,f2[2])
            deleteat!(f2,1)
            deleteat!(f2,1)
        end
    end
    return r
end


A = [[rand(PolAlg,-1:2,0:4,-10:6) for j=1:12] for i=1:30]
B = []
while length(B) != 30    x^8 + 143*x
    x= rand(PolAlg,-1:8,3:7,1:10)
    if x!= []
        push!(B,x)
    end
end
Ane = [[PolNeu(A[i][j],[0,0,0]) for j=1:12] for i = 1:30]
Bne = [PolNeu(B[i],[0,0,0]) for i= 1:30]
j = 1
for i=1:30
    j=1
    while j<= length(Ane[i])
        if Ane[i][j] == []
            deleteat!(Ane[i],j)
            deleteat!(A[i],j)
        else
            j+=1
        end
    end
end


function Test1()
    for i=1:30
        DIV2(Bne[i],Ane[i],3)
    end
end

function Test2()
    for i=1:30
        divrem(B[i],A[i])
    end
end
function Test3()
    for i=1:30
        DivRestFam(B[i],A[i])
    end
end

#Test1()
#Test2()
#Test3()