using Revise,Oscar,BenchmarkTools

R, (x,y,z) = polynomial_ring(QQ,[:x,:y,:z],internal_ordering=:degrevlex) 


function DivModular(f,G,ordering::MonomialOrdering=default_ordering(parent(f)))
    p = 3
    Ziel = 0*f
    z = 1
    QZiel = [parent(f)(0) for i in 1:length(G)]

    while true
        fPrim = QQdinFQM(f,p)
        GPrim = [parent(fPrim)(0) for i in 1:length(G)]
        for i in 1:length(G)
            GPrim[i] = QQdinFQM(G[i],p)
        end
        Loes = DivRestFam(fPrim,GPrim,ordering,true)
        FDiv = Loes[1]
        QDiv = Loes[2]
        a = 0
        if FDiv ==0
            return 0*f
        end
        if Ziel == 0
            Ziel = combine(FDiv,Ziel,p,z)
            a = Ziel
            for j in 1:length(G)
                QZiel[j] = combine(QDiv[j],QZiel[j],p,z)
                a+=QZiel[j]*G[j]
            end
            #println(p)
        elseif monomial(parent(f),exponent_vector(leading_monomial(FDiv,ordering= ordering),1)) == leading_monomial(Ziel,ordering=ordering)
            Ziel = combine(FDiv,Ziel,p,z)
            a = Ziel
            for j in 1:length(G)
                QZiel[j] = combine(QDiv[j],QZiel[j],p,z)
                a+=QZiel[j]*G[j]
            end
            #println(p)
        elseif leading_monomial(monomial(parent(f),exponent_vector(leading_monomial(FDiv,ordering = ordering),1))+Ziel,ordering=ordering) != leading_monomial(Ziel,ordering=ordering)
            Ziel = combine(FDiv,0*f,p,1)
            a = Ziel
            for j in 1:length(G)
                QZiel[j] = combine(QDiv[j],0*f,p,1)
                a+=QZiel[j]*G[j]
            end
            z = 1
            #println(p)
        end
        if a == f
            println(p) 
            break
        end
        z *= p
        p = next_prime(p)

        #println(f)
        #println(Ziel)
        if p == 53
            break
        end
    end
    return Ziel
end

G = [x+230y+z,13y+1+x^2+x+z]
f = 5y^2+7z+x^2+x^7
function combine(f1, f2,p,z)
    coeffs1 = Dict(m => c for (c, m) in zip(coefficients(f1), monomials(f1)))
    coeffs2 = Dict(m => c for (c, m) in zip(coefficients(f2), monomials(f2)))
    monoms = union(keys(coeffs1), keys(coeffs2))
    f=0
    R = parent(f2)
    for m in monoms
        a = haskey(coeffs1, m) ? coeffs1[m] : GF(2)(0)
        b = haskey(coeffs2, m) ? coeffs2[m] : 0
        c = crt(Int(lift(ZZ,a)),p,Int(b),z)
        if c> p*z/2
            c-=p*z
        end
        e=  exponent_vector(m,1)
        f += c*monomial(R,e)
    end
    if f==0
        f = parent(f2)(0)
    end
    return f
end

function QQdinFQM(f,p)
    ng = ngens(parent(f))
    varNam = ["xi" for i in 1:ng]
    R, vars = polynomial_ring(GF(p),varNam)
    f2 = R(0)
    for (c,m) in zip(coefficients(f),monomials(f))
        coeff = GF(p)(c)
        e = exponent_vector(m,1)
        f2+= coeff*monomial(R,e)
    end
    return f2
end

function Test()
    println(DivModular(f,G))
    A = DivRestFam(f,G,lex(R),true)
    println(A[1])
    #println(A[1]+A[2][1]*G[1])
end
Test()





