function my_spoly(f,g)
    if f == 0 || g == 0
        return 0
    else
        m12 = leading_term(f)/(gcd(leading_monomial(f), leading_monomial(g)))
        m21 = leading_term(g)/(gcd(leading_monomial(f), leading_monomial(g)))
    return m21*f - m12*g
    end
end

function div_mit_rest(L,f)
    t = length(L)
    r = 0
    while f != 0
           i = 1
           while divides(leading_term(f), leading_term(L[i])) == (false, 0)
                  i += 1
                  if i > t
                         break
                  end
           end
           if i > t
                  r += leading_term(f)
                  f -= leading_term(f)
           else #(i <= t)
                  f -= (leading_term(f)/leading_term(L[i]))*L[i]
           end
    end
    return r
end

function my_gb(L)
    i = length(L)
    i2 = 2
    L2 = deepcopy(L)
    while true
        S3 = []
        for k=1:i
            for j=min(k+1,i2):i
                if gcd(leading_monomial(L2[i]),leading_monomial(L2[k])) != 1
                    S = my_spoly(L2[i],L2[k])
                    S2 = div_mit_rest(L2,S)
                    if S2 != 0
                        S3 = push!(S3,S2)
                    end
                end
            end
        end
        L2 = push!(L2,S3)
        if i == length(L2)
            break
        end
        i2 = i
        i = length(L2)
    end
    return L2
end

function test_ord(G1)
    for i=1:length(G1)
        my_gb(G1[i])
    end
    return
end

function test()
    G = []
    for i = 1:4
        G1 = []
        for k = 1:1
            G1 = push!(G1,gens(katsura(i)))
        end
        G = push!(G,G1)
    end
    println("Zeit mit lex")
    R1, (x1,x2) = QQ[:x1,:x2]
    R2, (x1,x2,x3) = QQ[:x1,:x2,:x3]
    R3, (x1,x2,x3,x4) = QQ[:x1,:x2,:x3,:x4]
    R4, (x1,x2,x3,x4,x5) = QQ[:x1,:x2,:x3,:x4,:x5]
    R = [R1,R2,R3,R4]
    for k=1:4
        F = G[k]
        with_ordering(R[k],lex(R[k])) do
            @time test_ord(F)
        end
    end
    println("Zeit mit degrevlex")
    for k=1:4
        F = G[k]
        with_ordering(R[k],deginvlex(R[k])) do
            @time test_ord(F)
        end
    end
    println("Korrektheit")
    K = []
    for i = 1:4
        K1 = []
        for k = 1:10
            K1 = push!(K1,gens(katsura(i)))
        end
        K = push!(K,K1)
    end
    for k=1:4
        for i=1:length(K[k])
            if my_isgb(my_gb(K[k][i])) != true
                println("Nein")
                return false
            end
        end
    end
    println("Ja")
    return
end



