## Blatt 6, Aufgabe 1
## Joris Wirsen, Helena Belzl

using Oscar

function my_spoly(f,g)
    if f == 0 || g == 0
        return 0
    else
        m12 = leading_term(f)/(gcd(leading_monomial(f), leading_monomial(g)))
        m21 = leading_term(g)/(gcd(leading_monomial(f), leading_monomial(g)))
    return m21*f - m12*g
    end
end
