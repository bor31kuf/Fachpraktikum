## Blatt 6, Aufgabe 2
## Joris Wirsen, Helena Belzl

using Oscar

include("exercise61.jl")

R, (x,y,z) = QQ[:x,:y,:z]

# Division mit Rest
function div_mit_rest(A,L,f)
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

# div_mit_rest([x^2*y^2 - x^2*y + x, x^2*y + y^2], -x^2*y + x - y^3)


#G = [x^2*y^2 - x^2*y + x, x^2*y + y^2]
#G = [x^2*y^2 - x^2*y + x, x^2*y + y^2, x-y^3+y^2, y^8 - 2*y^7 + y^6 + y^3, y^7-2*y^6+y^5+y^2]
#my_isgb(G)