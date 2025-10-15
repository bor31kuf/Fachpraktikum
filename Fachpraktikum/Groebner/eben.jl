using Hecke
#using Nemo

#R, b = finite_field(5)
R, b = finite_field(23)
Rx, x = R[:x]
f = x^4+1
#println(f^2)
println(factor(x^4+1))