#Obliczenia naukowe lista 2, zadanie 4
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

using Polynomials

coefficients_array = [1, -210.0, 20615.0,-1256850.0,
      53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,
      11310276995381.0, -135585182899530.0,
      1307535010540395.0,     -10142299865511450.0,
      63030812099294896.0,     -311333643161390640.0,
      1206647803780373360.0,     -3599979517947607200.0,
      8037811822645051776.0,      -12870931245150988800.0,
      13803759753640704000.0,      -8752948036761600000.0,
      2432902008176640000.0]

roots_array = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]

new_coefficients_array = [1, -210.0 - 2.0^-23, 20615.0,-1256850.0,
      53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,
      11310276995381.0, -135585182899530.0,
      1307535010540395.0,     -10142299865511450.0,
      63030812099294896.0,     -311333643161390640.0,
      1206647803780373360.0,     -3599979517947607200.0,
      8037811822645051776.0,      -12870931245150988800.0,
      13803759753640704000.0,      -8752948036761600000.0,
      2432902008176640000.0]


P = Poly(flipdim(coefficients_array, 1)) #flipdim because Poly takes coefficient from the lowest order, e. g. a_0, a_1, ...
p = poly(roots_array)

P_roots = roots(P)
P_roots = flipdim(P_roots, 1)

println("Roots for P(x): ", P_roots)

println("Checking |P(z_k)|: ")
for i = 1 : size(P_roots)[1]
      println("P(",i,") = ", abs(polyval(P, P_roots[i])))
end

println("Checking |p(z_k)|: ")
for i = 1 : size(P_roots)[1]
      println("p(", i, ") = ", abs(polyval(p, P_roots[i])))
end

println("Checking |z_k - k|: ")
for i = 1 : size(P_roots)[1]
      println("z_", i, " - ", i, " = ", abs(P_roots[i] - i))
end

#=
for i = 1 : size(P_roots)[1]
      println(P_roots[i], " & ", abs(polyval(P, P_roots[i])), " & ", abs(polyval(p, P_roots[i])), " & ", abs(P_roots[i] - i), " \\\\")
end
=#
P_new = Poly(flipdim(new_coefficients_array, 1))

P_new_roots = roots(P_new)
P_new_roots = flipdim(P_new_roots, 1)


println("Checking |P'(z'_k)|: ")
for i = 1 : size(P_new_roots)[1]
      println("P'(",i,") = ", abs(polyval(P_new, P_new_roots[i])))
end

println("Checking |p(z'_k)|: ")
for i = 1 : size(P_new_roots)[1]
      println("p'(", i, ") = ", abs(polyval(p, P_new_roots[i])))
end

println("Checking |z'_k - k|: ")
for i = 1 : size(P_new_roots)[1]
      println("z'_", i, " - ", i, " = ", abs(P_new_roots[i] - i))
end

#=
for i = 1 : size(P_roots)[1]
      println(P_new_roots[i], " & ", abs(polyval(P, P_new_roots[i])), " & ", abs(polyval(p, P_new_roots[i])), " & ", abs(P_new_roots[i] - i), " \\\\")
end
=#
