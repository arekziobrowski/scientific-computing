#Obliczenia naukowe lista 3, zadanie 1-3 TEST
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728


workspace()

include("../zadanie1-3.jl")

using EquationSolving

δ = Float64(0.5) * 10.0e-5
ϵ = Float64(0.5) * 10.0e-5


println("x^2 - 1")
println("Bisection method:  ",mbisekcji(x -> x^2 - 1, -0.5, 0.5, δ, ϵ)) #expected error (no sign change)
println("Newton's method: ",mstycznych(x -> x^2 - 1, x -> 2*x, 0.0, δ, ϵ, 50)) #expected error (derivative value close to 0)
println("Secant method: ",msiecznych(x -> x^2 - 1, -0.5, 0.5, δ, ϵ, 50)) #expected error (value not found in the given number of iterations)

println("x^2 - 9")
println("Bisection method:  ",mbisekcji(x -> x^2 - 9, 1.5, 6.5, δ, ϵ))
println("Newton's method: ",mstycznych(x -> x^2 - 9, x -> 2*x, 1.5, δ, ϵ, 50))
println("Secant method: ",msiecznych(x -> x^2 - 9, 1.0, 2.0, δ, ϵ, 50))
