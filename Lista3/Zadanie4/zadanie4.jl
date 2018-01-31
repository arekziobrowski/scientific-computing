#Obliczenia naukowe lista 3, zadanie 4
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728


workspace()
include("../Zadanie1-3/zadanie1-3.jl")

using EquationSolving

δ = 5.0e-6
ϵ = 5.0e-6

f = x -> sin(x) - (0.5 * x)^2
df = x -> cos(x) - 0.5 * x

println("sin(x) - (1/2 * x)^2")
tic()
println("Bisection method:  ",mbisekcji(f, -1.0, 2.0, δ, ϵ))
toc()
tic()
println("Newton's method: ",mstycznych(f, df, 1.5, δ, ϵ, 50))
toc()
tic()
println("Secant method: ",msiecznych(f, 1.0, 2.0, δ, ϵ, 50))
toc()
