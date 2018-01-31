#Obliczenia naukowe lista 3, zadanie 6
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728


workspace()
include("../Zadanie1-3/zadanie1-3.jl")

using EquationSolving

δ = 10.0e-5
ϵ = 10.0e-5

f1 = x -> e^(1 - x) - 1
df1 = x -> -e^(1 - x)

f2 = x -> x * e^(-x)
df2 = x -> -e^(-x) * (x - 1)

println("e^(1 - x) - 1")
println("Bisection method:  ",mbisekcji(f1, 0.0, 1.5, δ, ϵ))
println("Newton's method: ",mstycznych(f1, df1, 7.5, δ, ϵ, 1140))
println("Secant method: ",msiecznych(f1, 0.0, 0.5, δ, ϵ, 40))

println("x * e^(-x)")
println("Bisection method:  ",mbisekcji(f2, -0.3, 0.4, δ, ϵ))
println("Newton's method: ",mstycznych(f2, df2, 1.0, δ, ϵ, 40))
println("Secant method: ",msiecznych(f2, -20.0, -19.0, δ, ϵ, 40))
