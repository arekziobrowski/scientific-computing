#Obliczenia naukowe lista 3, zadanie 5
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()
include("../Zadanie1-3/zadanie1-3.jl")

using EquationSolving

δ = 10.0e-4
ϵ = 10.0e-4

f = x -> 3*x - e^x
df = x -> 3 - e^x


println("3*x - e^x")
a = 0.0
b = 1.0
println("Bisection method for a =  ", a, ", b = ", b, " = ",mbisekcji(f, a, b, δ, ϵ))
a = 1.0
b = 2.0
println("Bisection method for a =  ", a, ", b = ", b, " = ",mbisekcji(f, a, b, δ, ϵ))
