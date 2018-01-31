#Obliczenia naukowe lista 4, zadanie 6
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()
include("../Zadanie1-4/zadanie1-4.jl")
using InterpolationUtils

#comment out the plots that shall not be on the graph

#rysujNnfx(x -> abs(x), -1.0, 1.0, 5)
#rysujNnfx(x -> abs(x), -1.0, 1.0, 10)
#rysujNnfx(x -> abs(x), -1.0, 1.0, 15)


#rysujNnfx(x -> 1 / (1 + x^2), -5.0, 5.0, 5);
#rysujNnfx(x -> 1 / (1 + x^2), -5.0, 5.0, 10);
rysujNnfx(x -> 1 / (1 + x^2), -5.0, 5.0, 15);
