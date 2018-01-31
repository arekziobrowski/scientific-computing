#Obliczenia naukowe lista 4, zadanie 5
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()
include("../Zadanie1-4/zadanie1-4.jl")
using InterpolationUtils

#comment out the plots that shall not be on the graph

rysujNnfx(x -> e^x, 0.0, 1.0, 5)
rysujNnfx(x -> e^x, 0.0, 1.0, 10)
rysujNnfx(x -> e^x, 0.0, 1.0, 15)


rysujNnfx(x -> (x^2) * sin(x), -1.0, 1.0, 5);
rysujNnfx(x -> (x^2) * sin(x), -1.0, 1.0, 10);
rysujNnfx(x -> (x^2) * sin(x), -1.0, 1.0, 15);
