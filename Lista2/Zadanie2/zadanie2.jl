#Obliczenia naukowe lista 2, zadanie 2
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()

using Gadfly
f(x) = e^x * log(1 + e^-x)
draw(PDF("gadfly.pdf", 654mm, 350mm), plot(f, -15, 45))
