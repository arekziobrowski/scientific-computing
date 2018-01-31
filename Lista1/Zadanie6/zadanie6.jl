#Obliczenia naukowe lista 1, zadanie 6
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

"""
    Compute sqrt(x^2 + 1) - 1
    Parameters: x - floating point number
    Returns: sqrt(x^2 + 1) - 1
"""
function f(x::Float64)
    return Float64(sqrt(x^2 + 1) - 1)
end

"""
    Compute x^2 / (sqrt(x^2 + 1) + 1)
    Parameters: x - floating point number
    Returns: x^2 / (sqrt(x^2 + 1) + 1)
"""
function g(x::Float64)
    return Float64(x^2/(sqrt(x^2 + 1) + 1))
end

for i = 1 : 15
    x = Float64(8)^-i
    println("f(8^", i, "):\t", f(x), ", g(8^", i, "):\t", g(x))
end
