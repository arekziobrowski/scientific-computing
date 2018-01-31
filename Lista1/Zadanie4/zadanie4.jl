#Obliczenia naukowe lista 1, zadanie 4
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

"""
    Find a number in range [a, b] that doesn't meet the x * (1.0 / x) == 1 condition
    Parameters: a - left parameter, b - right parameter in [a, b] range
    Returns: x - floating point number anomaly
"""
function find_anomaly(a::Float64, b::Float64)
    x = a
    while x * Float64(1.0 / x) == 1 && x < Float64(b)
        x = nextfloat(x)
    end
    return x
end

println("Anomaly in (1, 2):\t", find_anomaly(nextfloat(1.0), prevfloat(2.0)))
println("Smallest anomaly:\t", find_anomaly(nextfloat(typemin(Float64)), prevfloat(typemax(Float64))))
