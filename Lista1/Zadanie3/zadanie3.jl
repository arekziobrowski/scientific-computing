#Obliczenia naukowe lista 1, zadanie 3
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

"""
    Function to check density of numbers in the given range [a, b]
    Parameters: a - left parameter, b - right parameter in [a, b] range
    Returns:
"""
function check_density(a::Float64, b::Float64)
    x = a
    δ = Float64(2.0) ^ -52
    k = 0;
    while x < b
        k = k + 1
        x = a + Float64(k * δ)
    end
    if k == 2 ^ 52
        println("Numbers in range from ", a, " to ", b, " are evenly spread.")
    else
        println("Numbers are NOT evenly spread.")
    end
end

"""
    Check density of numbers in the given range [a, b] by showing the first 10 steps
    It is faster than check_density(::Float64, ::Float64)
    Parameters: a - left parameter, b - right parameter in [a, b] range
    Returns:
"""
function check_density_bits(a::Float64, b::Float64)
    x = a
    δ = 2.0 ^ -52
    k = 0
    for i = 1:10
        println(a, " + ", k, " * δ:\t",  bits(x))
        k = k + 1
        x = a + Float64(k * δ)

    end
    println(b, " - δ:\t", bits(prevfloat(b)))
end

check_density_bits(1.0, 2.0)
