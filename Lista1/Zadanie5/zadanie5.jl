#Obliczenia naukowe lista 1, zadanie 5
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728


"""
    Cast vectors to the given type
    Parameters: T - type
    Returns: x - vector x in the given arithmetic, y - vector y in the given arithmetic
"""
function cast_vectors(T)
    x = [T(2.718281828), T(-3.141592654), T(1.414213562), T(0.5772156649), T(0.3010299957)]
    y = [T(1486.2497), T(878366.9879), T(-22.37492), T(4773714.647), T(0.000185049)]
    return x, y
end

"""
    Split an array to arrays of its nonnegative and negative values
    Parameters: x - array
    Returns: nonnegative - array with nonnegative values of x, negative - array with negative values of x
"""
function split_sign(x)
    nonnegative = []
    negative = []
    for i = 1 : length(x)
        if x[i] < 0
            push!(negative, x[i])
        else
            push!(nonnegative, x[i])
        end
    end
    return nonnegative, negative
end

"""
    Compute approximation error
    Parameters: x - rounded value, y - exact value
    Returns: approximation error
"""
function approximation_error(x, y)
    return Float64(abs(x - y)) / Float64(abs(y))
end

"""
    Algorithm (a)
    Parameters: x - vector, y - vector, T - type
    Returns: S - scalar product
"""
function algorithm1(x, y, T)
    S = T(0.0)
    for i = 1 : length(x)
        S = S + T(x[i]) * T(y[i])
    end
    return S
end

"""
    Algorithm (b)
    Parameters: x - vector, y - vector, T - type
    Returns: S - scalar product
"""
function algorithm2(x, y, T)
    S = T(0.0)
    for i = length(x) : - 1 : 1
        S = S + T(x[i]) * T(y[i])
    end
    return S
end

"""
    Algorithm (c)
    Parameters: x - vector, y - vector, T - type
    Returns: S - scalar product
"""
function algorithm3(x, y, T)
    xy = []
    for i = 1 : length(x)
        push!(xy, T(x[i]) * T(y[i]))
    end
    xy_nonnegative, xy_negative = split_sign(xy)
    sort(xy_nonnegative, rev = true) #sorting nonnegative numbers in reversed order (from largest to smallest)
    sort(xy_negative) #sorting negative numbers (from smallest to largest)

    S_plus = T(0.0)
    for i = 1 : length(xy_nonnegative)
        S_plus = S_plus + T(xy_nonnegative[i])
    end

    S_minus = T(0.0)
    for i = 1 : length(xy_negative)
        S_minus = S_minus + T(xy_negative[i])
    end

    S = S_plus + S_minus
    return S
end

"""
    Algorithm (d)
    Parameters: x - vector, y - vector, T - type
    Returns: S - scalar product
"""
function algorithm4(x, y, T)
    xy = []
    for i = 1 : length(x)
        push!(xy, T(x[i]) * T(y[i]))
    end
    xy_nonnegative, xy_negative = split_sign(xy)
    sort(xy_nonnegative) #sorting nonnegative numbers (from smallest to largest)
    sort(xy_negative, rev = true) #sorting negative numbers in reversed order (from largest to smallest)

    S_plus = T(0.0)
    for i = 1 : length(xy_nonnegative)
        S_plus = S_plus + T(xy_nonnegative[i])
    end

    S_minus = T(0.0)
    for i = 1 : length(xy_negative)
        S_minus = S_minus + T(xy_negative[i])
    end

    S = S_plus + S_minus
    return S
end

value = Float64(-1.00657107000000e-11) #exact value of scalar product

x, y = cast_vectors(Float32)
algorithm1_value = algorithm1(x, y, Float32)
algorithm2_value = algorithm2(x, y, Float32)
algorithm3_value = algorithm3(x, y, Float32)
algorithm4_value = algorithm4(x, y, Float32)
println("algorithm1(x, y, Float32):\t", algorithm1_value, ",\t\t\tapproximation error: \t", approximation_error(algorithm1_value, value))
println("algorithm2(x, y, Float32):\t", algorithm2_value, ",\t\t\tapproximation error: \t", approximation_error(algorithm2_value, value))
println("algorithm3(x, y, Float32):\t", algorithm3_value, ",\t\t\t\tapproximation error: \t", approximation_error(algorithm3_value, value))
println("algorithm4(x, y, Float32):\t", algorithm4_value, ",\t\t\t\tapproximation error: \t", approximation_error(algorithm4_value, value))

x, y = cast_vectors(Float64)
algorithm1_value = algorithm1(x, y, Float64)
algorithm2_value = algorithm2(x, y, Float64)
algorithm3_value = algorithm3(x, y, Float64)
algorithm4_value = algorithm4(x, y, Float64)
println("algorithm1(x, y, Float64):\t", algorithm1_value, ",\t\tapproximation error: \t", approximation_error(algorithm1_value, value))
println("algorithm2(x, y, Float64):\t", algorithm2_value, ",\tapproximation error: \t", approximation_error(algorithm2_value, value))
println("algorithm3(x, y, Float64):\t", algorithm3_value, ",\t\t\t\tapproximation error: \t", approximation_error(algorithm3_value, value))
println("algorithm4(x, y, Float64):\t", algorithm4_value, ",\t\t\t\tapproximation error: \t", approximation_error(algorithm4_value, value))
