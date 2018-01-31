#Obliczenia naukowe lista 1, zadanie 7
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

"""
    Compute sin(x) + cos(3x)
    Parameters: x - floating point number
    Returns: sin(x) + cos(3x)
"""
function f(x::Float64)
    return Float64(sin(x)) + Float64(cos(3.0*x))
end

"""
    Compute derivative(sin(x) + cos(3x))
    Parameters: x - floating point number
    Returns: derivative(sin(x) + cos(3x))
"""
function derivative_f(x::Float64)
    return Float64(cos(x)) - 3.0 * Float64(sin(3.0*x))
end

"""
    Compute approximate derivative(sin(x) + cos(3x))
    Parameters: x - floating point number, h - floating point number
    Returns: approximate derivative(sin(x) + cos(3x))
"""
function approximate_derivative_f(x::Float64, h::Float64)
    return Float64((f(x + h) - f(x)) / h)
end

"""
    Compute absolute error |x - rd(x)|
    Parameters: x - floating point number, y - floating point number
    Returns: |x - rd(x)|
"""
function absolute_error(x::Float64, y::Float64)
    return abs(x - y)
end


x0 = 1.0

for i = 0 : 54
    h = Float64(2.0) ^ -i
    derivative = derivative_f(x0)
    approximate = approximate_derivative_f(x0, h)
    println("i = ", i, ", derivative: ", derivative, ", approximate: ", approximate, ", absolute error: ", absolute_error(derivative, approximate), ", 1 + h = ", 1 + h)
end
