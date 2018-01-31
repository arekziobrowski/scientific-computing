#Obliczenia naukowe lista 2, zadanie 6
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

using PyPlot

"""
    Compute x_n, using the following recurrence formula: x_n = x_{n - 1}^2 + c
    Parameters: n - integer with n-th index of p, c - constant from the formula, x_0 - bounding value
    Returns: x - x_n value, A - array with all x values from x_0 to x_n
"""
function x(n, c, x_0)
    A = fill(0.0, n + 1)
    x = x_0
    A[1] = x
    for i = 1 : n
        x = x^2 + c
        A[i + 1] = x
    end
    return x, A
end

x_0_args = [1.0, 2.0, 1.99999999999999, 1.0, -1.0, 0.75, 0.25]
c_args = [-2.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0]

for i = 1 : size(x_0_args)[1]
    println("x_0 = ", x_0_args[i], ", c = ", c_args[i], ", x(40, ", c_args[i], ", ", x_0_args[i], ") = ", x(40, c_args[i], x_0_args[i])[1])
end


plot(x(40, c_args[7], x_0_args[7])[2])
