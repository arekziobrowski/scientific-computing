#Obliczenia naukowe lista 2, zadanie 5
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()
using PyPlot

"""
    Compute p_n, using the following recurrence formula: p_n = p_{n - 1} + r * p_{n - 1} * (1 - p_{n - 1})
    Parameters: n - integer with n-th index of p, r - constant from the formula, p_0 - bounding value, T - type
    Returns: p - p_n value, A - array with all p values from p_1 to p_n
"""
function p(n, r, p_0, T)
    A = fill(T(0.0), n)
    p = T(p_0)
    #A[1] = p
    for i = 1 : n
        p = p + r * p * (T(1.0) - p)
        A[i] = p
    end
    return p, A
end

"""
    Compute p_n, using the following recurrence formula: p_n = p_{n - 1} + r * p_{n - 1} * (1 - p_{n - 1})
    Values are truncated to 3 decimal places after the 10th iteration
    Parameters: n - integer with n-th index of p, r - constant from the formula, p_0 - bounding value, T - type
    Returns: p - p_n value, A - array with all p values from p_1 to p_n
"""
function p_truncate(n, r, p_0, T)
    A = fill(T(0.0), n)
    p = T(p_0)
    #A[1] = p
    for i = 1 : n
        p = p + r * p * (T(1.0) - p)
        if i == 10
            p = trunc(p, 3)
        end
        A[i] = p
    end
    return p, A
end

#println(p(40, 3, 0.01, Float32)[1])
#println(p_truncate(40, 3, 0.01, Float32)[1])

A = p(40, 3, 0.01, Float32)[2]
B = p_truncate(40, 3, 0.01, Float32)[2]

for i  = 1 : size(A)[1]
    println("n = ", i, ", p_", i, " = ", A[i], ", p_", i, "trunc = ", B[i])
end

plot(p(40, 3, 0.01, Float32)[2])
plot(p(40, 3, 0.01, Float64)[2])
