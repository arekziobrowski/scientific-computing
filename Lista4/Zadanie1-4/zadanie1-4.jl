#Obliczenia naukowe lista 4, zadanie 1-4
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728


module InterpolationUtils

export ilorazyRoznicowe, warNewton, naturalna, rysujNnfx

using PyPlot


"""
    Function that finds the divided differences
    Input:
        x - Vector of size n + 1 containing the nodes
        f - Vector of size n + 1 containing the values of interpolated function in nodes
    Output:
        fx - Vector of size n + 1 containing the divided differences
"""
function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64})
    n = length(x)
    if(n != length(f))
        error("Length of x::Vector{Float64} isn't equal to length of f::Vector{Float64}.")
        return
    end

    fx = Vector{Float64}(n)
    for i = 1 : n
        fx[i] = f[i]
    end
    for j = 1 : n
        for i = n : -1 : j + 1
            fx[i] = (fx[i] - fx[i - 1]) / (x[i] - x[i - j])
        end
    end
    return fx;
end

"""
    Function that finds the value of interpolating polynomial in point
    Input:
        x - Vector of size n + 1 containing the nodes
        fx - Vector of size n + 1 containing the divided differences
        t - point in which we want to compute the value of the interpolating polynomial
    Output:
        nt - value of the interpolating polynomial in the point t
"""
function warNewton(x::Vector{Float64}, fx::Vector{Float64}, t::Float64)
    w = fx[length(x)]
    for i = length(x) - 1 : -1 : 1
        w = fx[i] + w * (t - x[i])
    end
    return w
end

"""
    Function that finds the coefficients of the power form of the Newton form interpolating polynomial
    Input:
        x - Vector of size n containing the nodes
        fx - Vector of size n + 1 containing the divided differences
    Output:
        a - Vector of size n + 1 containing coefficients of the power form
"""
function naturalna(x::Vector{Float64}, fx::Vector{Float64})
    n = length(x)
    a = fx
    z = 0
    for i = 1 : n
        for j = n : -1 : i
            a[j] = a[j] + (z - x[j]) * a[j + 1]
        end
        for i = n : -1 : 2
            x[i] = x[i - 1]
        end
        x[1] = z
    end
    return a
end

"""
    Function that interpolates a function and plots the interpolating polynomial alongside with the function
    Input:
        f - anonymous function that we want to interpolate
        a - left end of the interpolating range
        b - right end of the interpolating range
        n - degree of the interpolating polynomial
    Output:
        graphs the function f and interpolating polynomial in the given range [a, b]
"""
function rysujNnfx(f, a::Float64, b::Float64, n::Int)
    h = (b - a) / n
    x = Vector{Float64}(n + 1)
    y = Vector{Float64}(n + 1)

    for i = 1 : n + 1
        x[i] = a + (i - 1) * h
        y[i] = f(x[i])
    end

    fx = ilorazyRoznicowe(x, y)

    pointsNumber = n*100
    y = Vector{Float64}(pointsNumber)
    interpolation = Vector{Float64}(pointsNumber)
    error = Vector{Float64}(pointsNumber)

    h1 = (b - a) / pointsNumber
    for i = 1 : pointsNumber
        t = a + (i - 1) * h1
        y[i] = f(t)
        interpolation[i] = warNewton(x, fx, t)
        error[i] = y[i] - interpolation[i]
    end

    #plot(linspace(a, b, pointsNumber), error)
    plot(linspace(a, b, pointsNumber), y, color="blue")
    plot(linspace(a, b, pointsNumber), interpolation, color="red")
    #plot(linspace(a, b, pointsNumber), y, label=string(n))
    #plot(linspace(a, b, pointsNumber), interpolation, label=string(n, " exact"))
    #legend(bbox_to_anchor=[1,1],loc=2,borderaxespad=0)
end



end
