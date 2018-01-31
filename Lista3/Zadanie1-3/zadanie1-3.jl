#Obliczenia naukowe lista 3, zadanie 1-3
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728


module EquationSolving

export mbisekcji, mstycznych, msiecznych
export INCORRECT_MAXIT, INCORRECT_DELTA_EPSILON, BISECTION_NO_ERR, BISECTION_NO_SGN_CHANGE, NEWTON_CONVERGENCE, NEWTON_NOT_PRECISE, NEWTON_DERIVATIVE_CLOSE_ZERO,
        SECANT_CONVERGENCE, SECANT_NOT_PRECISE

#Error flag returned when δ or ϵ from input is less or equal to 0
global INCORRECT_DELTA_EPSILON = 3
#Error flag returned when maxit from input is less or equal to 0
global INCORRECT_MAXIT = 4

#Error flag returned when bisection algorithm returns a correct answer
global BISECTION_NO_ERR = 0
#Error flag returned when bisection algorithm meets an input with no change of sign for f(a) and f(b), where a and b are ends of [a,b]
global BISECTION_NO_SGN_CHANGE = 1

#Error flag returned when the Newton algorithm returns a correct answer
global NEWTON_CONVERGENCE = 0
#Error flag returned when the Newton algorithm's result doesn't meet the set precision in maxit iterations
global NEWTON_NOT_PRECISE = 1
#Error flag returned when a derivative has a value close to 0 in Newton's algorithm
global NEWTON_DERIVATIVE_CLOSE_ZERO = 2

#Error flag returned when the Secant algorithm returns a correct answer
global SECANT_CONVERGENCE = 0
#Error flag returned when the Secant algorithm's result doesn't meet the set precision in maxit iterations
global SECANT_NOT_PRECISE = 1

"""
    Function that finds the roots of the function f using the bisection method
    Input:
        f - anonymous function
        a - left side of the input interval [a, b]
        b - right side of the input interval [a, b]
        delta - precision of the zero of the function x
        epsilon - precision of the f(x)
    Output:
        (r, v, it, err) - tuple of 4
        r - zero of the function in the given interval
        v - value of f(r)
        it - number of iterations that have been made
        err - error flag
"""
function mbisekcji(f, a::Float64, b::Float64, delta::Float64, epsilon::Float64)

    if a > b
        temp = a
        a = b
        b = temp
    end

    if delta <= 0.0 || epsilon <= 0.0
        return "error", "error", "error", INCORRECT_DELTA_EPSILON
    end

    u = f(a)
    v = f(b)
    e = b - a
    if sign(u) == sign(v)
        return "error", "error", "error", BISECTION_NO_SGN_CHANGE
    end
    e = e / 2
    c = a + e
    w = f(c)
    it = 1
    while abs(e) > delta && abs(w) > epsilon
        it = it + 1
        if sign(w) != sign(u)
            b = c
            v = w
        else
            a = c
            u = w
        end
        e = e / 2
        c = a + e
        w = f(c)
    end
    return c, w, it, BISECTION_NO_ERR
end

"""
    Function that finds the roots of the function f using the Newton method
    Input:
        f - anonymous function
        pf - anonymous function, being the derivative of f
        x0 - approximation of the zero of the function
        delta - precision of the zero of the function x
        epsilon - precision of the f(x)
        maxit - maximum number of iterations that can be made
    Output:
        (r, v, it, err) - tuple of 4
        r - zero of the function
        v - value of f(r)
        it - number of iterations that have been made
        err - error flag
"""
function mstycznych(f, pf, x0::Float64, delta::Float64, epsilon::Float64, maxit::Int)

    local x1 = 0.0

    if delta <= 0.0 || epsilon <= 0.0
        return "error", "error", "error", INCORRECT_DELTA_EPSILON
    end

    if maxit <= 0
        return "error", "error", "error", INCORRECT_MAXIT
    end

    v = f(x0)
    if abs(v) < epsilon
        return x0, v, 1, NEWTON_CONVERGENCE
    end

    for k = 1 : maxit
        if abs(pf(x0)) < eps(Float64)
            return "error", "error", k, NEWTON_DERIVATIVE_CLOSE_ZERO
        end
        x1 = x0 - v / pf(x0)
        v = f(x1)
        if abs(x1 - x0) < delta || abs(v) < epsilon
            return x1, v, k, NEWTON_CONVERGENCE
        end
        x0 = x1
    end
    return x1, v, maxit + 1, NEWTON_NOT_PRECISE
end


"""
    Function that finds the roots of the function f using the secant method
    Input:
        f - anonymous function
        x0 - approximation of the zero of the function f
        x1 - approximation of the zero of the function f
        delta - precision of the zero of the function x
        epsilon - precision of the f(x)
        maxit - maximum number of iterations that can be made
    Output:
        (r, v, it, err) - tuple of 4
        r - zero of the function
        v - value of f(r)
        it - number of iterations that have been made
        err - error flag
"""
function msiecznych(f, x0::Float64, x1::Float64, delta::Float64, epsilon::Float64, maxit::Int)

    if delta <= 0.0 || epsilon <= 0.0
        return "error", "error", "error", INCORRECT_DELTA_EPSILON
    end

    if maxit <= 0
        return "error", "error", "error", INCORRECT_MAXIT
    end

    fx0 = f(x0)
    fx1 = f(x1)
    for k = 1 : maxit
        if abs(fx0) > abs(fx1)
            temp = x0
            x0 = x1
            x1 = temp

            temp = fx0
            fx0 = fx1
            fx1 = temp
        end
        s = (x1 - x0) / (fx1 - fx0)
        x1 = x0
        fx1 = fx0
        x0 = x0 - fx0 * s
        fx0 = f(x0)
        if abs(fx0) < epsilon || abs(x1 - x0) < delta
            return x0, fx0, k, SECANT_CONVERGENCE
        end
    end
    return x0, fx0, maxit + 1, SECANT_NOT_PRECISE
end

end
