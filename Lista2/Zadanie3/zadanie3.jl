#Obliczenia naukowe lista 2, zadanie 3
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()
include("hilb.jl")
include("matcond.jl")

"""
    Pretty print matrix
    Parameters: A - matrix
    Returns: nothing
"""
function print_matrix(A)
    for i = 1 : size(A)[1]
        print("| ")
        for j = 1 : size(A)[2]
            @printf("%.20f ", A[i, j])
        end
        println("|")
    end
end

"""
    Pretty print vector
    Parameters: b - vector
    Returns: nothing
"""
function print_vector(b)
    for i = 1 : size(b)[1]
        println("| ", b[i], " |")
    end
end

"""
    Solve Ax = b for Hilbert's matrix of size n
    Parameters: n - matrix size
    Returns: nothing
"""
function solve_hilbert(n)
    A = hilb(n)
    x = ones(n)
    b = zeros(n)
    b = A*x

    gauss = Array{Float64, 1}(A \ b)
    inverse = Array{Float64, 1}(inv(A) * b)

    gaussError = norm(x - gauss) / norm(x)
    inverseError = norm(x - inverse) / norm(x)
    println("n = ", n, ", rank(A) = ", rank(A), ", cond(A) = ", cond(A), ", gaussError = ", gaussError, ", inverseError = ", inverseError)
    #println(n, " & ", rank(A), " & ", cond(A), " & ", gaussError, " & ", inverseError, " \\\\")
end

"""
    Solve Ax = b for random matrix of size n and condition number c
    Parameters: n - matrix size, c - condition number
    Returns: nothing
"""
function solve_random(n, c)
    A = matcond(n, c)
    while det(A) == 0
        A = matcond(n, c)
    end
    x = ones(n)
    b = zeros(n)
    b = A * x

    gauss = A \ b
    inverse = inv(A) * b

    gaussError = norm(x - gauss) / norm(x)
    inverseError = norm(x - inverse) / norm(x)
    println("n = ", n, ", rank(A) = ", rank(A), ", cond(A) = ", c, ", gaussError = ", gaussError, ", inverseError = ", inverseError)
    #println(n, " & ", rank(A), " & ", c, " & ", gaussError, " & ", inverseError, " \\\\")
end

println("GENERATING MATRIX WITH HILB(n)")
for i = 1 : 20
    solve_hilbert(i)
end

println("GENERATING MATRIX WITH MATCOND(n, c)")
n = [5, 10, 20]
c = [1, 10, 10^3, 10^7, 10^12, 10^16]
for i = 1 : size(n)[1]
    for j = 1 : size(c)[1]
        solve_random(n[i], Float64(c[j]))
    end
end
