#Obliczenia naukowe lista 5, zadanie 1-3
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728
module Blocksys

export load_matrix, load_vector, compute_vector, gaussian_elimination!, gaussian_elimination_pivoting!, lufactorization!, lufactorization_pivoting!, permuted, solve_linear_system,
        solve_linear_system!, save_vector

#variable storing the boolean value indicating whether the vector b was computed or read from a file
bcomputed = true

#=
    Function that loads a sparse matrix A from a file
    Input:
        filename::String - name of the file containing sparse matrix
    Output:
        A::SparseMatrixCSC - sparse matrix
=#
function load_matrix(filename::String)
    open(filename) do f
        line = readline(f)
        (n, l) = split(line)
        (n, l) = (parse(Int64, n), parse(Int64, l))
        Is = Array{Int64, 1}()
        Js = Array{Int64, 1}()
        Vs = Array{Float64, 1}()
        while !eof(f)
            line = split(readline(f))
            push!(Is, parse(Int64, line[1]))
            push!(Js, parse(Int64, line[2]))
            push!(Vs, parse(Float64, line[3]))
        end
        return (sparse(Is, Js, Vs), n, l)
    end
end

#=
    Function that loads a vector b from a file
    Input:
        filename::String - name of the file containing vector
    Output:
        b::Vector{Float64} - vector b
=#
function load_vector(filename::String)
    open(filename) do f
        line = readline(f)
        n = parse(Int64, line)
        b = Vector{Float64}()
        while !eof(f)
            line = readline(f)
            push!(b, parse(Float64, line))
        end
        global bcomputed = false
        return b
    end
end

#=
    Function that computes a vector b assuming that x is (1, 1, ..., 1)^T
    Input:
        A - sparse matrix
        n - dimensions of the sparse matrix
    Output:
        b::Vector{Float64} - vector b
=#
function compute_vector(A::SparseMatrixCSC{Float64, Int64}, n::Int64)
    b = zeros(Float64, n)
    vals = nonzeros(A)
    rows = rowvals(A)
    for i = 1:n
        for j in nzrange(A, i)
            b[rows[j]] += vals[j]
        end
    end
    global bcomputed = true
    return b
end

#=
    Function that solves the linear system Ax = b
    Input:
        A - sparse matrix
        b - vector
        n - dimensions of the sparse matrix
        l - dimensions of the block of the sparse matrix
    Output:
        x::Vector{Float64} - vector containing solution of the linear system Ax = b
=#
function gaussian_elimination!(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    for k = 1:n-2
        if k % l <= l - 2
            rows_to_set = l - k % l
        else
            rows_to_set = 2*l - k % l
        end

        for i = k+1:k+rows_to_set
            if abs(A[k, k]) <= eps(Float64)
                error("Zero as the pivot")
            end
            multiplier = A[i, k] / A[k, k]
            A[i, k] = 0
            if k + l <= n
                columns_to_set = l
            else
                columns_to_set = l - k % l
            end
            for j = k+1:k+columns_to_set
                A[i, j] = A[i, j] - multiplier * A[k, j]
            end
            b[i] = b[i] - multiplier * b[k]
        end
    end
    multiplier = A[n, n - 1] / A[n - 1, n - 1]
    A[n, n - 1] = 0
    A[n, n] = A[n, n] - multiplier * A[n - 1, n]
    b[n] = b[n] - multiplier * b[n - 1]

    b[n] = b[n] / A[n, n]
    for i in n-1:-1:1
		s = b[i]
        offset = min(n, i + l)
		for j = offset:-1:i+1
			s = s - A[i, j] * b[j]
		end
		b[i] = s / A[i, i]
	end

    return b
end

#=
    Function that solves the linear system Ax = b using partial pivoting
    Input:
        A - sparse matrix
        b - vector
        n - dimensions of the sparse matrix
        l - dimensions of the block of the sparse matrix
    Output:
        x::Vector{Float64} - vector containing solution of the linear system Ax = b
=#
function gaussian_elimination_pivoting!(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    s = zeros(n)
    p = Vector{Int64}(n)

    for i = 1:n
        p[i] = i
        vals = nonzeros(A)
        rows = rowvals(A)
        for j in nzrange(A, i)
            val = abs(vals[j])
            if val > s[rows[j]]
                s[rows[j]] = val
            end
        end
    end
    for k = 1:n-2
        if k % l <= l - 2
            rows_to_set = l - k % l
        else
            rows_to_set = 2*l - k % l
        end
        mainrow = k
        max = abs(A[p[k], k]) / s[p[k]]
        for iter = k+1:k+rows_to_set
            val = abs(A[p[iter], k]) / s[p[iter]]
            if val > max
                mainrow = iter
                max = val
            end
        end

        p[k], p[mainrow] = p[mainrow], p[k]
        s[k], s[mainrow] = s[mainrow], s[k]

        for i = k+1:k+rows_to_set
            multiplier = A[p[i], k] / A[p[k], k]
            A[p[i], k] = 0
            if k + 2*l <= n
                columns_to_set = 2*l
            else
                columns_to_set = -k + n
            end
            for j = k+1:k+columns_to_set

                A[p[i], j] = A[p[i], j] - multiplier * A[p[k], j]
            end
            b[p[i]] = b[p[i]] - multiplier * b[p[k]]
        end
    end
    mainrow = n - 1
    max = abs(A[p[n - 1], n - 1]) / s[p[n - 1]]
    for iter = n - 1:n
        val = abs(A[p[iter], n - 1]) / s[p[iter]]
        if val > max
            mainrow = iter
            max = val
        end
    end
    p[n - 1], p[mainrow] = p[mainrow], p[n - 1]
    s[n - 1], s[mainrow] = s[mainrow], s[n - 1]
    multiplier = A[p[n], n - 1] / A[p[n - 1], n - 1]
    A[p[n], n - 1] = 0
    A[p[n], n] = A[p[n], n] - multiplier * A[p[n - 1], n]
    b[p[n]] = b[p[n]] - multiplier * b[p[n - 1]]

    b[p[n]] = b[p[n]] / A[p[n], n]
    for i in n-1:-1:1
		sum = b[p[i]]
        offset = min(n, p[i] + 2*l + 1)
		for j = offset:-1:i+1
			sum = sum - A[p[i], j] * b[p[j]]
		end
		b[p[i]] = sum / A[p[i], i]
	end
    return b
end

#=
    Function that computes an LU factorization of the matrix A
    Input:
        A - sparse matrix
        n - dimensions of the sparse matrix
        l - dimensions of the block of the sparse matrix
    Output:
        (L, U):
            L::SparseMatrixCSC - lower triangular matrix
            U::SparseMatrixCSC - upper triangular matrix
=#
function lufactorization!(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    Is = Vector{Int64}()
    Js = Vector{Int64}()
    Vs = Vector{Float64}()
    for k = 1:n-2
        if k % l <= l - 2
            rows_to_set = l - k % l
        else
            rows_to_set = 2*l - k % l
        end
        push!(Is, k)
        push!(Js, k)
        push!(Vs, 1.0)
        for i = k+1:k+rows_to_set
            if abs(A[k, k]) <= eps(Float64)
                error("Zero as the pivot")
            end
            multiplier = A[i, k] / A[k, k]
            A[i, k] = 0
            push!(Is, i)
            push!(Js, k)
            push!(Vs, multiplier)
            if k + l <= n
                columns_to_set = l
            else
                columns_to_set = l - k % l
            end
            for j = k+1:k+columns_to_set
                A[i, j] = A[i, j] - multiplier * A[k, j]
            end
        end
    end
    push!(Is, n - 1)
    push!(Js, n - 1)
    push!(Vs, 1.0)
    multiplier = A[n, n - 1] / A[n - 1, n - 1]
    A[n, n - 1] = 0
    push!(Is, n)
    push!(Js, n - 1)
    push!(Vs, multiplier)
    A[n, n] = A[n, n] - multiplier * A[n - 1, n]
    push!(Is, n)
    push!(Js, n)
    push!(Vs, 1.0)
    dropzeros!(A) #to improve the performance at the cost of memory, remove this line
    return (sparse(Is, Js, Vs), A)
end

#=
    Function that computes an LU factorization of the matrix A using partial pivoting
    Input:
        A - sparse matrix
        n - dimensions of the sparse matrix
        l - dimensions of the block of the sparse matrix
    Output:
        (L, U, p):
            L::SparseMatrixCSC - lower triangular matrix
            U::SparseMatrixCSC - upper triangular matrix
            p::Vector{Int64} - permutation
=#
function lufactorization_pivoting!(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
    s = zeros(n)
    p = Vector{Int64}(n)

    Is = Vector{Int64}()
    Js = Vector{Int64}()
    Vs = Vector{Float64}()

    for i = 1:n
        p[i] = i
        vals = nonzeros(A)
        rows = rowvals(A)
        for j in nzrange(A, i)
            val = abs(vals[j])
            if val > s[rows[j]]
                s[rows[j]] = val
            end
        end
    end

    for k = 1:n-2
        if k % l <= l - 2
            rows_to_set = l - k % l
        else
            rows_to_set = 2*l - k % l
        end

        mainrow = k
        max = abs(A[p[k], k]) / s[p[k]]
        for iter = k+1:k+rows_to_set
            val = abs(A[p[iter], k]) / s[p[iter]]
            if val > max
                mainrow = iter
                max = val
            end
        end

        p[k], p[mainrow] = p[mainrow], p[k]
        s[k], s[mainrow] = s[mainrow], s[k]

        push!(Is, p[k])
        push!(Js, k)
        push!(Vs, 1.0)

        for i = k+1:k+rows_to_set
            multiplier = A[p[i], k] / A[p[k], k]
            A[p[i], k] = 0
            push!(Is, p[i])
            push!(Js, k)
            push!(Vs, multiplier)
            if k + 2*l <= n
                columns_to_set = 2*l
            else
                columns_to_set = -k + n
            end
            for j = k+1:k+columns_to_set

                A[p[i], j] = A[p[i], j] - multiplier * A[p[k], j]
            end
        end
    end
    mainrow = n - 1
    max = abs(A[p[n - 1], n - 1]) / s[p[n - 1]]
    for iter = n - 1:n
        val = abs(A[p[iter], n - 1]) / s[p[iter]]
        if val > max
            mainrow = iter
            max = val
        end
    end
    p[n - 1], p[mainrow] = p[mainrow], p[n - 1]
    s[n - 1], s[mainrow] = s[mainrow], s[n - 1]
    push!(Is, p[n - 1])
    push!(Js, n - 1)
    push!(Vs, 1.0)
    multiplier = A[p[n], n - 1] / A[p[n - 1], n - 1]
    A[p[n], n - 1] = 0
    push!(Is, p[n])
    push!(Js, n - 1)
    push!(Vs, multiplier)
    A[p[n], n] = A[p[n], n] - multiplier * A[p[n - 1], n]
    push!(Is, p[n])
    push!(Js, n)
    push!(Vs, 1.0)
    dropzeros!(A) #to improve the performance at the cost of memory, remove this line
    return (sparse(Is, Js, Vs), A, p)
end

#=
    Function that permutes the rows of the matrix
    Input:
        A - sparse matrix
        n - dimensions of the sparse matrix
        p - permutation
    Output:
        A::SparseMatrixCSC - permuted matrix
=#
function permuted(A::SparseMatrixCSC{Float64, Int64}, n::Int64, p::Vector{Int64})
    Is = Vector{Int64}()
    Js = Vector{Int64}()
    Vs = Vector{Float64}()

    rows = rowvals(A)
    vals = nonzeros(A)

    perm = Vector{Int64}(n);
    for i = 1:n
        perm[p[i]] = i;
    end
    for i = 1:n
       for j in nzrange(A, i)
          push!(Is, perm[rows[j]])
          push!(Js, i)
          push!(Vs, vals[j])
       end
    end
    return sparse(Is, Js, Vs)
end

#=
    Function that permutes the vector b
    Input:
        b - vector
        n - dimensions of the vector
        p - permutation
    Output:
        b::Vector{Float64} - permuted vector
=#
function permuted(b::Vector{Float64}, n::Int64, p::Vector{Int64})
    bperm = Vector{Float64}(n)
    for i = 1:n
        bperm[i] = b[p[i]]
    end

    return bperm
end

#=
    Function that solves a linear system LUx = b
    Input:
        L - lower triangular matrix
        U - upper triangular matrix
        n - dimensions of the matrix L and U
        l - dimensions of the block of the matrix L and U
        b - vector
    Output:
        x::Vector{Float64} - vector with solution of the linear system LUx = b
=#
function solve_linear_system(L::SparseMatrixCSC{Float64, Int64}, U::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
    for k = 1:n-2
        if k % l <= l - 2
            rows_to_set = l - k % l
        else
            rows_to_set = 2*l - k % l
        end
        for i = k+1:k+rows_to_set
            b[i] = b[i] - L[i, k]*b[k]
        end
    end
    b[n] = b[n] - L[n, n-1] * b[n - 1]

    b[n] = b[n] / U[n, n]
    for i in n-1:-1:1
		s = b[i]
        offset = min(n, i + l)
		for j = offset:-1:i+1
			s = s - U[i, j] * b[j]
		end
		b[i] = s /U[i, i]
	end
    return b
end

#=
    Function that solves a linear system LUx = b where LU was computed using partial pivoting
    Input:
        L - lower triangular matrix
        U - upper triangular matrix
        n - dimensions of the matrix L and U
        l - dimensions of the block of the matrix L and U
        b - vector
        p - permutation
    Output:
        x::Vector{Float64} - vector with solution of the linear system LUx = b
=#
function solve_linear_system(L::SparseMatrixCSC{Float64, Int64}, U::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64}, p::Vector{Int64})
    for k = 1:n-2
        if k % l <= l - 2
            rows_to_set = l - k % l
        else
            rows_to_set = 2*l - k % l
        end
        for i = k+1:k+rows_to_set
            b[p[i]] = b[p[i]] - L[i, k]*b[p[k]]
        end
    end
    b[p[n]] = b[p[n]] - L[n, n-1] * b[p[n - 1]]

    b[p[n]] = b[p[n]] / U[n, n]
    for i in n-1:-1:1
		sum = b[p[i]]
        offset = min(n, p[i] + 2*l + 1)
		for j = offset:-1:i+1
			sum = sum - U[i, j] * b[p[j]]
		end
		b[p[i]] = sum / U[i, i]
	end
    return b
end

#=
    A convenience function that computes LU factoriation and uses it to solve the linear system Ax = b
    Input:
        A - sparse matrix
        b - vector
        n - dimensions of the sparse matrix
        l - dimensions of the block of the sparse matrix
        pivoting - indicator whether the partial pivoting shall be used to compute the LU factorization
        filename - name of the file where the solution will be saved
    Output:
        nothing
=#
function solve_linear_system!(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64, pivoting = false::Bool, filename = "solution.txt"::String)
    if !pivoting
        (L, U) = lufactorization!(A, n, l)
        x = solve_linear_system(L, U, n, l, b)
    else
        (L, U, p) = lufactorization_pivoting!(A, n, l)
        L = permuted(L, n, p)
        U = permuted(U, n, p)
        x = solve_linear_system(L, U, n, l, b, p)
    end
    save_vector(x, filename)
end

#=
    Function that saves a vector to a file, it adds the relative error as the first line in the file, when the vector is computed instead of being read from a file
    Input:
        x - vector
        filename - name of the file
    Output:
        nothing
=#
function save_vector(x::Vector{Float64}, filename::String)
    open(filename, "w") do f
        length = size(x, 1)
        if bcomputed
            exactx = ones(length)
            relative_error = norm(exactx - x) / norm(exactx)
            println(f, relative_error)
        end
        for i = 1:length
            println(f, x[i])
        end
    end
end
end
