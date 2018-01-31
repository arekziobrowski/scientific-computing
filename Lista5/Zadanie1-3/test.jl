#Obliczenia naukowe lista 5, program testujacy dla zadan 1-3
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

workspace()
include("blocksys.jl")
using Blocksys
include("matrixgen.jl")
using matrixgen

#show(IOContext(STDOUT), "text/plain", full(A)) #we can pretty print the matrix using this function (it takes O(n^2) memory, so be careful for large n-s)

#tests for n = 16

(A, n, l) = load_matrix("../data/A.txt")
b = compute_vector(A, n)
x = gaussian_elimination!(A, b, n, l)
save_vector(x, "../solutions/gaussian16.txt")

(A, n, l) = load_matrix("../data/A.txt")
b = compute_vector(A, n)
x = gaussian_elimination_pivoting!(A, b, n, l)
save_vector(x, "../solutions/gaussianpiv16.txt")

(A, n, l) = load_matrix("../data/A.txt")
b = compute_vector(A, n)
solve_linear_system!(A, b, n, l, false, "../solutions/lufact16.txt")

(A, n, l) = load_matrix("../data/A.txt")
b = compute_vector(A, n)
solve_linear_system!(A, b, n, l, true, "../solutions/lufactpiv16.txt")

#tests for n = 10000

(A, n, l) = load_matrix("../data/A10k.txt")
b = compute_vector(A, n)
x = gaussian_elimination!(A, b, n, l)
save_vector(x, "../solutions/gaussian10k.txt")

(A, n, l) = load_matrix("../data/A10k.txt")
b = compute_vector(A, n)
x = gaussian_elimination_pivoting!(A, b, n, l)
save_vector(x, "../solutions/gaussianpiv10k.txt")

(A, n, l) = load_matrix("../data/A10k.txt")
b = compute_vector(A, n)
solve_linear_system!(A, b, n, l, false, "../solutions/lufact10k.txt")

(A, n, l) = load_matrix("../data/A10k.txt")
b = compute_vector(A, n)
solve_linear_system!(A, b, n, l, true, "../solutions/lufactpiv10k.txt")

#tests for n = 50000

(A, n, l) = load_matrix("../data/A50k.txt")
b = compute_vector(A, n)
x = gaussian_elimination!(A, b, n, l)
save_vector(x, "../solutions/gaussian50k.txt")

(A, n, l) = load_matrix("../data/A50k.txt")
b = compute_vector(A, n)
x = gaussian_elimination_pivoting!(A, b, n, l)
save_vector(x, "../solutions/gaussianpiv50k.txt")

(A, n, l) = load_matrix("../data/A50k.txt")
b = compute_vector(A, n)
solve_linear_system!(A, b, n, l, false, "../solutions/lufact50k.txt")

(A, n, l) = load_matrix("../data/A50k.txt")
b = compute_vector(A, n)
solve_linear_system!(A, b, n, l, true, "../solutions/lufactpiv50k.txt")
