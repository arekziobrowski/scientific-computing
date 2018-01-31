workspace()
include("zadanie1-4.jl")
using InterpolationUtils

x = [3.0, 1.0, 5.0, 6.0]
f = [1.0, -3.0, 2.0, 4.0]

#println(ilorazyRoznicowe(x, f)) #should be [1.0, 2.0, -0.375, 0.175]
#println(warNewton(x, ilorazyRoznicowe(x, f) , 3.0)) #should be 1
#println(naturalna([-1.0, 0.0, 1.0], [3.0, -7.0, 8.0, -6.0])) #should be [-4.0, 7.0, 8.0, -6.0]
println(ilorazyRoznicowe([-1.0, 0.0, 1.0, 2.0], [-1.0, 0.0, -1.0, 2.0]))
println(naturalna([-1.0, 0.0, 1.0], ilorazyRoznicowe([-1.0, 0.0, 1.0, 2.0], [-1.0, 0.0, -1.0, 2.0])))
