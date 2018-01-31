#Obliczenia naukowe lista 1, zadanie 2
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

"""
    Calculate macheps using Kahan's macheps equation: 3.0 * (4.0/3.0 - 1.0) - 1.0
    Parameters: T - type
    Returns:
"""
function kahan_eps(T)
        return T(3.0) * (T(4.0) / T(3.0) - T(1.0)) - T(1.0)
end

"""
    Prints step by step execution of Kahan's equation using pretty_print_bits
    Parameters: T - type
    Returns:
"""
function kahan_eps_steps(T)
    print("4.0 / 3.0:\t\t\t")
    pretty_print_bits(T(4.0) / T(3.0))
    print("4.0 / 3.0 - 1.0:\t\t")
    pretty_print_bits(T(4.0) / T(3.0) - T(1.0))
    print("3.0 * (4.0 / 3.0 - 1.0):\t")
    pretty_print_bits(T(3.0) * (T(4.0) / T(3.0) - T(1.0)))
    print("3.0 * (4.0 / 3.0 - 1.0) - 1.0:\t")
    pretty_print_bits(T(3.0) * (T(4.0) / T(3.0) - T(1.0)) - T(1.0))
    print("eps():\t\t\t\t")
    pretty_print_bits(eps(T))
end

"""
    Function to pretty print floating type input bits
    Parameters: x - floating point number
    Returns:
"""
function pretty_print_bits(x)
    if(typeof(x) == Float16)
        str = String(bits(x))
        println(str[1], " ", str[2:6], " ", str[7:16])
    elseif (typeof(x) == Float32)
        str = String(bits(x))
        println(str[1], " ", str[2:9], " ", str[10:32])
    elseif (typeof(x) == Float64)
        str = String(bits(x))
        println(str[1], " ", str[2:12], " ", str[13:64])
    end
end

println("kahan_eps(Float16) = ", kahan_eps(Float16), ", eps(Float16) = ", eps(Float16))
println("kahan_eps(Float32) = ", kahan_eps(Float32), ", eps(Float32) = ", eps(Float32))
println("kahan_eps(Float64) = ", kahan_eps(Float64), ", eps(Float64) = ", eps(Float64))

kahan_eps_steps(Float32)
