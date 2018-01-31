#Obliczenia naukowe lista 1, zadanie 1
#Autor: Arkadiusz Ziobrowski, numer indeksu: 229728

"""
    Macheps functions for all floating point types
    Parameters: none
    Returns: macheps<FloatX>
"""
function macheps16()
    x = Float16(1.0)
    while Float16(1.0) + x / Float16(2.0) > Float16(1.0)
        x = x / Float16(2.0)
    end
    return x
end

function macheps32()
    x = Float32(1.0)
    while Float32(1.0) + x / Float32(2.0) > Float32(1.0)
        x = x / Float32(2.0)
    end
    return x
end

function macheps64()
    x = Float64(1.0)
    while Float64(1.0) + x / Float64(2.0) > Float64(1.0)
        x = x / Float64(2.0)
    end
    return x
end

"""
    Eta functions for all floating point types
    Parameters: none
    Returns: eta<FloatX>
"""
function eta16()
    x = Float16(1.0)
    while x / Float16(2.0) > Float16(0.0)
        x = x / Float16(2.0)
    end
    return x
end

function eta32()
    x = Float32(1.0)
    while x / Float32(2.0) > Float32(0.0)
        x = x / Float32(2.0)
    end
    return x
end

function eta64()
    x = Float64(1.0)
    while x / Float64(2.0) > Float64(0.0)
        x = x / Float64(2.0)
    end
    return x
end

"""
    Max functions for all floating point types
    Parameters: none
    Returns: max<FloatX>
"""
function max16()
    x = prevfloat(Float16(1.0))
    while !isinf(x * Float16(2.0))
        x = x * Float16(2.0)
    end
    return x
end

function max32()
    x = prevfloat(Float32(1.0))
    while !isinf(x * Float32(2.0))
        x = x * Float32(2.0)
    end
    return x
end

function max64()
    x = prevfloat(Float64(1.0))
    while !isinf(x * Float64(2.0))
        x = x * Float64(2.0)
    end
    return x
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

println("macheps16() = ", macheps16(), ", eps(Float16) = ", eps(Float16))
println("macheps32() = ", macheps32(), ", eps(Float32) = ", eps(Float32))
println("macheps64() = ", macheps64(), ", eps(Float64) = ", eps(Float64))

println("eta16() = ", eta16(), ", nextfloat(Float16(0.0)) = ", nextfloat(Float16(0.0)))
println("eta32() = ", eta32(), ", nextfloat(Float32(0.0)) = ", nextfloat(Float32(0.0)))
println("eta64() = ", eta64(), ", nextfloat(Float64(0.0)) = ", nextfloat(Float64(0.0)))

println("max16() = ", max16(), ", realmax(Float16) = ", realmax(Float16))
println("max32() = ", max32(), ", realmax(Float32) = ", realmax(Float32))
println("max64() = ", max64(), ", realmax(Float64) = ", realmax(Float64))
