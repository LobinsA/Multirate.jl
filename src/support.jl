
# Computes the dot product of a single column of a, specified by aColumnIdx, with the vector b.
# The number of elements used in the dot product determined by the size(A)[1].
# Note: bIdx is the last element of b used in the dot product.
function unsafedot( a::Matrix, aColIdx::Integer, b::Vector, bLastIdx::Integer )
    aLen     = size(a)[1]
    bBaseIdx = bLastIdx - aLen
    dotprod  = a[ 1, aColIdx ] * b[ bBaseIdx + 1 ]
    @simd for i in 2:aLen
        @inbounds dotprod += a[ i, aColIdx ] * b[ bBaseIdx + i ]
    end

    return dotprod
end

function unsafedot( a::Matrix, aColIdx::Integer, b::Vector{T}, c::Vector{T}, cLastIdx::Integer ) where T
    aLen = size(a)[1]
    bLen = length(b)
    bLen == aLen-1  || error( "length(b) must equal to length(a)[1] - 1" )
    cLastIdx < aLen || error( "cLastIdx but be < length(a)")

    dotprod = a[ 1, aColIdx] * b[ cLastIdx ]
    @simd for i in 2:aLen-cLastIdx
        @inbounds dotprod += a[ i, aColIdx ] * b[ i+cLastIdx-1 ]
    end
    @simd for i in 1:cLastIdx
        @inbounds dotprod += a[ aLen-cLastIdx+i, aColIdx ] * c[ i ]
    end

    return dotprod
end

function unsafedot( a::Vector, b::Vector, bLastIdx::Integer )
    aLen     = length(a)
    bBaseIdx = bLastIdx - aLen
    dotprod  = a[ 1 ] * b[ bBaseIdx + 1 ]
    @simd for i in 2:aLen
        @inbounds dotprod += a[ i ] * b[ bBaseIdx + i ]
    end

    return dotprod
end

function unsafedot( a::Vector, b::Vector{T}, c::Vector{T}, cLastIdx::Integer ) where T
    aLen    = length(a)
    dotprod = zero(a[1]*b[1])
    @simd for i in 1:aLen-cLastIdx
        @inbounds dotprod += a[ i ] * b[ i+cLastIdx-1 ]
    end
    @simd for i in 1:cLastIdx
        @inbounds dotprod += a[ aLen-cLastIdx+i ] * c[ i ]
    end

    return dotprod
end



# Shifts b into the end a.
# a = [ a, b ][1:length(a)]
function shiftin!( a::Vector{T}, b::Vector{T} ) where T
    aLen = length( a )
    bLen = length( b )

    if bLen >= aLen
        copyto!( a, 1, b, bLen - aLen + 1, aLen )
    else

        for i in 1:aLen-bLen
            @inbounds a[i] = a[i+bLen]
        end
        bIdx = 1
        for i in aLen-bLen+1:aLen
            @inbounds a[i] = b[bIdx]
            bIdx += 1
        end
    end

    return a
end

function shiftin!( dst::AbstractMatrix{T}, src::AbstractMatrix{T} ) where T
    @assert size( dst, 1 ) == size( src, 1 )
    width = size( dst, 1 )
    
    dstLen = size( dst, 2 )
    srcLen = size( src, 2 )

    if srcLen >= dstLen
        srcRange = CartesianIndices((1:width, srcLen-dstLen+1:srcLen))
        dstRange = CartesianIndices((1:width, 1:dstLen))
        copyto!( dst, dstRange, src, srcRange )
    else

        for i in 1:dstLen-srcLen
            @inbounds dst[:,i] = dst[:,i+srcLen]
        end
        srcIdx = 1
        for i in dstLen-srcLen+1:dstLen
            @inbounds dst[:,i] = src[:,srcIdx]
            srcIdx += 1
        end
    end

    return dst
end




function polyfit( y::AbstractVector, polyorder::Integer )
  A = [ x^p for x in 1:length(y), p = 0:polyorder ]
  Poly(A \ y)
end