import Multirate: PFB, shiftin!
using FFTW

# using a different coefficient layout from that used in Filters.jl
# this layout should allow for a simpler expression of the filter operation
# with less memory movement
function taps2pfb2( h::Vector{T}, N𝜙::Integer ) where T
    hLen     = length( h )
    stuffed  = [ h; zeros(T, mod(-hLen, N𝜙))]
    pfb       = reshape(reverse(stuffed, dims=1), N𝜙, :)
    
    return pfb
end

# Interpolator FIR kernel
mutable struct Channelizer{Th, Tx}
    pfb::PFB{Th}
    h::Vector{Th}
    Nchannels::Int
    tapsPer𝜙::Int
    history::Matrix{Tx}
    ifft!_plan
end

function Channelizer( Tx, h::Vector{Th}, Nchannels::Integer ) where Th
    pfb       = taps2pfb2( h, Nchannels )
    Nchannels = size( pfb )[1]
    tapsPer𝜙  = size( pfb )[2]
    history   = zeros(Tx, Nchannels, tapsPer𝜙 - 1)
    ifft!_plan = plan_ifft!(Array{promote_type(Tx,Th)}(undef, Nchannels))
    Channelizer( pfb, h, Nchannels, tapsPer𝜙, history, ifft!_plan )
end

function Channelizer( Th, Tx, Nchannels::Integer, tapsPer𝜙 = 20 )
    hLen = tapsPer𝜙 * Nchannels
    h    = firdes( hLen, 0.45/Nchannels, kaiser ) .* Nchannels
    Channelizer( Tx, Vector{Th}(h), Nchannels )
end




function filt!( output::Matrix{Tb}, kernel::Channelizer{Th, Tx}, x::AbstractMatrix{Tx} ) where {Tb,Th,Tx}
    Nchannels         = kernel.Nchannels
    tapsPer𝜙          = kernel.tapsPer𝜙
    pfb               = kernel.pfb
    outLen            = size(output,2)
    history           = kernel.history
    histLen           = tapsPer𝜙-1

    @assert size(x)         == size(output)
    @assert size(x,1)       == Nchannels
    @assert size(history)   == (Nchannels, histLen)
    @assert Tb              == promote_type(Th,Tx)

    # initial segment using history from previous call
    @inbounds xh = [history x[:,1:histLen]]
    @simd for s in 1:min(outLen, histLen)
        @inbounds output[:,s] = reverse(sum( view(xh,:,s:s+histLen) .* pfb, dims=2 ), dims=1)
    end

    # history-independent portion
    @simd for s in histLen+1:outLen
#        @inbounds output[:,s] = reverse(sum( view(x,:,s-histLen:s) .* pfb, dims=2 ), dims=1)
        @simd for k in 1:Nchannels
            fftElem = zero(Tb)
            
            @simd for m in 0:histLen
               @inbounds fftElem += x[k, s-histLen+m] * pfb[k, 1+m]
            end
            
            output[Nchannels+1-k,s] = fftElem
        end
    end

    for s in 1:outLen
        @inbounds kernel.ifft!_plan * view(output,:,s)
    end

    # set history for next call
    kernel.history = shiftin!( history, x );

    return output
end

function filt( kernel::Channelizer{Th, Tx}, x::AbstractVector{Tx} ) where {Th,Tx}
    xLen   = length( x )
    @assert xLen % kernel.Nchannels == 0
    xm = reshape(x, kernel.Nchannels, Int(xLen/kernel.Nchannels));

    output = similar(xm, promote_type(Th,Tx))
    return filt!( output, kernel, xm )
end
