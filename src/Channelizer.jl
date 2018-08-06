import Multirate: PFB, shiftin!

# using a different coefficient layout from that used in Filters.jl
# this layout should allow for a simpler expression of the filter operation
# with less memory movement
function taps2pfb2{T}( h::Vector{T}, N𝜙::Integer )
    hLen     = length( h )
    stuffed  = [ h; zeros(T, mod(-hLen, N𝜙))]
    pfb       = reshape(flipdim(stuffed, 1), N𝜙, :)
    
    return pfb
end

# Interpolator FIR kernel
type Channelizer{Th, Tx}
    pfb::PFB{Th}
    h::Vector{Th}
    Nchannels::Int
    tapsPer𝜙::Int
    history::Matrix{Tx}
end

function Channelizer{Th}( Tx, h::Vector{Th}, Nchannels::Integer )
    pfb       = taps2pfb2( h, Nchannels )
    Nchannels = size( pfb )[1]
    tapsPer𝜙  = size( pfb )[2]
    history   = zeros(Tx, Nchannels, tapsPer𝜙 - 1)
    Channelizer( pfb, h, Nchannels, tapsPer𝜙, history)
end

function Channelizer( Th, Tx, Nchannels::Integer, tapsPer𝜙 = 20 )
    hLen = tapsPer𝜙 * Nchannels
    h    = firdes( hLen, 0.45/Nchannels, kaiser ) .* Nchannels
    Channelizer( Tx, Vector{Th}(h), Nchannels )
end




function filt!{Tb,Th,Tx}( output::Matrix{Tb}, kernel::Channelizer{Th, Tx}, x::Matrix{Tx} )
    Nchannels         = kernel.Nchannels
    tapsPer𝜙          = kernel.tapsPer𝜙
    pfb               = kernel.pfb
    outLen            = size(output,2)
    history           = kernel.history
    histLen           = tapsPer𝜙-1
    fftBuf            = Array{Tb}(Nchannels)

    @assert size(x)         == size(output)
    @assert size(x,1)       == Nchannels
    @assert size(history)   == (Nchannels, histLen)
    @assert Tb              == promote_type(Th,Tx)

    # initial segment using history from previous call
    @inbounds xh = [history x[:,1:histLen]]
    @simd for s in 1:min(outLen, histLen)
        @inbounds fftBuf[:] = flipdim(sum( xh[:,s:s+histLen] .* pfb, 2 ), 1)
        @inbounds output[:,s] = fftshift(ifft(fftBuf))
    end

    # history-independent portion
    @simd for s in histLen+1:outLen
        # @inbounds fftBuf[:] = flipdim(sum( x[:,s-histLen:s] .* pfb, 2 ), 1)
        @simd for k in 1:Nchannels
            fftElem = zero(Tb)
            
            @simd for m in 0:histLen
               @inbounds fftElem += x[k, s-histLen+m] * pfb[k, 1+m]
            end
            
            fftBuf[Nchannels+1-k] = fftElem
        end
        
        @inbounds output[:,s] = fftshift(ifft(fftBuf))
    end

    # set history for next call
    kernel.history = shiftin!( history, x );

    return output
end

function filt{Th,Tx}( kernel::Channelizer{Th, Tx}, x::Vector{Tx} )
    xLen   = length( x )
    @assert xLen % kernel.Nchannels == 0
    xm = reshape(x, kernel.Nchannels, Int(xLen/kernel.Nchannels));

    output = similar(xm, promote_type(Th,Tx))
    return filt!( output, kernel, xm )
end
