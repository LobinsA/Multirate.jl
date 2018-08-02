import Multirate
import Multirate: NaiveResamplers
using Base.Test

Tx           = Float32
Th           = Float32
numFilters   = 32
x            = Array{Tx}(1:101)
resampleRate = 33/32

cutoffFreq      = 0.45
transitionWidth = 0.05
(hLen, β)       = Multirate.kaiserlength( transitionWidth, samplerate = numFilters )
hLen            = ceil( Int, hLen/numFilters ) .* numFilters
h               = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = numFilters, beta = β ) .* numFilters
h               = convert( Vector{Th}, h)

@time yNaive = NaiveResamplers.naivefilt( h, x, resampleRate, numFilters )
@time yArb   = Multirate.filt( h, x, resampleRate )

commonLen = min( length(yNaive), length(yArb) )

yNaiveCommon = yNaive[1:commonLen]
yArbCommon = yArb[1:commonLen]

isapprox( yNaiveCommon, yArbCommon, atol=eps() ) ||  display( [ (1:commonLen) yNaiveCommon yArbCommon (yNaiveCommon.-yArbCommon)] )
