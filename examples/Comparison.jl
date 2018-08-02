import DSP
import Multirate

sampleRate    = 48000
interpolation = 147
decimation    = 160
ratio         = interpolation  // decimation
numTaps       = 24*interpolation
x             = rand( 1_000_000 )
h             = Multirate.firdes( numTaps, 0.5/interpolation, Multirate.kaiser, beta = 7.8562  )

function naiveresampler{T}( h::Vector{T}, x::Vector{T}, ratio::Rational{Int} )
    upfactor   = numerator( ratio )
    downfactor = denominator( ratio )
    xStuffed   = zeros( T, length(x) * upfactor )

    for n = 0:length(x)-1;
        xStuffed[ n*upfactor+1 ] = x[ n+1 ]
    end

    yInterpolated = DSP.filt( h, xStuffed )
    y = [ yInterpolated[n] for n = 1:downfactor:length( yInterpolated ) ]
end

@time y = naiveresampler( h, x, 147//160 );
@time y = Multirate.filt( h, x, ratio );
