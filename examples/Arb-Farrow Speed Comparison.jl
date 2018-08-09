import Multirate
import DSP
using Printf

function time_firfarrow( self::Multirate.FIRFilter{Multirate.FIRFarrow{Th}}, x::Vector{Tx}, resampleRate ) where {Th,Tx}
    kernel = self.kernel
    xLen   = length( x )
    println( "\nFIRFarrow speed test" )
    @printf( "\tresampling rate  %f\n", resampleRate )
                                                             # @printf( "\tpolynomial order %d\n", polyorder )
    @printf( "\tx type           %s\n", string(Tx) )
    @printf( "\tx length         %d\n", xLen )
    @printf( "\th type           %s\n", string(Th) )
    @printf( "\th length         %d\n", kernel.tapsPer𝜙*kernel.N𝜙 )
    @printf( "\tN𝜙               %d\n", kernel.N𝜙 )
    @printf( "\ttaps per 𝜙       %d\n", kernel.tapsPer𝜙 )
    ( y, elapsed, allocated, z ) = @timed DSP.filt( self, x )
    @printf( "\telapsed time (s) %1.3f\n", elapsed )
    @printf( "\tinput samples/s  %1.3e\n", xLen/elapsed )
    @printf( "\toutput samples/s %1.3e\n", length(y)/elapsed )
end


function time_firarbitrary( self::Multirate.FIRFilter{Multirate.FIRArbitrary{Th}}, x::Vector{Tx}, resampleRate ) where {Th,Tx}
    println( "\nFIRArbitrary Speed Test" )
    @printf( "\tresampling rate  %f\n", resampleRate )
    @printf( "\tx type           %s\n", string(Tx) )
    @printf( "\tx length         %d\n", xLen )
    @printf( "\th type           %s\n", string(Th) )
    @printf( "\th length         %d\n", hLen )
    @printf( "\tN𝜙               %d\n", N𝜙 )
    @printf( "\ttaps per 𝜙       %d\n", tapsPer𝜙 )
    ( y, elapsed, allocated, z ) = @timed DSP.filt( self, x )    
    @printf( "\telapsed time (s) %1.3f\n", elapsed )
    @printf( "\tinput samples/s  %1.3e\n", xLen/elapsed )
    @printf( "\toutput samples/s %1.3e\n", length(y)/elapsed )
end

N𝜙           = 32                                            # Number of polyphase partitions
tapsPer𝜙     = 10                                            # N𝜙 * tapsPer𝜙 will be the length of out protoyTimepe filter taps
polyorder    = 4                                             # Our taps will tranformed into
Th           = Float32
hLen         = tapsPer𝜙*N𝜙                                   # Total number of filter taps
xLen         = 10_000_000                                    # Number of signal samples

for resampleRate in ( 1.0, 1/2.123456789 ), Tx in ( Float32, Float64, ComplexF32, ComplexF64 )
    cutoffFreq = min( 0.45/N𝜙, resampleRate/N𝜙 )               # N𝜙 is also the integer interpolation, so set cutoff frequency accordingly
    h          = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser ) .* N𝜙  # Generate filter taps and scale by polyphase interpolation (N𝜙)
    
    farrowfilt = Multirate.FIRFilter( h, resampleRate, N𝜙, polyorder ) # Construct a FIRFilter{FIRFarrow} object
    arbfilt    = Multirate.FIRFilter( h, resampleRate, N𝜙 )
    x          = rand( Tx, xLen )
    time_firfarrow( farrowfilt, x, resampleRate )
    time_firarbitrary( arbfilt, x, resampleRate )
end
