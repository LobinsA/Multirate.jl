using Multirate
using Printf
using PyPlot
using FFTW

plot = false

function wgn(n::Int; power::Real=1)
    return power * (randn(n) + randn(n)im)
end


function powerspectrum( x::Vector, window::Function = blackman )
    xLen = length( x )
    xWin = x .* window( xLen )
    10*log10.(fftshift(abs2.(fft( xWin ))))
end

function linspace(start,stop,len)
    collect(range(start,stop=stop,length=len))
end

function channelizerplots( signal::Vector, channelizedSignals::Matrix )
    (Nchannels,samplesPerChannel) = size( channelizedSignals )

    # Compute the spectrum of the original signal and add it to the table
    subplot(4,1,1)
    signalSpectrum     = powerspectrum( signal )
    plot( linspace(-0.5,0.5,length(signalSpectrum)), signalSpectrum )
    title( "Original Signal" )
    ylim(-10,60)
    
    subplot(4,1,2)
    plot( Base.axes(signal,1), real(signal), "-b", Base.axes(signal,1), imag(signal), "-r" )
    
    freqs = linspace( -0.5, 0.5, samplesPerChannel )
    t     = linspace( 0, Nchannels*samplesPerChannel-1, samplesPerChannel )
    
    # Plot the spectrum/time-domain traces of each channel and add them to the subtable
    for channel in 1:Nchannels
        thisSignal          = channelizedSignals[channel,:]
        spectrum            = powerspectrum( thisSignal )
        subplot(4,Nchannels,Nchannels*2+channel)
        plot( freqs, spectrum ); ylim(-10,60)
        title( "Channel $channel" )
        subplot(4,Nchannels,Nchannels*3+channel)
        plot( t, real(thisSignal), "-b", t, imag(thisSignal), "-r" )
    end
    
    show()
end

function segments(n; alpha=0.5)
    s0 = 1; n0 = n; segs = []
    
    while n > 0
        s=min(n, round(Int, alpha*rand()*n0))
        
        if s > 0
            push!(segs, s0:s0+s-1)
            n -= s; s0 += s
        end
    end
    
    return segs
end

const Th = Float32
const Tx = ComplexF32 # Datatype for x
const ƒs = 1.0        # Input sample rate

const Nchannels = 5
const samplesPerChannel = 5000

# Construct a linear chirp signal from ƒ = -0.5 to -0.5
const n                   = Nchannels * samplesPerChannel
const t                   = 0:n-1

const ψ = π * (t / n .- 1) .* t

const signal             = Array{Tx}(exp.(ψ * 1im / ƒs) + wgn( n, power=0.1))


for i in 1:10
    channelizer0        = Channelizer( Th, Tx, Nchannels, 32 )
    @time filt( channelizer0, signal )
end

const channelizer1        = Channelizer( Th, Tx, Nchannels, 32 )
@time const channelizedSignals1 = filt( channelizer1, signal )

const channelizer2        = Channelizer( Th, Tx, Nchannels, 32 )
global channelizedSignals2 = Array{Tx}(undef, Nchannels, 0)
for segment in segments(samplesPerChannel)
    expandedSegment = (segment.start-1)*Nchannels+1 : segment.stop*Nchannels
    @printf("Channelizing segment %s as %s\n", segment, expandedSegment)
    @time global channelizedSignals2 = [channelizedSignals2 filt( channelizer2, signal[expandedSegment] )]
end

if isapprox(channelizedSignals1, channelizedSignals2, atol=sqrt(Nchannels*samplesPerChannel*eps(Th)))
    @printf("Results match\n")
else
    @printf("Results don't match!\n")
end

# Create the table of plots
if plot
    table = channelizerplots( signal, fftshift(channelizedSignals1,1) )
end

