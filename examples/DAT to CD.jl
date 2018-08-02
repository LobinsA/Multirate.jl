# Borrowed from http://www.mathworks.com/help/signal/ref/upfirdn.html

import DSP: kaiser
import Multirate: firdes, filt
import PyPlot

sampleRate    = 48000
interpolation = 147
decimation    = 160
ratio         = interpolation  // decimation
Nt            = 24*interpolation

h = firdes( Nt, 0.5/interpolation, kaiser, beta = 7.8562  ) .* interpolation

n = 0:10239
x = sin.(2*pi*1e3/sampleRate*n)
y = filt( h, x, ratio )

tOriginal = n[1:49]./sampleRate
xOriginal = x[1:49]
tResamp   = n[1:45]./(sampleRate*ratio)
xResamp   = y[12:56]

PyPlot.stem( tOriginal, xOriginal, linefmt = "b-", markerfmt = "b." )
PyPlot.stem( tResamp, xResamp, linefmt = "r-", markerfmt = "r." )
PyPlot.show()