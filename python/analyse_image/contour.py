VERBOSE = True

###

p("\nimporting modules")

from ops import *
from filts import *

#actual script

p("\nimporting image")
camille = misc.imread("camille.png")[:1024, :1024,:]
psize(camille)
ptime()
pmem()

#p("\nfourier transforming the image")
#fftcamille = fft(camille)
#psize(fftcamille)
#ptime()
#pmem()
#show(sf(fftcamille, 10))

#p("\nwindows vectorizing the image")
#wcamille = winVect(camille, 4, 4)
#psize(wcamille)
#ptime()
#pmem()
#
#p("\nmulti fourier transforming the image")
#mfftcamille = mfft(wcamille)
#psize(mfftcamille)
#ptime()
#pmem()
#
#p("\nmeasuring response")
#resp = red(ubint(np.sum(np.abs(filt * mfftcamille), axis = (0,1))))

p("\nisp vectorizing the image")
ispcamille = ispVect(camille)
psize(ispcamille)
ptime()
pmem()

p("\nsetting filter")
f = expand(vgrad(8), 2)
psize(f)
ptime()
pmem()

p("\nfiltering")
ispgradcamille = ispFilt(ispcamille, f)
psize(ispgradcamille)
ptime()
pmem()
