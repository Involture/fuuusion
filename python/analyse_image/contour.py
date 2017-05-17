print("\nimporting modules")

from ops import *
from filts import *

#actual script

p("\nimporting image")
image = openIm("test_images/pq1.png")
psize(image)
ptime()
pmem()

#p("\nfourier transforming the image")
#fftimage = fft(image)
#psize(fftimage)
#ptime()
#pmem()
#show(sf(fftimage, 10))

p("\nwindows vectorizing the image")
wimage = winVect(image, 8, 8)
psize(wimage)
ptime()
pmem()

p("\nmulti fourier transforming the image")
mfftimage = fft(wimage)
psize(mfftimage)
ptime()
pmem()

p("\nmeasuring response")
trans = (-2, -1) + tuple(range(2, wimage.ndim - 2)) + (0, 1)
filt = vgauss(8, 50) * circleCut(8,8,100)
mult = np.transpose(filt * np.transpose(mfftimage, trans), trans)
resp = np.sum(np.abs(mult), axis = (0,1))
reducedResp = bToL(ubint(resp), 99)
show(resp)
psize(resp)
ptime()
pmem()

p("\nisp vectorizing the image")
ispimage = ispVect(image)
psize(ispimage)
ptime()
pmem()

p("\nsetting filter")
f = expand(vgrad(8), 2)
psize(f)
ptime()
pmem()

p("\nfiltering")
ispgradimage = ispFilt(ispimage, f)
psize(ispgradimage)
ptime()
pmem()
