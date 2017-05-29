print("\nimporting modules")

from ops import *

MAXPOW = 8
WINSIZE = 8
SIGMAS = [0.25, 0.5, 1, 2, 4, 8]
NDIR = 8
IMAGENAME = "pq1.png"
FILTSIZE = 5
EPS = 1

#actual script

p("\nimporting image")
image = openIm("test_images/" + IMAGENAME, MAXPOW)
#show(image)
psize(image)
ptime()
pmem()

p("\nto lab !!!")
labImage = RGBtoLAB(image)
#show(cint(labImage))
ptime()
pmem()

#p("\nfourier transforming")
#fimage = fft(labImage)
#show(sf(np.abs(fimage), 1000))
#
#p("\nreverse fourier transforming")
#rfimage = ifft(fimage)
#show(cint(np.abs(rfimage)))

p("\nwindows vectorizing the image")
wimage = winVect(labImage, WINSIZE, WINSIZE)
psize(wimage)
ptime()
pmem()

p("\nmulti fourier transforming the image")
mfftimage = np.abs(fft(wimage))
psize(mfftimage)
ptime()
pmem()

p("\nmeasuring responses")
resp = lambda f: np.sum(multDim(mfftimage, f, [0,1]), axis = (0,1))
respGene = lambda sigma: (resp(f) for f in fGaussGene(WINSIZE, NDIR, sigma))
summedResp = lambda sigma: sum(r for r in respGene(sigma))
summedRespList = [summedResp(sigma) for sigma in SIGMAS]
ptime()
pmem()

#p("\nreducing and displaying responses")
#red = lambda arr : cint(arr * 255 / np.max(arr))
#for r in summedRespList:
#    show(red(r))
#

p("\nconcatenating")
combinedResp = np.concatenate(summedRespList, axis = 2)
ptime()
pmem()

p("\nstacking with original image")
restoredCombinedResp = restoreShape(combinedResp, WINSIZE // 2, WINSIZE // 2)
combinedImage = np.concatenate((labImage, restoredCombinedResp), axis = 2)
p(combinedImage.shape)
psize(combinedImage)
ptime()
pmem()

#p("\nisp vectorizing the image")
#ispimage = ispVect(combinedImage)
#psize(ispimage)
#ptime()
#pmem()
#

#p("\nfiltering")
#ispgradimage = ispFilt(ispimage, f)
#psize(ispgradimage)
#ptime()
#pmem()

p("\nfiltering")
c = circleCut(FILTSIZE, FILTSIZE / 2)
loc = lambda arr, u: localise(arr, WINSIZE, u, EPS)
gaussFilt = lambda u, sigma: filt(combinedImage, gaussd(FILTSIZE, u, sigma))
detectEdge = lambda u, sigma: loc(gaussFilt(u, sigma), u)
it = itertools.product(dirGene(NDIR), SIGMAS)
filtImage = np.sum(detectEdge(u, sigma) for u, sigma in it)
psize(filtImage)
ptime()
pmem()

p("\nreducing")
for i in range(filtImage.shape[2]):
    comp = np.abs(filtImage[:,:,i])
    filtImage[:,:,i] = comp * 255 / comp.max()
    
p("\ndisplaying")
for i in range(filtImage.shape[2]):
    show(cint(filtImage[:,:,i]))
