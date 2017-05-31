print("\nimporting modules")

from ops import *

IMAGENAME = "pq1.png"

MAXPOW = 7

FOURSIGMAS = [1]
FOURWINSIZE = 4
FOURNDIR = 2
FOURCOUNT = len(FOURSIGMAS) * FOURNDIR

FILTSIGMAS = [1]
FILTSIZE = 5
FILTNDIR = 2
FILTCOUNT = len(FILTSIGMAS) * FILTNDIR

LOCALISESIZE = 4

TOTALCOUNT = 3 * (FOURCOUNT * FILTCOUNT + 1)

EPS = [1] * TOTALCOUNT
COMBINE = [1] * TOTALCOUNT

fourIt = itertools.product(dirGene(FOURNDIR), FOURSIGMAS)
filtIt = itertools.product(dirGene(FILTNDIR), FILTSIGMAS)

#actual script

#>>> Opening image
image = openIm("test_images/" + IMAGENAME, MAXPOW)
#show(image)
psize(image, 3)
ptime(3)
pmem(3)

#>>> To LAB!
labImage = RGBtoLAB(image)
#show(cint(labImage))
ptime(3)
pmem(3)

#>>> Windows vectorizing
wimage = winVect(labImage, FOURWINSIZE, FOURWINSIZE)
psize(wimage, 3)
ptime(3)
pmem(3)

#>>> Multi fourier transforming
mfftimage = np.abs(fft(wimage))
psize(mfftimage, 3)
ptime(3)
pmem(3)


#>>> Measuring responses
pr("\nmeasuring responses", 2)
resp = lambda f: np.sum(multDim(mfftimage, f, [0,1]), axis = (0,1))
filtGene = (fgauss(FOURWINSIZE, u, sigma) for u, sigma in fourIt)
respGene = (resp(f) for f in filtGene)
ptime(3)
pmem(3)

#>>> Concatenating
combinedResp = concatGene(respGene, FOURCOUNT, 2)
ptime(3)
pmem(3)

#>>> Concatenate with original image
pr("\nconcatenating with original image", 2)
restoredCombinedResp = restoreShape(combinedResp, FOURWINSIZE // 2, FOURWINSIZE // 2)
combinedImage = np.concatenate((labImage, restoredCombinedResp), axis = 2)
psize(combinedImage, 3)
ptime(3)
pmem(3)

#>>>filtering
pr("\nFiltering and localising all channels", 2)
c = circleCut(FILTSIZE, FILTSIZE / 2)
loc = lambda arr, u, eps: localiseAndRestore(arr, LOCALISESIZE, u, eps)
gaussFilt = lambda u, sigma : filt(combinedImage, gaussd(FILTSIZE, u, sigma))
detectEdge = lambda u, sigma, eps: loc(gaussFilt(u, sigma), u, eps)
detect = lambda u, sigma, i : detectEdge(u, sigma, EPS[i]) * COMBINE[i]
filtImage = np.sum(detect(u, sigma, i) for i, (u, sigma) in enumerate(filtIt))
psize(filtImage, 3)
ptime(3)
pmem(3)

#>>> Reducing
for i in range(filtImage.shape[2]):
    comp = filtImage[:,:,i]
    filtImage[:,:,i] = red(comp)
    
#>>> Displaying
for i in range(filtImage.shape[2]):
    show(cint(filtImage[:,:,i]))
