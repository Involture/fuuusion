from ops import *

fileList = []

MAXPOW = 10
   
FOURSIGMAS =  [8]
FOURNSIGMAS = len(FOURSIGMAS)
FOURWINSIZE = 16
FOURNDIR = 4
FOURCOUNT = FOURNSIGMAS  * FOURNDIR
   
FILTSIGMAS = [1, 2]
FILTNSIGMAS = len(FILTSIGMAS)
FILTSIZE = 5
FILTNDIR = 8
FILTCOUNT = FILTNSIGMAS * FILTNDIR
   
LOCALISESIZE = 8
   
TOTALCOUNT = 3 * ((FOURCOUNT + 1) * FILTCOUNT)
   
EPS = [1]
COMBINE = np.array([[[.1, .3, .3],
                    [.1, .3, .3],
                    [.1, .3, .3],
                    [.1, .3, .3],
                    [.3, 2., 2.]
                    ],
                   [[.05, .15, .15],
                    [.05, .15, .15],
                    [.05, .15, .15],
                    [.05, .15, .15],
                    [.15, 1., 1.]]
                      ])
   
fourIt = itertools.product(dirGene(FOURNDIR), FOURSIGMAS)


def __main__(IMAGENAME):
    
    image = openIm("test_images/" + IMAGENAME, MAXPOW)
    labImage = RGBtoLAB(image)
    np.save(IMAGENAME + "_lab", labImage)

    
    fftImage = fft(labImage)
    np.save(IMAGENAME + "_fft", fftImage)
    

    wimage = winVect(labImage, FOURWINSIZE, FOURWINSIZE)
    mfftimage = np.abs(fft(wimage))
    

    pr("\nmeasuring responses", 2)
    filtGene = (fgauss(FOURWINSIZE, u, sigma) for u, sigma in fourIt)
    resp = lambda f: np.sum(multDim(mfftimage, f, [0,1]), axis = (0,1))
    respGene = (resp(f) for f in filtGene)
    
    combinedResp = stackGene(respGene, FOURCOUNT, 2)
    np.save(IMAGENAME + "_resp", combinedResp[:,:,0,0])
    

    pr("\nconcatenating with original image", 2)
    t = FOURWINSIZE // 2
    restoredCombinedResp = restoreShape(combinedResp, t, t)
    expLabImage = np.expand_dims(labImage, axis =2)
    combinedImage = np.append(restoredCombinedResp, expLabImage,  axis = 2)
    
    
    pr("\nFiltering and localising all channels", 2)

    c = circleCut(FILTSIZE, FILTSIZE / 2)

    def detect(u, sigma):
        f = gaussd(FILTSIZE, u, sigma) * c
        filtered = filt(combinedImage, f)
        localised = localise(filtered, LOCALISESIZE, u, EPS)
        if u == (1., 0.) and sigma == 1.:
            np.save(IMAGENAME + "_filtered", filtered[:,:,0,0])
            np.save(IMAGENAME + "_localised", localised[:,:,0,0])
        return localised

    def summedDetect(sigma):
        return np.sum((detect(u, sigma) for u in dirGene(FILTNDIR)), axis = 2)

    filtImageGene = (summedDetect(sigma) for sigma in FILTSIGMAS)
    filtImage = stackGene(filtImageGene, FILTNSIGMAS, axis = 2)
    psize(filtImage, 3)
    ptime(3)
    pmem(3)
    
    for i,j,k in itertools.product(range(FILTNSIGMAS), range(FOURCOUNT + 1), range(3)):
        comp = filtImage[:,:,i,j,k]
        filtImage[:,:,i,j,k] = red(np.abs(comp))
    
    weightedImage = filtImage * COMBINE
    summed = np.sum(weightedImage, axis = (-3, -2, -1))
    res = red(summed)
    np.save(IMAGENAME + "res", res)
    
for IMAGENAME in fileList:
    __main__(IMAGENAME + ".png")

