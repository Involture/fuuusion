print("\nimporting modules")

from ops import *

fileList = ["roller2.png", "bureau.png", "pq1.png"]

def __main__(IMAGENAME):
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
    
    #actual script
    
    #>>> Opening image
    image = openIm("test_images/" + IMAGENAME, MAXPOW)
    #show(image)
    psize(image, 3)
    ptime(3)
    pmem(3)
    
    #>>> To LAB!
    labImage = RGBtoLAB(image)
    ptime(3)
    pmem(3)
    np.save(IMAGENAME + "lab", labImage)
    
    #>>> FourierTransforming
    fftImage = fft(labImage)
    np.save(IMAGENAME + "fft", fftImage)
    
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
    combinedResp = stackGene(respGene, FOURCOUNT, 2)
    np.save(IMAGENAME + "resp", combinedResp[:,:,0,0])
    ptime(3)
    pmem(3)
    
    #>>> Concatenate with original image
    pr("\nconcatenating with original image", 2)
    restoredCombinedResp = restoreShape(combinedResp, FOURWINSIZE // 2, FOURWINSIZE // 2)
    print(restoredCombinedResp.shape)
    print(labImage.shape)
    print(np.expand_dims(labImage, axis = 0).shape)
    combinedImage = np.append(restoredCombinedResp, np.expand_dims(labImage, axis =2),  axis = 2)
    psize(combinedImage, 3)
    ptime(3)
    pmem(3)
    
    #>>>filtering
    pr("\nFiltering and localising all channels", 2)
    c = circleCut(FILTSIZE, FILTSIZE / 2)
    
    def detect(u, sigma):
        f = gaussd(FILTSIZE, u, sigma) * c
        filtered = filt(combinedImage, f)
        localised = localise(filtered, LOCALISESIZE, u, EPS)
        if u == (1., 0.) and sigma == 1.:
            np.save(IMAGENAME + "filtered", filtered[:,:,0,0])
            np.save(IMAGENAME + "localised", localised[:,:,0,0])
        return localised
    
    def summedDetect(sigma):
        return np.sum((detect(u, sigma) for u in dirGene(FILTNDIR)), axis = 2)
    
    filtImageGene = (summedDetect(sigma) for sigma in FILTSIGMAS)
    filtImage = stackGene(filtImageGene, FILTNSIGMAS, axis = 2)
    psize(filtImage, 3)
    ptime(3)
    pmem(3)
    
    #>>> Reducing
    for i,j,k in itertools.product(range(FILTNSIGMAS), range(FOURCOUNT + 1), range(3)):
        comp = filtImage[:,:,i,j,k]
        filtImage[:,:,i,j,k] = red(np.abs(comp))
    
    weightedImage = filtImage * COMBINE
    summed = np.sum(weightedImage, axis = (-3, -2, -1))
    res = red(summed)
    np.save(IMAGENAME + "res", res)
    
for IMAGENAME in fileList:
    __main__(IMAGENAME)

