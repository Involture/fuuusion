from glob import *
from reader import read
from fourier import fft, ifft, sf, mfft
from scipy import misc

#type converters

def _2Disp(arr, size):
    isp = np.zeros(N, dtype = ubint)
    for i in range(N):
        isp[i] = np.sum(arr == i)
    return isp

def _isp(arr):
    shape = np.shape(arr)
    if len(shape) >= 3:
        size = shape[:2]
        return (_2Disp(arr[:,:,0], size), _2Disp(arr[:,:,1], size), _2Disp(arr[:,:,0], size))
    else:
        size = shape
        return _2Disp(arr, size)

def bToL(arr, loss):
    shape = np.shape(arr)
    nbpix = 1
    for size in shape:
        nbpix *= size
    perpix = 100
    uplim = R
    downlim = 0
    while perpix != loss:
        lim = (uplim + downlim) // 2
        perpix = (np.sum(arr <= lim) * 100) // nbpix
        if perpix > loss:
            uplim = lim
        else:
            downlim = lim
    lim = (uplim + downlim) // 2
    res = arr.copy()
    res = np.minimum(res, lim)
    if np.max(arr) * 255 < R:
        res = (res * 255) // lim
    else:
        res = (res // lim) * 255
    return ulint(res)

red = lambda x : bToL(x, 98)

#ploting functions

def _plotSpec(subp, spectrum, c):
    """Plot the intensity spectrum to the specified subplot
    c is the color of the plot
    """
    ind = np.arange(N)
    width = 1

    subp.bar(ind, spectrum, width, color = c, edgecolor = c)

    subp.set_xlabel("pixel value")
    subp.set_ylabel("nb of pixels")

def plotIsp(arr):
    isp = _isp(arr)
    n = len(isp)
    fig = plt.figure()
    if n != 3:
        sub = fig.add_subplot(111)
        _plotSpec(sub, isp, "black")
    else:
        sub1 = fig.add_subplot(131)
        sub2 = fig.add_subplot(132)
        sub3 = fig.add_subplot(133)
        _plotSpec(sub1, isp[0], "red")
        _plotSpec(sub2, isp[1], "green")
        _plotSpec(sub3, isp[2], "blue")
    fig.show()

def show(arr):
    """Display the image."""
    size = arr.shape
    array = arr.copy()
    if np.max(arr) == 1:
        array = array * 255
    fig = plt.figure()
    subp = fig.add_subplot(111)
    subp.imshow(array)
    subp.axis('off')
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
    fig.show()

#grey level and binary

def greyAv(arr):
    return ulint(np.sum(ubint(arr), dtype = ubint, axis = 2) // 3)

def greyMax(arr):
    return  np.max(arr, axis = 2)

def binary(arr, cutIntensity):
    return ulint(arr > cutIntensity)

#algebra operations

def binNot(arr):
    cisul(arr)
    return ulint(bin + 1 == 1)

def expend(bin, filter):
    return ulint(filterArr(bin, filter) > 0)

def erode(bin, filter):
    s = np.sum(filter)
    return ulint(filterArr(bin, filter) == s)

def match(bin, filter):
    return ulint((erode(bin, filter) + erode(binNot(bin), binNot(filter))) == 2)

def open(bin, filter):
    return expend(erode(bin, filter), filter.T)

def close(bin,filter):
    return erode(expend(bin, filter), filter.T)

def diminish(bin, filter):
    return bin - match(bin, filter)

def seqDiminish(bin, filter):
    filterList = [np.rot90(filter, i) for i in range(4)]
    for f in filterList:
        bin = diminish(bin, f).copy()
    return bin

def fullDiminish(bin, filter):
    nextbin = seqDiminish(bin, filter)
    while (nextbin != bin).any():
        bin = nextbin
        nextbin = seqDiminish(bin,filter)
    return bin

def squeletize(bin):
    nextbin = seqDiminish(seqDiminish(bin, filters["flat"][0]), filters["lcorner"][0])
    while nextbin != bin:
        bin = nextbin
        nextbin = seqDiminish(seqDiminish(nextbin, filters["flat"][0]), filters["lcorner"][0])
    return bin 

#filtering functions

def filter(arr, filter):
    cisfilt(filter)
    shape = arr.shape
    filterShape = filter.shape

    dx = (filterShape[0] - 1) // 2
    dy = (filterShape[1] - 1) // 2
    bigShape = list(shape)
    bigShape[0] += 2 * dx
    bigShape[1] += 2 * dy
    bigShape = tuple(bigShape)
    res = np.full(bigShape, 0, dtype = bint)
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            if filter[dx - i, dy - j] != 0:
                #ind = [slice(dx + i, dx + i + shape[0]), slice(dy + j, dy + j + shape[1])] + [slice(None)] * (arr.ndim - 2)
                shiftedArr = np.zeros(bigShape, dtype = bint)
                shiftedArr[dx + i : dx + i + shape[0], dy + j : dy + j + shape[1]] = filter[dx - i, dy - j] * bint(arr)
                res += shiftedArr
    #ind = [slice(dx, shape[0] + dx), slice(dy, shape[1] + dy)] + [slice(None)] * (arr.ndim - 2)
    res = res[dx : shape[0] + dx, dy : shape[1] + dy]
    return np.abs(res)

f = lambda x, y : red(filter(x, y))

#color space conversion

def RGBtoLAB(arr):
    return np.stack([ulint(np.sum(ubint(arr), axis = 2) // 3), arr[:,:,1] - arr[:,:,0], arr[:,:,1] - arr[:,:,2]], axis = 2)

def LABtoRGB(arr):
    G = np.sum(arr, axis = 2)
    return ulint(np.stack([G - arr[:,:,1], G, G - arr[:,:,2]], axis = 2))

camille = misc.imread("camille.png")[:1024, :1024,:]
im1 = misc.imread("image1.png")[:512, :1024,:]
bin = binary(greyAv(im1),100)
filters = read("filter")
