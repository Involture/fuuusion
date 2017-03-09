from glob import *
from reader import read
from fourier import fft, ifft, sf
from scipy import misc

#array conversion and type checking functions

def isul(array):
    return array.dtype == ulint
def isub(array):
    return array.dtype == ubint
def issl(array):
    return array.dtype == lint
def issb(array):
    return array.dtype == bint
def isbin(array):
    return array.dtype == np.bool

def cisul(array):
    if not isul(array):
        raise TypeError("not an unsigned little integer array")
def cisub(array):
    if not isub(array):
        raise TypeError("not an unsigned big integer array")
def cissl(array):
    if not issl(array):
        raise TypeError("not a signed little integer array")
def cissb(array):
    if not issb(array):
        raise TypeError("not a signed big integer array")
def cisbin(array):
    if not isbin(array):
        raise TypeError("not a boolean array")

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

def bToL(arr):
    shape = np.shape(arr)
    nbpix = 1
    for size in shape:
        nbpix *= size
    perpix = 100
    uplim = R
    downlim = 0
    while perpix != 98:
        lim = (uplim + downlim) // 2
        perpix = (np.sum(arr <= lim) * 100) // nbpix
        if perpix > 98:
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

def binNot(bin):
    return ulint(bin + 1 == 1)

def expend(bin, filter):
    return ulint(filterArr(bin, filter) > 0)

def erode(bin, filter):
    s = np.sum(filter)
    return ulint(filterArr(bin, filter) == s)

def match(bin, filter):
    notFilter = binNot(filter)
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

def isFilter(filterArr):
    filterSize = np.shape(filterArr)
    if (filterSize[0] % 2 != 1) or (filterSize[1] % 2 != 1):
        raise ValueError("can't filter with an even sized filter array")

def filterArr(arr, filter):
    isFilter(filter)
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
                ind = [slice(dx + i, dx + i + shape[0]), slice(dy + j, dy + j + shape[1])] + [slice(None)] * (arr.ndim - 2)
                shiftedArr = np.zeros(bigShape, dtype = bint)
                shiftedArr[ind] = filter[dx - i, dy - j] * bint(arr)
                res += shiftedArr
    ind = [slice(dx, shape[0] + dx), slice(dy, shape[1] + dy)] + [slice(None)] * (arr.ndim - 2)
    res = res[ind]
    return np.abs(res)

def filterSum(arr, filterList):
    res = np.zeros(arr.shape, dtype = ubint)
    for filter in filterList:
        cissl(filter)
        res += filterArr(arr, filter)
    return res

def filter(arr, filterList):
    cisul(arr)
    res = bint(arr.copy)
    for filter in filterList:
        if type(filter) == list:
            filterSum(res, filter)
        else:
            res = filterArr(res, filter)
    return BtoL(res)

#fourier transform filters

def squareCut(shape, percentage):
    ilim = shape[0] * percentage // 100
    jlim = shape[1] * percentage // 100
    littleShape = list(shape)
    littleShape[0] -= 2 * ilim
    littleShape[1] -= 2 * jlim
    littleShape = tuple(littleShape)
    ia1 = np.array(range(ilim, shape[0] - ilim))
    ia1 = np.array(range(jlim, shape[1] - jlim))
    ind = [slice(ilim, shape[0] - ilim), slice(jlim, shape[1] - jlim)] + [slice(None)] * (shape.ndim - 2)
    arr = np.zeros(shape, dtype = ulint)
    arr[ind] = np.ones(littleShape, dtype = ulint)
    return arr

def ellipse(shape):
    p, q = shape[:2]
    colorized = len(shape) == 3
    arr1 = np.stack([np.arange(p//2) for i in range(q)], axis = 1)
    arr1 = np.concatenate((arr1, np.flipud(arr1)), axis = 0)
    arr2 = np.stack([np.arange(q//2) for i in range(p)], axis = 0)
    arr2 = np.concatenate((arr2, np.fliplr(arr2)), axis = 1)
    arr = (arr1 * q) ** 2 + (arr2 * p) ** 2
    if colorized:
        arr = np.stack([arr for i in range(3)], axis = 2)
    return arr

def circleCut(shape, percentage):
    return 20000 * ellipse(shape) <= (shape[0] * shape[1] * percentage) ** 2

def BW(shape, n, p0):
    return (1 / (1 + (ellipse(shape) * 100 ** 2 / (shape[0] * shape[1] * p0) ** 2)))

def gauss(shape, p0):
    return np.exp(-ellipse(shape) * 100 ** 2 / (shape[0] * shape[1] * p0) ** 2)

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
