from glob import *
from fourier import fft, ifft, sf
from scipy.misc import imread
import matplotlib.pyplot as plt

#image size processing

def log2(n):
    if n == 1:
        return 0
    else:
        return log2(n // 2) + 1

def openIm(imName):
    im = imread(imName)
    p, q = im.shape[:2]
    r, s = log2(p), log2(q)
    pp, qq = 2 ** r, 2 ** s
    im = im[:pp, :qq]
    reduceFactor = log2(r + s // 20) // 2
    im = im[::4 ** reduceFactor, ::4 ** reduceFactor]
    return im

#type converters

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
        res = (res / lim) * 255
    return ulint(res)

#isp computing

def ispVect(arr):
    return ulint(np.stack([arr == i for i in range(256)], axis = -1))

def _isp(arr):
    return np.sum(ispVect(arr), axis = (0,1), dtype = ubint)

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
    return ulint(np.sum(arr, dtype = ubint, axis = 2) // 3)

def greyMax(arr):
    return  np.max(arr, axis = 2)

def binary(arr, cutIntensity):
    return ulint(arr > cutIntensity)

#algebra operations

def binNot(arr):
    cisul(arr)
    return ulint(arr + 1 == 1)

def expend(arr, f):
    return ulint(_filt(arr, f) > 0)

def erode(arr, f):
    s = np.sum(f)
    return ulint(_filt(arr, f) == s)

def match(arr, f):
    return ulint((erode(arr, f) + erode(binNot(arr), binNot(f))) == 2)

def open(arr, f):
    return expend(erode(arr, f), f.T)

def close(arr,f):
    return erode(expend(arr, f), f.T)

def diminish(arr, f):
    return arr - match(arr, f)

def seqDiminish(arr, f):
    flist = [np.rot90(f, i) for i in range(4)]
    for rf in flist:
        arr = diminish(arr, rf)
    return arr

def fullDiminish(arr, f):
    nextbin = seqDiminish(arr, f)
    while (nextbin != arr).any():
        arr = nextbin
        nextbin = seqDiminish(arr,f)
    return arr

def squeletize(arr):
    nextarr = seqDiminish(seqDiminish(arr, filters["flat"]), filters["lcorner"])
    while nextarr != arr:
        arr = nextarr
        nextarr = seqDiminish(seqDiminish(nextbin, filters["flat"]), filters["lcorner"])
    return arr

#window vectorization

def winVect(arr, p, q):
    powCheck(p)
    powCheck(q)
    cisul(arr)
    a, b = np.shape(arr)[:2]
    la = a - p
    lb = b - q
    return np.stack([np.stack([arr[i : i + p, j : j + q] for i in range(la)], axis = 2) for j in range(lb)], axis = 3)

#filtering functions

def ispFilt(arr, f):
    fpos = np.maximum(f, 0)
    fneg = np.minimum(f, 0)
    return ubint(np.sqrt(np.sum((_filt(arr, fpos)  +  _filt(arr, fneg)) ** 2, axis = -1, dtype = ubint)))

def filt(arr, fs):
    if type(fs) == list:
        for f in fs:
            arr = _filt(arr, f)
    else:
        arr = _filt(arr, fs)

def _filt(arr, f):
    cisfilt(f)
    shape = arr.shape
    filtShape = f.shape

    dx = (filtShape[0] - 1) // 2
    dy = (filtShape[1] - 1) // 2
    bigShape = list(shape)
    bigShape[0] += 2 * dx
    bigShape[1] += 2 * dy
    bigShape = tuple(bigShape)
    res = np.full(bigShape, 0, dtype = bint)
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            print(str(i) + str(j))
            if f[dx - i, dy - j] != 0:
                shiftedArr = np.zeros(bigShape, dtype = bint)
                shiftedArr[dx + i : dx + i + shape[0], dy + j : dy + j + shape[1]] = f[dx - i, dy - j] * bint(arr)
                res += shiftedArr
    res = res[dx : shape[0] + dx, dy : shape[1] + dy]
    return res

#color space conversion

def RGBtoLAB(arr):
    arr =  np.stack([ulint(np.sum(ubint(arr), axis = 2) // 3), arr[:,:,1] - arr[:,:,0], arr[:,:,1] - arr[:,:,2]], axis = 2)

def LABtoRGB(arr):
    G = np.sum(arr, axis = 2)
    arr =  ulint(np.stack([G - arr[:,:,1], G, G - arr[:,:,2]], axis = 2))
