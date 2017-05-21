from glob import *
from fourier import fft, ifft, sf
from scipy.misc import imread
import matplotlib.pyplot as plt

#image size processing

def log2(n):
    """Return the floor value of log2(n)."""
    if n == 1:
        return 0
    else:
        return log2(n // 2) + 1

def openIm(imName):
    """Open the file named imName as a nparray and cut it so its \ 
    dimensions are powers of two."""
    im = imread(imName)
    p, q = im.shape[:2]
    r, s = log2(p), log2(q)
    pp, qq = 2 ** r, 2 ** s
    im = im[:pp, :qq]
    reduceFactor = log2(r + s // 20) // 2
    im = im[::4 ** reduceFactor, ::4 ** reduceFactor]
    return im

#isp computing

def ispVect(arr):
    """Replace each pixel by an array of length 256 full of zero excepted a \
    one at the index corresponding to the value of the pixel."""
    return ulint(np.stack([arr == i for i in range(256)], axis = -1))

def _isp(arr):
    return np.sum(ispVect(arr), axis = (0,1), dtype = ubint)

#ploting functions

def _plotSpec(subp, spectrum, c):
    """Plot the intensity spectrum to the specified subplot
    c is the color of the plot."""
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
    """Return a grey level image from a coloured image taking the average \
    of the three color channel."""
    return cint(np.sum(arr, dtype = cflt, axis = 2) // 3)

def greyMax(arr):
    """Return a grey level image from a coloured image taking the maximum \
    of the three color channel."""
    return  np.max(arr, axis = 2)

def binary(arr, cutIntensity):
    """Return a binary image from a grey level image where only pixels with \
    an intensity superior to cutIntensity are set to 1."""
    return cint(arr > cutIntensity)

#algebra operations

def binNot(arr):
    """The not operator on a binary image."""
    return cint(arr + 1 == 1)

def expand(arr, f):
    return cint(_filt(arr, f) > 0)

def erode(arr, f):
    s = np.sum(f)
    return cint(_filt(arr, f) == s)

def match(arr, f):
    return cint((erode(arr, f) + erode(binNot(arr), binNot(f))) == 2)

def open(arr, f):
    return expand(erode(arr, f), f.T)

def close(arr,f):
    return erode(expand(arr, f), f.T)

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

#window vectorization

def winVect(arr, p, q):
    powCheck(p)
    powCheck(q)
    a, b = np.shape(arr)[:2]
    la = a - p
    lb = b - q
    makeWin = lambda i, j : arr[i: i + p, j: j + q]
    makeWinGene = lambda l, j : (makeWin(i, j) for i in range(l))
    makeWinColumn = lambda j : np.stack(makeWinGene(la, j), axis = -1)
    winColumnGene = (makeWinColumn for j in range(lb))
    return np.stack(winColumnGene, axis = -1)

#filtering functions

def filt(arr, fs):
    if type(fs) == list:
        for f in fs:
            arr[...] = _filt(arr, f)
    else:
        arr[...] = _filt(arr, fs)

def _filt(arr, f):
    """Filter arr with f on his first two dimensions."""
    cisfilt(f)
    shape = arr.shape
    a, b = shape[: 2]
    filtShape = f.shape

    dx = (filtShape[0] - 1) // 2
    dy = (filtShape[1] - 1) // 2
    bigShape = list(shape)
    bigShape[0] += 2 * dx
    bigShape[1] += 2 * dy
    bigShape = tuple(bigShape)

    res = np.zeros(bigShape, dtype = cflt)
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            pp(str(i) + str(j))
            if f[dx - i, dy - j] != 0:
                shiftedArr = np.zeros(bigShape, dtype = cflt)
                shiftedArr[dx + i: dx + i + a, dy + j: dy + j + b] = f[dx - i, dy - j] * cflt(arr)
                res += shiftedArr
    return res[dx : a + dx, dy : b + dy]
    

#color space conversion

def RGBtoLAB(arr):
    assert(arr.dtype == cint)
    lchanFloat = np.sum(arr, axis = 2, dtype = cflt) // 3
    lchan = cint(lchanFloat)
    achan = arr[: , : , 1] - arr[: , : , 0]
    bchan = arr[: , : , 1] - arr[: , : , 2]
    arr =  np.stack(lchan, achan, bchan, axis = 2)

def LABtoRGB(arr):
    assert(arr.dtype == cint)
    g = np.sum(arr, axis = 2)
    r = g - arr[: , : , 1]
    b = g - arr[: , : , 2]
    arr =  ulint(np.stack(r, g, b, axis = 2))
