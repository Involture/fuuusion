from scipy.misc import imread
import matplotlib.pyplot as plt

from glob import *
from fourier import fft, ifft, sf
from filts import *

#image size processing

def log2(n):
    """Return the floor value of log2(n)."""
    if n <= 1:
        return 0
    else:
        return log2(n // 2) + 1

def openIm(imName, maxPow):
    """Open the file named imName as a nparray and cut it so its \ 
    dimensions are powers of two."""
    pr("opening", 1)
    im = imread(imName)
    p, q = im.shape[:2]
    r, s = log2(p), log2(q)
    pp, qq = 2 ** r, 2 ** s
    im = im[: pp, : qq]
    reduceFactor = (r + s) // 2 + (r + s) % 2 - maxPow
    im = im[::2 ** reduceFactor, ::2 ** reduceFactor]
    return im

#isp computing

def ispVect(arr):
    """Replace each pixel by an array of length 256 full of zero excepted a \
    one at the index corresponding to the value of the pixel."""
    pr("isp vectorizing", 1)
    return cint(np.stack([arr == i for i in range(256)], axis = -1))

#ploting functions

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
    return cint(filt(arr, f) > 0)

def erode(arr, f):
    s = np.sum(f)
    return cint(filt(arr, f) == s)

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
    pr("windows vectorizing", 2)
    powCheck2(p)
    powCheck2(q)
    a, b = np.shape(arr)[:2]
    la = a - p
    lb = b - q
    makeWin = lambda i, j : arr[i: i + p, j: j + q]
    makeWinGene = lambda l, j : (makeWin(i, j) for i in range(l))
    makeWinColumn = lambda j : np.stack(makeWinGene(la, j), axis = 2)
    winColumnGene = (makeWinColumn(j) for j in range(lb))
    return np.stack(winColumnGene, axis = 3)

def restoreShape(arr, dp, dq):
    a, b = arr.shape[:2]
    bigShape = (a + 2 * dp, b + 2 * dq) + arr.shape[2:]
    res = np.zeros(bigShape, dtype = cflt)
    res[dp: -dp, dq: -dq] = arr
    res[dp: -dp, : dq] = np.expand_dims(arr[:, 0], 1)
    res[dp: -dp, -dq:] = np.expand_dims(arr[:, -1], 1)
    res[: dp, dq: -dq] = np.expand_dims(arr[0, :], 0)
    res[-dp:, dq: -dq] = np.expand_dims(arr[-1, :], 0)
    res[: dp, : dq] = arr[0, 0]
    res[: dp, -dq:] = arr[0, -1]
    res[-dp: , : dq] = arr[-1, 0]
    res[-dp: , -dq:] = arr[-1, -1]
    return res

#filtering functions

def filtl(arr, fs):
    if type(fs) != np.ndarray:
        for f in fs:
            arr[...] = _filt(arr, f)
    else:
        arr[...] = _filt(arr, fs)

def filt(arr, f):
    pr("filtering", 2)
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
            pr(str(i) + str(j), 4)
            if f[dx - i, dy - j] != 0:
                shiftedArr = np.zeros(bigShape, dtype = cflt)
                shiftedArr[dx + i: dx + i + a, dy + j: dy + j + b] = f[dx - i, dy - j] * cflt(arr)
                res += shiftedArr
    return res[dx : a + dx, dy : b + dy]
    
#parabolic approximation

def parabolicApprox(xArr, yArr):
    pr("parabolic approximating", 2)
    x4sum = np.sum(xArr ** 4)
    x3sum = np.sum(xArr ** 3)
    x2sum = np.sum(xArr ** 2)
    x1sum = np.sum(xArr)
    x0sum = xArr.shape[0]
    yx2sum = np.sum(yArr * xArr ** 2, axis = -1)
    yx1sum = np.sum(yArr * xArr, axis = -1)
    yx0sum = np.sum(yArr, axis = -1)
    mat = np.array([[x4sum, x3sum, x2sum],
                    [x3sum, x2sum, x1sum],
                    [x2sum, x1sum, x0sum]])
    vect = np.stack([yx2sum, yx1sum, yx0sum], axis = -2)
    return np.linalg.solve(mat, vect)

#localisator

def stackPoints(winArr, u):
    pr("projecting points", 2)
    t = winArr.shape[0]
    c = t / 2
    circle = circleCut(t, t / 2).astype(np.bool)
    count = np.sum(circle)
    ind1 = itertools.product(range(t), repeat = 2)
    ind2 = itertools.product(range(t), repeat = 2)
    coord = lambda i, j: (scal((i, j), u), winArr[i,j])
    xGene = (scal((i, j), u) for i, j in ind1 if circle[i, j])
    yGene = (winArr[i, j] for i, j in ind2 if circle[i, j])
    xArr = np.fromiter(xGene, dtype = cflt)
    yArr = stackGene(yGene, count, winArr.ndim - 2)
    return (xArr, yArr)

def smooth(a, b, c, eps):
    pr("smoothing", 2)
    aplus = np.maximum(a, 0)
    cplus = np.maximum(c, 0)
    return -(2 * aplus * cplus) / (abs(b) + eps)

def localise(arr, t, u, eps):
    pr("localising", 2)
    wArr = winVect(arr, t, t)
    pointsCouple = stackPoints(wArr, u)
    parabArr = parabolicApprox(*pointsCouple)
    a = parabArr[..., 0]
    b = parabArr[..., 1]
    c = parabArr[..., 2]
    smoothArr = smooth(a, b, c, eps)
    return smoothArr

def localiseAndRestore(arr, t, u, eps):
    return restoreShape(localise(arr, t, u, eps), t // 2, t // 2)

#color space conversion

def RGBtoLAB(arr):
    assert(arr.dtype == cint)
    r = arr[: , : , 0].astype(cflt)
    g = arr[: , : , 1].astype(cflt)
    b = arr[: , : , 2].astype(cflt)
    lchan = (r + g + b) / 3
    achan = (g - r + 255) / 2
    bchan = (g - b + 255) / 2
    return np.stack((lchan, achan, bchan), axis = 2)
