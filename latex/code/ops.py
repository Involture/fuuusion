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
    r, g, b = im[:,:,0], im[:,:,1], im[:,:,2]
    mr, mg, mb = r.max(), g.max(), b.max()
    nr, ng, nb = r / mr, g / mg, b / mb
    return np.stack((nr, ng, nb), axis = 2)

#isp computing

def ispVect(arr):
    """Replace each pixel by an array of length 256 full of zero excepted a \
    one at the index corresponding to the value of the pixel."""
    pr("isp vectorizing", 1)
    return cflt(np.stack([arr == i for i in range(256)], axis = -1))

#ploting functions

def show(arr):
    """Display the image."""
    fig = plt.figure()
    subp = fig.add_subplot(111)
    subp.imshow(arr)
    subp.axis('off')
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
    fig.show()

#grey level and binary

def greyAv(arr):
    """Return a grey level image from a coloured image taking the average \
    of the three color channel."""
    return cflt(np.sum(arr, dtype = cflt, axis = 2) // 3)

def greyMax(arr):
    """Return a grey level image from a coloured image taking the maximum \
    of the three color channel."""
    return  np.max(arr, axis = 2)

def binary(arr, cutIntensity):
    """Return a binary image from a grey level image where only pixels with \
    an intensity superior to cutIntensity are set to 1."""
    return cflt(arr > cutIntensity)

#window vectorization

def winVect(arr, p, q):
    pr("windows vectorizing", 2)
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
    """Filter arr with f on his first two dimensions."""
    shape = arr.shape
    a, b = shape[: 2]
    filtShape = f.shape

    #Create a bigger shape to fit the shifted images
    dx = (filtShape[0] - 1) // 2
    dy = (filtShape[1] - 1) // 2
    bigShape = list(shape)
    bigShape[0] += 2 * dx
    bigShape[1] += 2 * dy
    bigShape = tuple(bigShape)
    
    res = np.zeros(bigShape, dtype = cflt)
    #iteration over the filter
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            #checking if the filter value is non-zero to avoid useless calculus
            if f[dx - i, dy - j] != 0:
                #creating the shifted image
                shiftedArr = np.zeros(bigShape, dtype = cflt)
                shiftedArr[dx + i: dx + i + a, dy + j: dy + j + b] = f[dx - i, dy - j] * cflt(arr)
                #adding it to the result
                res += shiftedArr
    return res[dx : a + dx, dy : b + dy]
    
#parabolic approximation

def parabolicApprox(xArr, yArr):
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
    t = winArr.shape[0]
    c = t / 2 - 0.5
    circle = circleCut(t, t / 2).astype(np.bool)
    count = np.sum(circle)
    #Product returns the cartesian product of two geneators, usefull to avoid nested loops.
    ind1 = itertools.product(range(t), repeat = 2)
    ind2 = itertools.product(range(t), repeat = 2)
    xGene = (scal((i - c, j - c), u) for i, j in ind1 if circle[i, j])
    yGene = (winArr[i, j] for i, j in ind2 if circle[i, j])
    #The vector of coordinates in the image plane.
    xArr = np.fromiter(xGene, dtype = cflt)
    #The vector of intensity of the pixels.
    yArr = stackGene(yGene, count, winArr.ndim - 2)
    return (xArr, yArr)

def smooth(a, b, c, eps):
    aminus = np.maximum(-a, 0)
    cplus = np.maximum(c, 0)
    return (2 * aminus * cplus) / (np.abs(b) + eps)

def localise(arr, t, u, eps):
    wArr = winVect(arr, t, t)
    pointsCouple = stackPoints(wArr, u)
    parabArr = parabolicApprox(*pointsCouple)
    a = parabArr[..., 0, :]
    b = parabArr[..., 1, :]
    c = parabArr[..., 2, :]
    smoothArr = smooth(a, b, c, eps)
    return smoothArr

def localiseAndRestore(arr, t, u, eps):
    return restoreShape(localise(arr, t, u, eps), t // 2, t // 2)

#color space conversion

def RGBtoLAB(arr):
    r = arr[: , : , 0]
    g = arr[: , : , 1]
    b = arr[: , : , 2]
    l = (r + g + b) / 3
    a = (g - r + 1.) / 2
    b = (g - b + 1.) / 2
    return np.stack((l, a, b), axis = 2)

#algebra operations

def binNot(arr):
    """The not operator on a binary image."""
    return cflt(arr + 1. == 1.)

def expand(arr, f):
    return cflt(filt(arr, f) > 0)

def erode(arr, f):
    s = np.sum(f)
    return cflt(filt(arr, f) == s)

def match(arr, f):
    fpos = np.maximum(f, 0.)
    fneg = np.minimum(f, 0.)
    fneg = np.abs(fneg)
    ps = np.sum(fpos)
    present = cflt(filt(arr, fpos) == ps)
    absent = cflt(filt(arr, fneg) == 0.)
    return cflt(present + absent == 2)

def open(arr, f):
    return expand(erode(arr, f), f)

def close(arr,f):
    return erode(expand(arr, f), f)

def seqClose(arr, t, ndims):
    for theta in np.linspace(0., np.pi * 2, ndims, endpoint = False):
        arr = close(arr, line(t, theta))
    return arr

#polygone creation

south = [1, 0]
north = [-1, 0]
east = [0, 1]
west = [0, -1]
convert = [south, east, north, west]

def safe(bim, cell):
    x, y = cell
    p, q = bim.shape
    return (x < p and y < q and bim[x, y])

def next(bim, pos, head):
    indTryOrder = [(i - 1) % 4 for i in range(head, head + 4)]
    tryOrder = [(i, pos + convert[i]) for i in indTryOrder]
    for i, cell in tryOrder:
        if safe(bim, cell):
            return i, cell
    raise Exception

def run(bim):
    xs, ys = np.nonzero(bim)
    x = np.min(xs)
    i = np.nonzero(xs == x)[0][0]
    y = ys[i]
    p0 = np.array([x, y])
    head = 0
    head, p = next(bim, p0, head)
    res = np.stack((p0, p))
    while (p != p0).any():
        head, p = next(bim, p, head)
        res = np.append(res, np.expand_dims(p, axis = 0), axis = 0)
    return res

def segdist(shape, point, u):
    x, y = point
    a, b = u / norm(u)
    p, q = shape
    ind = np.mgrid[:p, :q]
    ind = np.transpose(ind, (1,2,0))
    ind = ind - [x, y]
    scal = ind[:,:,0] * b - ind[:,:,1] * a
    return scal

def further(shape, contour, start, stop, threshold):
    p1 = contour[start]
    p2 = contour[stop]
    u = p2 - p1
    d = segdist(shape, p1, u)
    d = np.abs(d)
    contourSlice = contour[start + 1:stop]
    dist = d[contourSlice[:,0], contourSlice[:,1]]
    m = dist.max()
    if m > threshold:
        ind = np.nonzero(dist == m)
        return ind[0][0] + start + 1
    else:
        return -1

def updateList(shape, contour, indList, threshold):
    res = []
    it = zip(indList[:-1], indList[1:])
    changed = False
    for i1, i2 in it:
        res.append(i1)
        i = further(shape, contour, i1, i2, threshold)
        if i != -1:
            res.append(i)
            changed = True
    res.append(len(contour) - 1)
    return res, changed

def initSeg(contour):
    return [0, len(contour) // 2, len(contour) - 1]

def polygonize(shape, contour, threshold):
    l = initSeg(contour)
    notFinished = True
    while notFinished:
        l, notFinished = updateList(shape, contour, l, threshold)
    return contour[l[:-1]]

def showList(shape, contour, l):
    p, q = shape
    arr = np.zeros((p + 10, q + 10))
    for i in l:
        point = contour[i]
        arr[point[0]: point[0] + 10, point[1]: point[1] + 10] = 1.
    return arr
