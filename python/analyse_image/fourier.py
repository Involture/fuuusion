from glob import *

#vector of exponantial

def exp(n, inverse = False):
    """Generate a vector of exp(2i*pi n) if inverse is False of \
    exp(-2i*pi/n) either."""
    if inverse:
        eps = 1
    else:
        eps = -1
    expRow = np.exp((eps * 2 * 1j * np.pi / n) * np.arange(n, dtype = ccpl))
    return expRow

#efficient fourier transform

def hVectFastFourier2D(arr, inverse = False):
    """Compute the fourier transform of arr on the second dimension \
    IN PLACE."""
    cispow(arr)
    shape = np.shape(arr)
    el = len(shape) - 2
    n = shape[1]
    m = n // 2
    e = exp(n, inverse)
    j = n // 2
    while j >= 1:
        for k in range(j):
            p(str(j) + " " + str(k))
            temp = arr[:, k: k + j + n: j].copy()
            a = temp[:, : -1: 2]
            b = temp[:, 1:: 2]
            b1 = multDim(b, e[0: n // 2: j], [1])
            b2 = multDim(b, e[n // 2: n: j], [1])
            arr[:, k: k + m: j] = a + b1
            arr[:, k + m: k + n: j] = a + b2
            pptime()
            ppmem()
        j = j // 2
    arr[...] =  arr / np.sqrt(n)

def vVectFastFourier2D(arr, inverse = False):
    """Compute the fourier transform of arr on the first dimension IN PLACE."""
    cispow(arr)
    shape = np.shape(arr)
    el = len(shape) - 2
    n = shape[0]
    m = n // 2
    e = exp(n, inverse)
    j = n // 2
    while j >= 1:
        for k in range(j):
            p(str(j) + " " + str(k))
            temp = arr[k: k + j + n: j,:].copy()
            a = temp[: -1: 2, :]
            b = temp[1:: 2, : ]
            b1 = multDim(b, e[0: n // 2: j], [0])
            b2 = multDim(b, e[n // 2: n: j], [0])
            arr[k: k + m: j,:] = a + b1
            arr[k + m: k + n: j,:] = a + b2
            pptime()
            ppmem()
        j = j // 2
    arr[...] =  arr / np.sqrt(n)

def fastFastFourier2D(arr, inverse = False):
    """Compute the fourier transform of arr on the first \
    two dimensions IN PLACE."""
    res = arr.astype(np.complex64)
    hVectFastFourier2D(res, inverse)
    vVectFastFourier2D(res, inverse)
    return res

fft = fastFastFourier2D
ifft = lambda x: fastFastFourier2D(x, True)

#function to render fourier transform

def showFourier(arr2D, amp = 1):
    """Switch quadrant and convert arr2D in integer.\
    Makes a fourier transform result displayable."""
    amplified = np.minimum(arr2D * 255 * amp / np.max(arr2D), 255.)
    return switchQuad(cint(amplified))

sf = showFourier

#inefficient fourier transform

def fourier1D(arr1D, inverse = False):
    n = np.shape(arr1D)[0]
    e = exp(n, inverse)
    res = np.zeros((n,), dtype = ccpl)
    for i in range(n):
        eTab = np.array([e[(j * i) % n] for j in range(n)])
        res[i] = np.sum(arr1D * eTab)
    return res / np.sqrt(n)

def fastFourier1D(arr1D, inverse = False):
    n = np.shape(arr1D)[0]
    m = n // 2
    e = exp(n, inverse)
    res = np.copy(arr1D)
    temp = np.zeros((n,), dtype = ccpl)
    j = n // 2
    while j >= 1:
        for k in range(j):
            temp[k: k + m: j] = res[k: k + n: 2 * j] + res[k + j: k + j + n: 2 * j] * e[0: n // 2: j]
            temp[k + m: k + n: j] = res[k: k + n: 2 * j] + res[k + j: k + j + n: 2 * j] * e[(n // 2): n: j]
        j = j // 2
        res = np.copy(temp)
    return res / np.sqrt(n)

def fastFourier2D(arr2D, inverse = False):
    shape = np.shape(arr2D)
    fourierH = np.row_stack([fastFourier1D(arr2D[i,:], inverse) for i in range(shape[0])])
    fourierHV = np.column_stack([fastFourier1D(fourierH[:,j], inverse) for j in range(shape[1])])
    return fourierHV


