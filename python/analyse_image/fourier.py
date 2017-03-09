from glob import *

#vector of exponantial

def exp(n, inverse = False):
    if inverse:
        eps = 1
    else:
        eps = -1
    return np.exp((eps * 2 * 1j * np.pi / n) * np.fromiter([i for i in range(n)], dtype = np.complex))

#inefficient fourier transform

def fourier1D(arr1D, inverse = False):
    n = np.shape(arr1D)[0]
    e = exp(n, inverse)
    res = np.zeros((n,), dtype = np.complex)
    for i in range(n):
        eTab = np.array([e[(j * i) % n] for j in range(n)])
        res[i] = np.sum(arr1D * eTab)
    return res / np.sqrt(n)

def fastFourier1D(arr1D, inverse = False):
    n = np.shape(arr1D)[0]
    m = n // 2
    e = exp(n, inverse)
    res = np.copy(arr1D)
    temp = np.zeros((n,), dtype = np.complex)
    j = n // 2
    while j >= 1:
        for k in range(j):
            temp[k : k + m : j] = res[k : k + n : 2 * j] + res[k + j : k + j + n : 2 * j] * e[0 : n // 2 : j]
            temp[k + m : k + n : j] = res[k : k + n : 2 * j] + res[k + j : k + j + n : 2 * j] * e[(n // 2) : n : j]
        j = j // 2
        res = np.copy(temp)
    return res / np.sqrt(n)

def fastFourier2D(arr2D, inverse = False):
    shape = np.shape(arr2D)
    fourierH = np.row_stack([fastFourier1D(arr2D[i,:], inverse) for i in range(shape[0])])
    fourierHV = np.column_stack([fastFourier1D(fourierH[:,j], inverse) for j in range(shape[1])])
    return fourierHV

#efficient fourier transform

def hVectFastFourier2D(arr2D, inverse = False):
    shape = np.shape(arr2D)
    colorized = (len(shape) == 3)
    n = shape[1]
    m = n // 2
    e1D = exp(n, inverse)
    e = np.row_stack([e1D for i in range(shape[0])])
    if colorized:
        e = np.stack([e for i in range(3)], axis=2)
    res = np.copy(arr2D)
    temp = np.zeros((shape), dtype = np.complex)
    j = n // 2
    if colorized:
        while j >= 1:
            for k in range(j):
                temp[:, k : k + m : j,:] = res[:, k : k + n : 2 * j,:] + res[:, k + j : k + j + n : 2 * j,:] * e[:, 0 : n // 2 : j,:]
                temp[:, k + m : k + n : j,:] = res[:, k : k + n : 2 * j,:] + res[:, k + j : k + j + n : 2 * j,:] * e[:, (n // 2) : n : j,:]
            j = j // 2
            res = np.copy(temp)
    else:
        while j >= 1:
            for k in range(j):
                temp[:, k : k + m : j] = res[:, k : k + n : 2 * j] + res[:, k + j : k + j + n : 2 * j] * e[:, 0 : n // 2 : j]
                temp[:, k + m : k + n : j] = res[:, k : k + n : 2 * j] + res[:, k + j : k + j + n : 2 * j] * e[:, (n // 2) : n : j]
            j = j // 2
            res = np.copy(temp)
    return res / np.sqrt(n)

def fastFastFourier2D(arr2D, inverse = False):
    res = hVectFastFourier2D(arr2D, inverse)
    res = np.rot90(hVectFastFourier2D(np.rot90(res), inverse), 3)
    return res

fft = fastFastFourier2D
ifft = lambda x : fastFastFourier2D(x, True)

#function to render fourier transform

def flip(arr2D):
    return arr2D[::-1, ::-1]

def switchQuad(arr2D):
    m1 = np.shape(arr2D)[0] // 2
    m2 = np.shape(arr2D)[1] // 2
    res = np.copy(arr2D)
    res[:m1, :m2] = flip(res[:m1, :m2])
    res[:m1, m2:] = flip(res[:m1, m2:])
    res[m1:, :m2] = flip(res[m1:, :m2])
    res[m1:, m2:] = flip(res[m1:, m2:])
    return res

def ubintFromComplex(complexArr):
    return ubint(np.abs(complexArr))

def reduceFourier(arr2D, cut):
    return ulint(np.minimum(ubintFromComplex(arr2D), cut) * 255 // cut)

def ulintFourier(arr2D, cut):
    return switchQuad(reduceFourier(arr2D, cut))

sf = ulintFourier
