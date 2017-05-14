from glob import *

#vector of exponantial

def exp(n, l, ax, inverse = False):
    if inverse:
        eps = 1
    else:
        eps = -1
    expRow = np.exp((eps * 2 * 1j * np.pi / n) * np.fromiter([i for i in range(n)], dtype = np.complex64))
    exp = np.expand_dims(expRow, axis = ax)
    exp = expand(exp, l)
    return exp

#efficient fourier transform

def hVectFastFourier2D(arr, inverse = False):
    cispow(arr)
    shape = np.shape(arr)
    el = len(shape) - 2
    n = shape[1]
    m = n // 2
    e = exp(n, el, 0, inverse)
    j = n // 2
    while j >= 1:
        for k in range(j):
            print(str(j) + " " + str(k))
            temp = arr[:, k : k + j + n : j].copy()
            arr[:, k : k + m : j] = arr[:, k : k + n : 2 * j] + arr[:, k + j : k + j + n : 2 * j] * e[:, 0 : n // 2 : j]
            arr[:, k + m : k + n : j] = temp[:, : -1 : 2] + temp[:, 1 :: 2] * e[:, n // 2 : n : j]
            ptime()
            pmem()
        j = j // 2
    arr[...] =  arr / np.sqrt(n)

def vVectFastFourier2D(arr, inverse = False):
    cispow(arr)
    shape = np.shape(arr)
    el = len(shape) - 2
    n = shape[0]
    m = n // 2
    e = exp(n, el, 1, inverse)
    j = n // 2
    while j >= 1:
        for k in range(j):
            print(str(j) + " " + str(k))
            temp = arr[k : k + j + n : j, :].copy()
            arr[k : k + m : j, :] = arr[k : k + n : 2 * j, :] + arr[k + j : k + j + n : 2 * j, :] * e[0 : n // 2 : j, :]
            arr[k + m : k + n : j, :] = temp[: -1 : 2, :] + temp[1 :: 2, :] * e[n // 2 : n : j, :]
            ptime()
            pmem()
        j = j // 2
    arr[...] =  arr / np.sqrt(n)

def fastFastFourier2D(arr, inverse = False):
    res = arr.astype(np.complex64)
    hVectFastFourier2D(res, inverse)
    vVectFastFourier2D(res, inverse)
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

