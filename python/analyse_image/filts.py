from glob import *

#standard filters

def blur(d):
    t = 2 * d + 1
    return np.ones((t,t), dtype = lint)

def dgrad(d):
    t = 2 * d + 1
    onesArr = np.ones((1,t), dtype = lint)
    row = np.array(np.concatenate((onesArr, -onesArr), axis = 1))
    row[0][0] = 0
    res = np.array(row)
    for j in range (t - 1):
        permute(row)
        res = np.concatenate((res,row), axis = 0)
    return res[:, :7]

def vgrad(d):
    t = 2 * d + 1
    return np.concatenate((-np.ones((t, d), dtype = lint), np.zeros((t, 1), dtype = lint), np.ones((t,d), dtype = lint)), axis = 1)

def pdgrad(d):
    t = 2 * d + 1
    row = np.array(np.concatenate((np.arange(t), -np.arange(t - 1, 0, -1)), axis = 0), ndmin = 2)
    res = np.array(row)
    for j in range(t - 1):
        permute(row)
        res = np.concatenate((res, row), axis = 0)
    return res[:,:t]

def pvgrad(d):
    t = 2 * d + 1
    arr1 = np.stack([np.arange(1, d + 1) for i in range(t)], axis = 1)
    arr2 = np.stack([-np.arange(d + 1) for i in range(t)], axis = 1)
    return np.concatenate((np.flipud(arr1), arr2), axis = 0)

#topological elements

lcorner = np.array([[0,0,0],[1,1,0],[1,1,0]], dtype = ulint)
bcorner = np.array([[1,0,0],[1,1,0],[1,1,1]], dtype = ulint)

#fourier transform filters functions

def vdist(p, q):
    arr1 = np.stack([np.arange(p // 2) for i in range(q)], axis = 1)
    arr2 = np.stack([np.arange(p // 2 + (p % 2)) for i in range(q)], axis = 1)
    return np.concatenate((arr1, np.flipud(arr2)), axis = 0)

def permute(arr):
    row = arr[0]
    temp = row[-1]
    arr[0][1:] = row[:-1]
    arr[0][0] = temp

def ddist(t):
    row = np.array(np.concatenate((np.arange(t), np.arange(t - 1, 0, -1)), axis = 0), ndmin = 2)
    res = np.array(row)
    for j in range(t - 1):
        permute(row)
        res = np.concatenate((res, row), axis = 0)
    return res[:,:t]

def ellipse(p, q):
    return (vdist((p, q)) * q) ** 2 + (vdist((q, p)).T * p) ** 2

def squareCut(p, q, p0):
    ilim = p * p0 // 200
    jlim = p * p0 // 200
    arr = np.zeros((p, q), dtype = ulint)
    arr[ilim : -ilim, jlim : -jlim] = np.ones((p - 2 * ilim, q - 2 * jlim), dtype = ulint)
    return arr

def circleCut(p, q, p0):
    return 20000 * ellipse(p, q) <= (p * q * p0) ** 2

def BW(p, q, p0):
    return (1 / (1 + (ellipse(p, q) * 160000 / (p * q * p0) ** 2)))

def cgauss(p, q, p0):
    return np.exp(-ellipse(p, q) * 160000 / (p * q * p0) ** 2)

def vgauss(t, p0):
    return np.exp(-vdist(t, t) ** 2 * 40000 / (t * p0)  ** 2)

def hgauss(t, p0):
    return vgauss(t, p0).T

def d1gauss(t, p0):
    return np.exp(-ddist(t) ** 2 * 40000 / (t * p0) ** 2)

def d2gauss(t, p0):
    return np.flipud(d1gauss(t, p0))
