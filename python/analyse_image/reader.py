from glob import *

def charListToVal(clist):
    i = 0
    res = []
    while i < len(clist):
        if clist[i] == "-":
            res.append(-int(clist[i+1]))
            i += 2
        else:
            res.append(int(clist[i]))
            i += 1
    return res

def toArr(lines):
    rows = []
    for line in lines:
        valList = charListToVal(list(line))
        rows.append(np.fromiter(valList, dtype = lint))
    return np.row_stack(rows)

def fromSingle(string):
    split = string.split("\n")
    return (split[0], toArr(split[1:]))

def fromGroup(dic, l):
    split = l.split("\n")
    for i in range(1, len(split)):
        split[i] = dic[split[i]]
    return (split[0], split[1:])

def read(dbaseName):
    file = open(dbaseName + ".dbase", "r")
    string0 = file.read()[:-1]
    funcs = string0.split("\n\n")
    dic = {}
    i = 0
    while funcs[i] != '###' : 
            couple = fromSingle(funcs[i])
            dic[couple[0]] = couple[1]
            i += 1
    for j in range(i + 1, len(funcs)):
        filtList = fromGroup(dic, funcs[j])
        dic[filtList[0]] = filtList[1]
    return dic

#fourier transform filters functions

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
    arr1 = np.stack([np.arange(p//2) for i in range(q)], axis = 1)
    arr1 = np.concatenate((arr1, np.flipud(arr1)), axis = 0)
    arr2 = np.stack([np.arange(q//2) for i in range(p)], axis = 0)
    arr2 = np.concatenate((arr2, np.fliplr(arr2)), axis = 1)
    arr = (arr1 * q) ** 2 + (arr2 * p) ** 2
    for i in range(2, len(shape)):
        arr = np.stack([arr for i in shape[i]], axis = -1)
    return arr

def circleCut(shape, p0):
    return 20000 * ellipse(shape) <= (shape[0] * shape[1] * p0) ** 2

def BW(shape, n, p0):
    return (1 / (1 + (ellipse(shape) * 100 ** 2 / (shape[0] * shape[1] * p0) ** 2)))

def gaussCircle(shape, p0):
    return np.exp(-ellipse(shape) * 100 ** 2 / (shape[0] * shape[1] * p0) ** 2)

def gauss(shape, p0, angle):
    quad = angle // 90
    if quad == 0:
        t = (1,1)
    elif quad == 1:
        t = (-1, 1)
    elif quad == 2:
        t = (-1, -1)
    else:
        t = (1, -1)
    ang = angle % 180
    u = 

