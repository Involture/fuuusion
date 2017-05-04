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

def toArg(l):
    split = l
    for i in range(len(split)):
        if ',' in split[i]:
            splitsplit = split[i].split(',')
            split[i] = tuple([int(x) for x in splitsplit])
        else:
            split[i] = int(split[i])
    return tuple(split)

def fromFunc(func):
    split = func[1:].split(" ")
    return funcNames[split[0]](*toArg(split[1:]))
    
def fromSingle(string):
    split = string.split("\n")
    if split[1][0] == '>':
        return (split[0], fromFunc(split[1]))
    else:
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

def vdist(shape):
    p, q = shape[:2]
    arr1 = np.stack([np.arange(p // 2) for i in range(q)], axis = 1)
    arr2 = np.stack([np.arange(p // 2 + (p % 2)) for i in range(q)], axis = 1)
    return np.concatenate((arr1, np.flipud(arr2)), axis = 0)

def permute(arr):
    row = arr[0]
    temp = row[-1]
    arr[0][1:] = row[:-1]
    arr[0][0] = temp

def ddist(size):
    row = np.array(np.concatenate((np.arange(size), np.arange(size - 2, 0, -1)), axis = 0), ndmin = 2)
    res = np.array(row)
    for j in range(size - 1):
        permute(row)
        res = np.concatenate((res, row), axis = 0)
    return res[:,:size]

def ellipse(shape):
    p, q = shape[:2]
    arr = (vdist((p,q)) * q) ** 2 + (vdist((q,p)).T * p) ** 2
    for i in range(2, len(shape)):
        arr = np.stack([arr for i in shape[i]], axis = -1)
    return arr

def circleCut(shape, p0):
    return 20000 * ellipse(shape) <= (shape[0] * shape[1] * p0) ** 2

def BW(shape, n, p0):
    return (1 / (1 + (ellipse(shape) * 160000 / (shape[0] * shape[1] * p0) ** 2)))

def cgauss(shape, p0):
    return np.exp(-ellipse(shape) * 160000 / (shape[0] * shape[1] * p0) ** 2)

def vgauss(shape, p0):
    return np.exp(-vdist(shape) ** 2 * 40000 / (shape[0] * p0)  ** 2)

def hgauss(shape, p0):
    p, q = shape[:2]
    return hgauss((q, p), p0).T

def d1gauss(shape, p0):
    assert(shape[0] == shape[1])
    return np.exp(-ddist(shape[0]) ** 2 * 40000 / (shape[0] * p0) ** 2)

def d2gauss(shape, p0):
    return np.flipud(d1gauss(shape, p0))

funcNames = locals()
