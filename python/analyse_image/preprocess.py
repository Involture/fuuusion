import png
import numpy as np
import matplotlib.pyplot as plt
Bcalc = 32 #number of bits of the hardware
Bim = 8 #number of bits of the colors
R = 2 ** Bcalc #range of the integers used for calculus
N = 2 ** Bim #number of intensity levels in picture
ulint = np.uint8 #numpy unsigned integer type corresponding
ubint = np.uint32 #numpy unsigned integer type corresponding
lint = np.int8 #numpy signed integer type corresponding
bint = np.int32 #numpy signed integer type corresponding

#opening fonction

def open(fileName):
    tuple = png.Reader(fileName).asRGB8()
    flatSize = (tuple[1], tuple[0] * 3)
    iter = tuple[2]
    flatArray = np.zeros(flatSize, dtype = ulint)
    i = 0
    for arr in iter:
        flatArray[i] = np.fromiter(arr, dtype = ulint)
        i += 1
    size = (tuple[1], tuple[0], 3)
    array = np.zeros(size, dtype = ulint)
    for colour in range(3):
        array[:,:,colour] = flatArray[:,colour::3]
    return array

#array conversion and type checking functions

def isul(array):
    return array.dtype == ulint
def isub(array):
    return array.dtype == ubint
def issl(array):
    return array.dtype == lint
def issb(array):
    return array.dtype == bint
def isbin(array):
    return array.dtype == np.bool

def cisul(array):
    if not isul(array):
        raise TypeError("not an unsigned little integer array")
def cisub(array):
    if not isub(array):
        raise TypeError("not an unsigned big integer array")
def cissl(array):
    if not issl(array):
        raise TypeError("not a signed little integer array")
def cissb(array):
    if not issb(array):
        raise TypeError("not a signed big integer array")
def cisbin(array):
    if not isbin(array):
        raise TypeError("not a boolean array")

def reduce(array):
    """Convert an unsigned big integer array into an unsigned little integer using the maximum intensity range"""
    cisub(array)
    s = np.min(array)
    e = np.max(array)
    r = e - s
    newArr = np.copy(array)
    if r != N:
        if np.float64(e) * N < R: 
            newArr = ulint((newArr * (N - 1)) // r)
        else:
            newArr = ulint(newArr * ((N - 1) // r))
    return newArr

def expend(array):
    cisul(array)
    s = np.min(array)
    e = np.max(array)
    r = e - s
    newArr = np.copy(array)
    if r < N - 1:
        newArr = ulint((ubint(newArr) * (N - 1)) // r)
    return newArr

#ploting functions

def plotIsp(subp, spectrum, c, n):
    """Plot the intensity spectrum to the specified subplot
    c is the color of the plot
    n is its name
    """
    ind = np.arange(N)
    width = 1

    subp.bar(ind, spectrum, width, color = c, edgecolor = c)

    subp.set_title(n)
    subp.set_xlabel("pixel value")
    subp.set_ylabel("nb of pixels")

def show(arr):
    """Display the image."""
    size = np.shape(arr)
    if len(size) == 2:
        array = np.zeros((size[0], size[1], 3), dtype = arr.dtype)
        for colour in range(3):
            array[:,:,colour] = arr
    else:
        array = np.copy(arr)
    if isul(array):
        array = array / 255
    elif not isbin(array):
        raise TypeError("must be an unsigned little integer or a boolean array to be plot")
    fig = plt.figure()
    subp = fig.add_subplot(111)
    subp.imshow(array)
    subp.axis('off')
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
    fig.show()

#type converting functions

def greyAv(arr):
    return ulint(np.sum(ubint(arr), dtype = ubint, axis = 2) // 3)

def greyMax(arr):
    return  np.max(arr, axis = 2)

def binary(arr, cutIntensity):
    return arr > cutIntensity

#base operations on images

def operate(arr1, arr2, func):
    if not ((isul(arr1) and isul(arr2)) or (isub(arr1) and isub(arr2))) :
        raise TypeError("not same type or not integer array")
    if np.shape(arr1) != np.shape(arr2) :
        raise ValueError("can't operate on images with different sizes")
    return func(arr1, arr2)

def add(arr1, arr2):
    return operate(ubint(arr1), ubint(arr2), (lambda x,y : x + y))

def sub(arr1, arr2):
    return operate(ubint(arr1), ubint(arr2), (lambda x,y : x - y))

def max(arr1, arr2):
    return operate(ubint(arr1), ubint(arr2), np.maximum)

def min(arr1, arr2):
    return operate(ubint(arr1), ubint(arr2), np.minimum)

def diff(arr1, arr2):
    return operate(ubint(arr1), ubint(arr2), (lambda x,y : np.maximum(arr1, arr2) - np.minimum(arr1, arr2)))

def neg(arr):
    cisul(arr)
    array = np.copy(arr)
    return -(arr + 1)

def operateBin(bin1, bin2, func):
    cisbin(bin1)
    cisbin(bin2)
    if np.shape(bin1) != np.shape(bin2):
        raise ValueError("can't operate on images with different size")
    return func(bin1, bin2)

def algOr(bin1, bin2):
    return operateBin(bin1, bin2, (lambda x,y : x + y))

def algAnd(bin1, bin2):
    return operateBin(bin1, bin2, (lambda x,y : np.uint8(x) * np.uint8(y) == 1))

def algXor(bin1, bin2):
    return operateBin(bin1, bin2, (lambda x,y : np.uint8(x) + np.uint8(y) == 1))

def algNot(bin):
    cisbin(bin)
    return bin == False

#algebra operations

def makeFilter(set):
    cisbin(set)
    shape = np.shape(set)
    dx, dy = shape[0], shape[1]
    bigShape =(dx * 2 - 1, dy * 2 - 1)
    result = np.full(bigShape, False,  dtype = np.bool)
    for i in range(dx):
        for j in range(dy):
            if set[dx - i - 1, dy - j - 1]:
                result[i : dx + i, j : dy + j] += set 
    return result

def algClose(bin, set):
    filter = makeFilter(set)
    return filterArr(bin, filter, (lambda x,y : y), algOr, False, np.bool)

def algOpen(bin, set):
    filter = (makeFilter(set))
    return filterArr(bin, filter, (lambda x,y : y), algAnd, True, np.bool)

def path(ind):
    if ind <= 2:
        return (0,ind)
    elif ind <= 4:
        return (ind - 2, 2)
    elif ind <= 6:
        return (2, 6 - ind)
    elif ind <= 8:
        return (8 - ind, 0)
    else:
        raise ValueError("path index out of range (0-8)")

def isFrontier(bin):
    frontierFilter = np.array([[False, True, False],[True, True, True], [False, True, False]])
    return algXor(bin, filterArr(bin, frontierFilter, (lambda x,y : y), algAnd, True, np.bool))

def isLink(bin):
    shape = np.shape(bin)
    Dx, Dy = shape
    bigShape = (Dx + 2, Dy + 2)
    changeCount = np.zeros(shape, dtype = np.uint8)
    currentVal = np.full(bigShape, True, dtype = np.bool)
    currentVal[2:, 2:] = np.copy(bin)
    for i in range(1, 9):
        ind = path(i)
        changeCount += np.uint8(algXor(currentVal[2 - ind[0] : 2 - ind[0] + Dx, 2 - ind[1] : 2 - ind[1] + Dy], currentVal[1 : 1 + Dx, 1 : 1 + Dy]))
        currentVal[2 - ind[0] : 2 - ind[0] + Dx, 2 - ind[1] : 2 - ind[1] + Dy] = np.copy(bin)
    return (currentVal >= 3)[1 : 1 + Dx, 1 : 1 + Dy]

def diminish(bin):
    return algAnd(algOr(algNot(Frontier), isLink), bin)
    

#filtering functions

def isFilter(filterArr):
    filterSize = np.shape(filterArr)
    if (filterSize[0] % 2 != 1) or (filterSize[1] % 2 != 1):
        raise ValueError("can't filter with an even sized filter array")

def filterArr(arr, filterArr, filterFunc, addFunc, baseValue, arrType, finalFunc = None) :
    isFilter(filterArr)
    shape = np.shape(arr)
    filterSize = np.shape(filterArr)

    dx = (filterSize[0] - 1) // 2
    dy = (filterSize[1] - 1) // 2
    Dx = shape[0] + 2 * dx
    Dy = shape[1] + 2 * dy
    if len(shape) == 3:
        bigShape = (Dx, Dy, 3)
    else:
        bigShape = (Dx, Dy)

    result = np.full(bigShape, baseValue, dtype = arrType)
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            if ulint(filterArr[dx - i, dy - j]) != 0:

                shiftedArr = np.zeros(bigShape, dtype = arrType)
                if len(shape) == 3:
                    shiftedArr[dx + i : dx + i + shape[0], dy + j : dy + j + shape[1], :] = filterFunc(filterArr[dx - i, dy - j], arr)
                else:
                    shiftedArr[dx + i : dx + i + shape[0], dy + j : dy + j + shape[1]] = filterFunc(filterArr[dx - i, dy - j], arr)
                result = addFunc(result, shiftedArr)

    if finalFunc != None:
        result = finalFunc(result)
    if len(shape) == 3:
        return result[dx : shape[0] + dx, dy : shape[1] + dy, :]
    else:
        return result[dx : shape[0] + dx, dy : shape[1] + dy]

def filterSum(arr, filterList):
    """filter the array wit heach filter and add the result.
    The array must be a signed big integer array."""
    newFilterList = reduceFilterList(filterList)
    result = filterArr(arr, filterList[0], (lambda x,y : x * y), np.add, 0, np.abs)
    for filter in filterList[1:]:
        cissl(filter)
        result += filterArr(arr, filter, (lambda x,y : x * y), np.add, 0, np.abs)
    return Image("filtered_" + im.name, reduce(result), im.shape)

def addFilter(firstFilter, secondFilter):
    """Return a new filter equivalent to applying successively the first and the second filter.
    Filters must be signed little integers arrays of odd size"""
    isFilter(firstFilter)
    isFilter(secondFilter)
    cissl(firstFilter)
    cissl(secondFilter)
    shape1 = np.shape(firstFilter)
    shape2 = np.shape(secondFilter)
    dx = (shape1[0] - 1) // 2
    dy = (shape1[1] - 1) // 2
    bigShape =(shape2[0] + 2 * dx, shape2[1] + 2 * dy)
    result = np.zeros(bigShape, dtype = lint)
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            shiftedArr = np.zeros(bigShape, dtype = lint)
            shiftedArr[dx + i : dx + i + shape2[0], dy + j : dy + j + shape2[1]] = firstFilter[dx - i, dy - j] * secondFilter
            result += shiftedArr
    return result

def reduceFilterList(filterList):
    """add consecutive positive filters in order to optimize calculus time.
    operate on a signed little integer filterList"""
    newFilterList = []
    i = 1
    filter = filterList[0]
    while i < len(filterList):
        nextFilter = filterList[i]
        if type(filter) != list:
            cissl(filter)
            if type(nextFilter) != list and np.min(nextFilter) >= 0 and np.min(filter) >=0:
                nextFilter = addFilter(filter, nextFilter)
            else:
                newFilterList.append(filter)
        else:
            newFilterList.append(filter)
        i += 1
        filter = nextFilter
    newFilterList.append(filter)
    return newFilterList
                

def filter(arr, filterList):
    """Filter the array with all the filters contained in filterList.
    Array must be an unsigned little integer array."""
    cisul(arr)
    array = bint(np.copy(arr))
    for filter in reduceFilterList(filterList):
        if type(filter) == list:
            filterSum(array, filter)
        else:
            array = filterArr(array, filter, (lambda x,y : x * y), np.add, 0, bint, np.abs)
    return reduce(ubint(np.abs(array)))

#maximize intensity while keeping color

def colorize(arr):
    cisul(arr)
    array = np.copy(arr)
    shape = np.shape(array)
    maxCompArray = np.max(array, axis = 2)
    divideArray = np.zeros(shape, dtype = ulint)
    for colour in range(3):
        divideArray[:,:,colour] = maxCompArray
    divideArray = np.maximum(divideArray, np.ones(shape, dtype = ubint))
    array = (ubint(array) * 255) // divideArray
    return array

#im1 = open("image1.png")
#im2 = open("image2.png")
#blur = np.ones((5,5), dtype = lint)
#corner = np.array([[1,0,1],[0,0,0],[1,0,1]], dtype = lint)
#grad = [np.array([[1,1,1], [0,0,0], [-1, -1, -1]], dtype = lint), np.array([[1,0,-1], [1,0,-1], [1,0,-1]], dtype = lint)]
#im3 = filter(im1, grad)
#im4 = greyAv(im3)
#bin1 = binary(im4, 20)
testbin = np.full((1000,1000), False, dtype = np.bool)
testbin[500] = np.full((1000), True, dtype = np.bool)
#set1 = np.array([[True,True],[True,True]])

