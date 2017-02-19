from filterBase import *

def equalize(im):
    arr = np.copy(im.array)
    ispRef = [[],[],[]] #list containing indexes of pixels sorted by intensity value
    if im._ispExist() :
        isp = np.copy(im.isp())
    else:
        isp = np.zeros((N,3), dtype = ubint)
    #initializing ispRef lists
    for i in range(N):
        for colour in range(3):
            ispRef[colour].append([])
    #initializing ispRef
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            for colour in range(3):
                intensity = arr[i,j,colour]
                ispRef[colour][intensity].append((i,j))
                if not im._ispExist():
                    isp[intensity, colour] += 1

    nbPix = im.shape[0] * im.shape[1]
    q = nbPix // N
    r = nbPix % N
    equalRep = np.full((N), q, dtype = ubint) #desired intensity spectrum
    #intializing equalRep
    for i in range(r):
        equalRep += 1

    newIsp = np.zeros((N,3), dtype = ubint) #the isp of the futur equalized image
    for colour in range(3):
        remaining = equalRep[0] 
        newI = 0
        for intensity in range(1):
            currentIPixNb = isp[intensity, colour]
            #loop until desired newI is found
            while  currentIPixNb > remaining:
                currentIPixNb -= remaining
                newI += 1
                remaining = equalRep[newI]
            remaining -= currentIPixNb
            #updating newIsp and changing pixels value to match the desired repartition
            newIsp[newI, colour] += isp[intensity, colour]
            for pix in ispRef[colour][intensity]:
                arr[pix[0], pix[1], colour] = newI

    #setting the argument image isp in case it wasn't initialized
    if not im._ispExist():
        im._isp = isp
        
    return Image(im.name + '_equalized', arr, im.shape, newIsp)

#base operations on images

def operate(im1, im2, f, opName):
    if im1.shape != im2.shape:
        raise ValueError("can't add images with different sizes")
    shape = im1.shape
    arr1 = bint(np.copy(im1.array))
    arr2 = bint(np.copy(im2.array))
    return Image(opName + "_" + im1.name + "_" + im2.name, f(arr1, arr2), shape)

def add(im1, im2):
    return operate(im1, im2, (lambda x,y : x + y), "+")

def sub(im1, im2):
    return operate(im1, im2, (lambda x,y : x - y), "-")

def max(im1, im2):
    return operate(im1, im2, np.maximum, "max")

def min(im1, im2):
    return operate(im1, im2, np.minimum, "min")

def diff(im1, im2):
    return operate(im1, im2, (lambda x,y : np.maximum(im1, im2) - np.minimum(im1, im2)), "diff")

def neg(im):
    shape = im.shape
    arr = np.copy(im.array)
    isul(arr)
    return Image("neg_" + im.name, -(arr + 1), shape)

#algebra operations

def makeFilter(set):
    if set.dtype != np.bool:
        raise AttributeError("shape must be a boolean array")
    shape = np.shape(set)
    dx, dy = shape[0], shape[1]
    bigShape =(dx * 2 - 1, dy * 2 - 1)
    result = np.full(bigShape, False,  dtype = np.bool)
    for i in range(dx):
        for j in range(dy):
            if set[dx - i - 1, dy - j - 1]:
                result[i : dx + i, j : dy + j] += set 
    return result

def algOp(im, func, baseValue, set):
    setFilter  = lint(makeFilter(set))
    arr = ubint(im.array)
    arr = filterArr(arr, setFilter, func, baseValue)
    return Image("closed_" + im.name, reduce(arr), im.shape)

def algClose(im, set):
    return algOp(im, np.maximum, 0, set)

def algOpen(im, set):
    return algOp(im, np.minimum, 255, set)

def isLink(boolArr):
    

def squeletize(binIm):
    arr = np.bool(binIm.array)


#filtering functions

def isFilter(filterArr):
    issl(filterArr)
    filterSize = np.shape(filterArr)
    if (filterSize[0] % 2 != 1) or (filterSize[1] % 2 != 1):
        raise ValueError("can't filter with an even sized filter array")

def addFilter(firstFilter, secondFilter):
    """Return a new filter equivalent to applying successively the first and the second filter WITHOUT TAKING THE ABSOLUTE VALUE"""
    isFilter(firstFilter)
    isFilter(secondFilter)
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


def filterArr(sourceArr, filterArr, func, baseValue):
    """Take a unsigned big integer array as argument. Return the absolute value of the filtered array."""
    isFilter(filterArr)
    isub(sourceArr)
    arr = np.copy(sourceArr)
    shape = np.shape(sourceArr)
    filterSize = np.shape(filterArr)
    dx = (filterSize[0] - 1) // 2
    dy = (filterSize[1] - 1) // 2
    bigShape =(shape[0] + 2 * dx, shape[1] + 2 * dy, 3)
    result = np.full(bigShape, baseValue, dtype = bint)
    for i in range(-dx, dx + 1):
        for j in range(-dy, dy + 1):
            if filterArr[dx - i, dy - j] != 0:
                shiftedArr = np.zeros(bigShape, dtype = bint)
                shiftedArr[dx + i : dx + i + shape[0], dy + j : dy + j + shape[1], :] = filterArr[dx - i, dy - j] * arr
                result = func(result, shiftedArr)
    if np.min(filterArr) < 0:
        result = np.abs(result)
    return ubint(result[dx : shape[0] + dx, dy : shape[1] + dy, :])

def filter(im, filterList):
    """Filter the array of the image successively with all the filter contained in filterList"""
    arr = ubint(im.array)
    for filter in filterList:
        arr = filterArr(arr, filter, np.add, 0)
    return Image("filtered_" + im.name, reduce(arr), im.shape)

def filterSum(im, filterList):
    arr = ubint(im.array)
    result = filterArr(arr, filterList[0], np.add, 0)
    for filter in filterList[1:]:
        result += filterArr(arr, filter, np.add, 0)
    return Image("filtered_" + im.name, reduce(result), im.shape)

#transforming an image to binary image

def greyAv(im):
    array = im.array
    arrayAv = np.sum(ubint(array), dtype = ubint, axis = 2) // 3
    greyArr = np.zeros(im.shape, dtype = ubint)
    for colour in range(3):
        greyArr[:,:,colour] = np.copy(arrayAv)
    return Image("greyAv_" + im.name, reduce(greyArr), im.shape)

def greyMax(im):
    array = im.array
    arrayMax = np.max(array, axis = 2)
    greyArr = np.zeros(im.shape, dtype = ulint)
    for colour in range(3):
        greyArr[:,:,colour] = np.copy(arrayMax)
    return Image("greyMax_" + im.name, greyArr, im.shape)

def binary(im, cutIntensity):
    arr = np.zeros(im.shape, dtype = ulint)
    for colour in range(3):
        arr[:,:,colour] = ulint(im.array[:,:,0] > cutIntensity) * 255
    return Image("bin_" + im.name, arr, np.shape)

#maximize intensity while keeping color

def colorize(im):
    array = np.copy(im.array)
    maxCompArray = np.max(array, axis = 2)
    divideArray = np.zeros(im.shape, dtype = ulint)
    for colour in range(3):
        divideArray[:,:,colour] = maxCompArray
    divideArray = np.maximum(divideArray, np.ones(im.shape, dtype = ubint))
    array = (ubint(array) * 255) // divideArray
    return Image("colorized_" + im.name, reduce(array), im.shape)













