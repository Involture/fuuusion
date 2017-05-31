import numpy as np
import time as t
import itertools
from sys import getsizeof
from subprocess import check_output

#custom types

cint = np.uint8
cflt = np.float32
ccpl = np.complex64

#Prints and verbose

start = t.time()

V = 3

def pr(s, depth):
    if V >= depth:
        print(s)

def ptime(depth):
    pr("time " + str(t.time() - start), depth)

def pmem(depth):
    pr("mem " + str(
        check_output(
            "ps aux|grep python|awk '{sum=sum+$6}; \
            END {print sum/1024}'",
            shell = True
            )
        )[2:-3],
        depth)

def psize(arr, depth):
    pr(arr.nbytes / 1024, depth)

#ype checkers

def cisfilt(filterArr):
    """raise an error if one of the first two dimension of filterArr is even"""
    filterSize = np.shape(filterArr)
    if (filterSize[0] % 2 != 1) or (filterSize[1] % 2 != 1):
        raise ValueError("can't filter with an even sized filter array")

def powCheck2(n):
    """raise an error if n is not a power of two"""
    if n != 1:
        if n % 2 == 1:
            raise ValueError("array dimensions are not powers of two")
        else:
            powCheck2(n // 2)

def cispow(arr):
    """raise an error if one of the first dimension of arr is not a power of two"""
    p, q = np.shape(arr)[:2]
    powCheck2(p)
    powCheck2(q)

#basic geometric transformation on arrays
def permute(l, i, j):
    if i != j:
        l[i], l[j] = l[j], l[i]

def transposeEnd(ndim, dims):
    dimList = list(range(ndim))
    dims.sort(reverse = True)
    l = len(dims)
    for i in range(l):
        ind = dims[i]
        lastInd = ndim - i - 1
        permute(dimList, ind, lastInd)
    return tuple(dimList)

def multDim(arr, multArr, dims):
    """Multiply arr by multArr on its dimension dims, element-wise,\
    IN PLACE."""
    transTuple = transposeEnd(arr.ndim, dims)
    transArr = np.transpose(arr, transTuple)
    multArr = transArr * multArr
    return np.transpose(multArr, transTuple)

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

#basic data transformation on array

def red(arr, amp = 1.):
    amplified = arr * 255. * amp / arr.max()
    cut = np.minimum(arr, 255.)
    return cut

#concate or stack on generator

def stackInd(axis, naxis, i):
    if axis >= 0:
        if axis < naxis - 1:
            return (slice(None),) * axis + (i,) + (Ellipsis,)
        else:
            return (slice(None),) * axis + (i,)
    else:
        if -axis < naxis:
            return (Ellipsis,) + (i,) + (slice(None),) * (- axis - 1)
        else:
            return (i,) + (slice(None),) * (- axis - 1)

def stackBigShape(shape, geneLen, axis):
    if axis >= 0:
        return shape[:axis] + (geneLen,) + shape[axis:]
    else:
        if axis < -1:
            return shape[:axis + 1] + (geneLen,) + shape[axis + 1:]
        else:
            return shape + (geneLen,)

def concatInd(concatDimLen, axis, naxis,  i):
    if axis >= 0:
        if axis < naxis - 1:
            return (
                    (slice(None),) * axis + 
                    (slice(i * concatDimLen, (i + 1) * concatDimLen),) + 
                    (Ellipsis,)
                    )
        else:
            return (
                    (slice(None),) * axis + 
                    (slice(i * concatDimLen, (i + 1) * concatDimLen),) 
                    )
    else:
        if -axis < naxis:
            return (
                    (Ellipsis,) + 
                    (slice(i * concatDimLen, (i + 1) * concatDimLen),) + 
                    (slice(None),) * (- axis - 1)
                    )
        else:
            return (
                    (slice(i * concatDimLen, (i + 1) * concatDimLen),) + 
                    (slice(None),) * (- axis - 1)
                    )

def concatBigShape(shape, geneLen, axis):
    concatDimLen = shape[axis]
    if axis >= 0:
        if axis < len(shape) - 1:
            return (
                    shape[:axis] + 
                    (concatDimLen * geneLen,) + 
                    shape[axis + 1:]
                    )
        else:
            return (
                    shape[:axis] + 
                    (concatDimLen * geneLen,) 
                    )
    else:
        if axis < -1:
            return (
                    shape[:axis] + 
                    (concatDimLen * geneLen,) + 
                    shape[axis + 1:]
                    )
        else:
            return (
                    shape[:axis] + 
                    (concatDimLen * geneLen,)
                    )

def stackGene(gene0, geneLen, axis):

    first = next(gene0)
    shape = first.shape
    gene = itertools.chain((first,), gene0)

    naxis = len(shape) + 1
    assert((axis < naxis and axis >= 0) 
            or 
            (-axis <= naxis)
            )

    bigShape = stackBigShape(shape, geneLen, axis)
    bigArr = np.empty(bigShape)

    for i, arr in enumerate(gene):
        bigArr[stackInd(axis, naxis, i)] = arr
    return bigArr

def concatGene(gene0, geneLen, axis):

    first = next(gene0)
    shape = first.shape
    gene = itertools.chain((first,), gene0)

    naxis = len(shape)
    assert((axis < naxis and axis >= 0) or (-axis <= naxis))

    concatDimLen = shape[axis]
    bigShape = concatBigShape(shape, geneLen, axis)
    bigArr = np.empty(bigShape)

    for i, arr in enumerate(gene):
        bigArr[concatInd(concatDimLen, axis, naxis, i)] = arr
    return bigArr
