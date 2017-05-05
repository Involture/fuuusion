import numpy as np
import time as t
from sys import getsizeof
from subprocess import check_output

Bcalc = 32 #number of bits of the hardware
Bim = 8 #number of bits of the colors
R = 2 ** Bcalc #range of the integers used for calculus
N = 2 ** Bim #number of intensity levels in picture

ulint = np.uint8 #numpy unsigned integer type corresponding
ubint = np.uint32 #numpy unsigned integer type corresponding
lint = np.int8 #numpy signed integer type corresponding
bint = np.int32 #numpy signed integer type corresponding

#measuring time and memory functions

start = t.time()

def p(s):
    if VERBOSE:
        print(s)

def ptime():
    p("time " + str(t.time() - start))
def pmem():
    p("mem " + str(check_output("ps aux|grep python|awk '{sum=sum+$6}; END {print sum/1024}'", shell = True))[2:-3])
def psize(a):
    p("size " + str(getsizeof(a)/1024**2))

#expand to match

def expand(arr, n):
    res = arr
    for i in range(n):
        res = np.expand_dims(res, axis = -1)
    return res

#type checking functions

def isul(array):
    return array.dtype == ulint
def isub(array):
    return array.dtype == ubint
def issl(array):
    return array.dtype == lint
def issb(array):
    return array.dtype == bint

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

def cisfilt(filterArr):
    filterSize = np.shape(filterArr)
    if (filterSize[0] % 2 != 1) or (filterSize[1] % 2 != 1):
        raise ValueError("can't filter with an even sized filter array")

def powCheck(n):
    if n != 1:
        if n % 2 == 1:
            raise ValueError("array dimensions are not powers of two")
        else:
            powCheck(n // 2)
def cispow(arr):
    p, q = np.shape(arr)[:2]
    powCheck(p)
    powCheck(q)

