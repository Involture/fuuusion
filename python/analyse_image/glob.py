import numpy as np
import time as t
from sys import getsizeof
from subprocess import check_output

#custom types

cint = np.uint8
cflt = np.float16
ccpl = np.complex64

#Prints and verbose

start = t.time()

V = True
VV = True

def p(s):
    """Print only if V (standing for a verbose parametre) is True."""
    if V:
        print(s)

def pp(s):
    """Print only if VV (standing for a very verbose parameter) is True."""
    if VV:
        print(s)

def ptime():
    """Print the total time since the beginning of the program if V is True."""
    p("time " + str(t.time() - start))

def pptime():
    """Print the total time since the beginning of the program if VV is True."""
    pp("time" + str(t.time() - start))

def pmem():
    """Print the total amount of memory used by python processes in MB \
    if V is True."""
    p("mem " + str(check_output("ps aux|grep python|awk '{sum=sum+$6}; END {print sum/1024}'", shell = True))[2:-3])

def ppmem():
    """Print the total amount of memory used by python processes in MB \
    if VV is True."""
    pp("mem " + str(check_output("ps aux|grep python|awk '{sum=sum+$6}; END {print sum/1024}'", shell = True))[2:-3])


def psize(a):
    """Print the size of and obect in MB if V is True."""
    p("size " + str(getsizeof(a)/1024**2))

def ppsize(a):
    """Print the size of and obect in MB if VV is True"""
    pp("size " + str(getsizeof(a)/1024**2))

#type checkers

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
            powCheck(n // 2)

def cispow(arr):
    """raise an error if one of the first dimension of arr is not a power of two"""
    p, q = np.shape(arr)[:2]
    powCheck(p)
    powCheck(q)

#???

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
