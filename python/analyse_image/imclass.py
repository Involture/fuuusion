import PIL.Image as img
import numpy as np
import matplotlib.pyplot as plt
emptya = np.array(None)
Bcalc = 32 #number of bits of the hardware
Bim = 8 #number of bits of the colors
R = 2 ** Bcalc #range of the integers used for calculus
N = 2 ** Bim #number of intensity levels in picture
ulint = np.uint8 #numpy unsigned integer type corresponding
ubint = np.uint32 #numpy unsigned integer type corresponding
lint = np.int8 #numpy signed integer type corresponding
bint = np.int32 #numpy signed integer type corresponding

#ploting functions

def _plotIsp(subp, spectrum, c, n):
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

def _plotArr(arr):
    """Display the image."""
    fig = plt.figure()
    subp = fig.add_subplot(111)
    subp.imshow(np.float64(arr) / 255)
    subp.axis('off')
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
    fig.show()

#array of integers converting functions

def isul(array):
    if not array.dtype == ulint:
        raise AttributeError("not an unsigned little integer array")
def isub(array):
    if not array.dtype == ubint:
        raise AttributeError("not an unsigned big integer array")
def issl(array):
    if not array.dtype == lint:
        raise AttributeError("not a signed little integer array")
def issb(array):
    if not array.dtype == bint:
        raise AttributeError("not a signed big integer array")
    """Take an unsigned little array as argument.
    Return an unsigned little array which intensity spectrum as been expended at the maximum range."""
    isul(array)
    s = np.min(array)
    e = np.max(array)
    r = e - s
    newArr = np.copy(array)
    if r < N:
        newArr = ulint((ubint(newArr) * N - 1) // r)
    return newArr

def reduce(array):
    """Convert an unsigned big integer array into an unsigned little integer using the maximum intensity range"""
    isub(array)
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
    isul(array)
    s = np.min(array)
    e = np.max(array)
    r = e - s
    newArr = np.copy(array)
    if r < N - 1:
        newArr = ulint((ubint(newArr) * (N - 1)) // r)
    return newArr

#the actual class

class Image:
    """A custom Image class. Either a RGB or grey levels image"""

    def __init__(self, name, array, isColour):
        """Takes an unsigned little integer array as an argument."""
        isul(array)

        self.shape = np.shape(array)
        dim = len(self.shape)
        self.isColour = False
        if dim == 3:
            self.isColour = True
        elif dim != 2:
            raise AttributeError("dimension of the array isn't 2 or 3")

        self.array = expend(array)
        self.name = name

    def __str__(self):
        return ("%s ; size %s"%(self.name, self.shape))

    def __repr__(self):
        return ("%s ; size %s"%(self.name,self.shape))

    def isp(self):
        """Return the intensity spectrum"""
        if self.isColour:
            isp = np.zeros((N,3), dtype = ubint)
        else:
            isp = np.zeros((N), dtype = ubint)
        for i in range(N):
            isp[i] = np.sum(self.array == i, axis = (0, 1))
        return isp

    def show(self):
        _plotArr(self.array)

    def showc(self, colour):
        """Display only one chromatic color of the image."""
        if not self.isColour:
            raise TypeError("not a coloured image")
        _plotArr(self.array[:,:,colour])

    def showisp(self):
        """Display the intensity spectrum of the image."""
        isp = self.isp()
        fig = plt.figure()
        if self.isColour:
            sub1 = fig.add_subplot(131)
            sub2 = fig.add_subplot(132)
            sub3 = fig.add_subplot(133)
            _plotIsp(sub1, isp[:,0], "red", self.name)
            _plotIsp(sub2, isp[:,1], "green", self.name)
            _plotIsp(sub3, isp[:,2], "blue", self.name)
        else:
            sub = fig.add_subplot(111)
            _plotIsp(sub, isp, "black", self.name)
        fig.show()

    #def save(self, fileName = ""):
        #if fileName == "":
            #fileName = self.name
        #img.fromarray(self.array).convert('RGB').save(fileName)

class BinImage:

    def __init__(self, name, array):
        self.name = name
        
        shape = np.shape(array)
        if shape != 2:
            raise AttributeError("array dimension must be two")
        self.shape = shape

        if array.dtype != np.bool:
            raise AttributeError("array must be boolean")
        self.array = array

    def __str__(self):
        return ("%s ; size %s"%(self.name, self.shape))

    def __repr__(self):
        return ("%s ; size %s"%(self.name,self.shape))

    def show(self):
        _plotArr(np.float(self.array))

#opening fonction

def open(fileName):
    im = img.open(fileName)
    arr = np.array(im, dtype = ulint)
    return Image(fileName, arr, True)
