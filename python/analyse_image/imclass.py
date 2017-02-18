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

def _plotisp(subp, spectrum, c, n):
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
    """A custom Image class, mainly an RGB array."""

    def __init__(self, name, array, shape = None, _isp = emptya, _fsp = emptya):
        """Takes an unsigned little integer array as an argument."""
        isul(array)
        self.array = expend(array)
        if shape == None:
            self.shape = np.shape(array)
        else:
            self.shape = shape
        self.name = name
        if (_isp != emptya).any():
            self._isp = _isp
        if (_fsp != emptya).any():
            self._fsp = _fsp

    def __str__(self):
        return ("%s ; size %s"%(self.name, self.shape))

    def __repr__(self):
        return ("%s ; size %s"%(self.name,self.shape))

    def _initIsp(self):
        """Initialize the itensity spectrum."""
        self._isp = np.zeros((N,3), dtype = ubint)
#        for i in range(self.shape[0]):
#            for j in range(self.shape[1]):
#                for colour in range(3):
#                    self._isp[self.array[i,j,colour], colour] += 1
        for i in range(N):
            self._isp[i] = np.sum(self.array == i, axis = (0, 1))

    def isp(self):
        """Return the intensity spectrum of the image"""
        try:
            return self._isp
        except AttributeError:
            self._initIsp()
            return self._isp

    def _ispExist(self):
        """Return True only if the intensity spectrum is already initialized."""
        try:
            ignore = self._isp
            return True
        except AttributeError:
            return False

    def _fspExist(self):
        """Return True only if the frenquency spectrum is already initialized."""
        try:
            ignore = self._fsp
            return True
        except AttributeError:
            return False

    def extract(self, colour):
        """Return an image with only one chromatic componant."""
        arr = np.zeros(self.shape, dtype = ulint)
        arr[:,:,colour] = np.copy(self.array[:,:,colour])
        isp = emptya
        fsp = emptya
        if self._ispExist():
            isp = np.zeros((N,3), dtype = ubint)
            isp[:,colour] = np.copy(self._isp[:,colour])
        if self._fspExist():
            fsp = np.zeros(self.shape, dtype = ulint)
            fsp[:,:,colour] = np.copy(self._fsp[:,:,colour])
        return Image(self.name + "%d"%colour, arr, self.shape, isp, fsp)

    def show(self):
        """Display the image."""
        fig = plt.figure()
        subp = fig.add_subplot(111)
        subp.imshow(self.array)
        subp.axis('off')
        fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
        fig.show()

    def showc(self, colour):
        """Display only one chromatic color of the image."""
        self.extract(colour).show()

    def showisp(self):
        """Display the intensity spectrum of the image."""
        isp = self.isp()
        fig = plt.figure()
        sub1 = fig.add_subplot(131)
        sub2 = fig.add_subplot(132)
        sub3 = fig.add_subplot(133)
        _plotisp(sub1, isp[:,0], "red", self.name)
        _plotisp(sub2, isp[:,1], "green", self.name)
        _plotisp(sub3, isp[:,2], "blue", self.name)
        fig.show()

    def save(self, fileName):
        img.fromarray(self.array).convert('RGB').save(fileName)

    def copy(self, name):
        return Image(name, np.copy(self.array), self.shape, np.copy(self._isp), np.copy(self._fsp))

#opening fonction

def open(fileName):
    im = img.open(fileName)
    arr = np.array(im, dtype = ulint)
    return Image(fileName, arr)
