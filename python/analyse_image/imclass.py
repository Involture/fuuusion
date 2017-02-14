import PIL.Image as img
import numpy as np
import matplotlib.pyplot as plt
emptya = np.array(None)
B = 32 #number of bits of the hardware
R = 2 ** B #range of the integers used
N = 256 #number of intensity levels in picture

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

#the actual class

class Image:
    """A custom Image class, mainly an RGB array."""

    def __init__(self, name, array, shape = None, _isp = emptya, _fsp = emptya):
        if shape == None:
            self.shape = np.shape(array)
        else:
            self.shape = shape
        self.name = name
        s = np.min(array)
        e = np.max(array)
        self.array = (array - s) * (N // e - s)
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
        self._isp = np.zeros((N,3))
        for pix in self.array:
            for colour in range(3):
                self._isp[, colour] += 1

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
        arr = np.zeros(np.shape(self.array))
        arr[:,:,colour] = np.copy(self.array[:,:,colour])
        isp = emptya
        fsp = emptya
        if self._ispExist():
            isp = np.zeros((N,3))
            isp[:,colour] = np.copy(self._isp[:,colour])
        if self._fspExist():
            fsp = np.zeros(self.shape)
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
    arr = np.array(im)
    return Image(fileName, arr)
