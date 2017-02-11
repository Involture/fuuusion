import PIL.Image as img
import numpy as np
import matplotlib.pyplot as plt

"ploting functions"

def _plotisp(subp, spectrum, c, n):
    """Plot the intensity spectrum to the specified subplot
    c is the color of the plot
    n is its name
    """
    N = 256
    ind = np.arange(N)
    width = 1

    subp.bar(ind, spectrum, width, color = c, edgecolor = c)

    subp.set_title(n)
    subp.set_xlabel("pixel value")
    subp.set_ylabel("nb of pixels")

"the actual class"

class Image:
    """A custom Image class, mainly an RGB array"""
    def __init__(self, name, array):
        s = np.shape(array)
        sx = s[1]
        sy = s[0]
        self.name = name
        self.array = array
        self.size = (sx, sy)
        self._initsp = np.zeros((2,1))
        self._isp = np.zeros((256, 3))
        self._fsp = np.zeros((sy, sx, 3))
    def __str__(self):
        return ("%s ; size %s")
    def _initIsp(self):
        if self._initsp[0] == 0:
            for colour in range(3):
                for pix in self.array[:,:,colour]:
                    self.isp[pix, colour] += 1
            self._initsp[0] = 1
    def show(self):
        fig = plt.figure()
        subp = fig.add_subplot(111)
        subp.imshow(self.array)
        subp.axis('off')
        fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
        fig.show()
    def __repr__(self):
        return ("%s ; size %s")
    def extract(self, colour):
        arr = np.zeros(np.shape(self.array))
        arr[:,:,colour] = np.copy(self.array[:,:,colour])
        colorIm = Image(self.name + "%d"%colour, arr)
        if self._initsp[0] == 1:
            colourIm._initsp[0] = 1
            colourIm._isp[:,:,colour] = np.copy(self._ips[:,:,colour])
        if self._initsp[1] == 1:
            colourIm._initsp[1] = 1
            colourIm._fsp[:,:,colour] = np.copy(self._ifs[:,:,colour])
        for i in range(2):
            for j in range(3):
                if j != colour:
                    colorIm._initsp = 1
        return colorIm
    def showc(self, colour):
        self.extract(colour).show()
    def isp(self):
        self._initIsp()
        return self._isp
    def showisp(self):
        self._initIsp()
        fig = plt.figure()
        sub1 = fig.add_subplot(131)
        sub2 = fig.add_subplot(132)
        sub3 = fig.add_subplot(133)
        _plotisp(sub1, self._isp[:,0], "red", self.name)
        _plotisp(sub2, self._isp[:,1], "green", self.name)
        _plotisp(sub3, self._isp[:,2], "blue", self.name)
    def save(self, fileName):
        img.fromarray(imarray).convert('RGB').save(fileName)

"opening fonction"

def open(fileName):
    im = img.open(fileName)
    arr = np.array(im)
    return Image(fileName, arr)
