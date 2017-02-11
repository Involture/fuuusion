from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

"open and save"

def opena(file_name) :
    """Return a numpy array of a png image."""
    im = Image.open(file_name)
    return np.array(im)

def _r(array):
    """Return the array of red values of pixels."""
    return array[:,:,0]
def _g(array):
    """Return the array of green values of pixels."""
    return array[:,:,1]
def _b(array):
    """Return the array of blue values of pixels."""
    return array[:,:,2]

def save(imarray, file_name):
    """save an image array"""
    image = Image.fromarray(imarray).convert('RGB')
    image.save(file_name)

"display"

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

def plotisp(spectrums, n):
    """Return a figure containing the 3 plots of the three intensity spectrum
    n is the name of the plot
    """
    fig = plt.figure()
    sub1 = fig.add_subplot(131)
    sub2 = fig.add_subplot(132)
    sub3 = fig.add_subplot(133)
    _plotisp(sub1, spectrums[0], "red", n)
    _plotisp(sub2, spectrums[1], "green", n)
    _plotisp(sub3, spectrums[2], "blue", n)
    return fig

def plotim(imarray):
    """Return a figure containing the image of an image array."""
    fig = plt.figure()
    subp = fig.add_subplot(111)
    subp.imshow(imarray)
    subp.axis('off')
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1)
    return fig

def comp_plotim(imarray, color):
    """Return a figure containing only one color component of the image of an image array."""
    shape = np.shape(imarray)
    arr = np.zeros(shape)
    ind = 0
    if color == "red":
        ind = 0
    elif color == "green":
        ind = 1
    else:
        ind = 2
    arr[:,:,ind] = imarray[:,:,ind]
    return plotim(arr)

"spectrum creation"

def _isp(imarray):
    """Return the intensity spectrum of a single componant array image."""
    spectrum = np.zeros((256,), dtype = np.int)
    for i in imarray:
        spectrum[i] += 1
    return spectrum

def isp(imarray3):
    """Return an array containing the three intensity spectrums of an array image"""
    return np.array([_isp(_r(array)), _isp(_g(array)), _isp(_b(array))]) 

"intensity spectrum transformations"

def _expend(imarray):
    """Expended the spectrum of a single componant image array"""
    start = np.min(imarray)
    end = np.max(imarray)
    width = end - start
    t = (lambda x : ((x - start) / width) * 155)
    imarray = t(imarray).astype(int)

def expend(imarray):
    """Return the spectrum-expended image array"""
    arr = np.copy(imarray)
    for i in range(3):
        _expend(arr[:,:,i])
    return arr
