from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

"ouverture d'image"
def opena(file_name) :
    im = Image.open(file_name)
    return np.array(im)

def r(array):
    return array[:,:,0]
def g(array):
    return array[:,:,1]
def b(array):
    return array[:,:,2]

"création du spectre"

def sp(array):
    spectrum = np.zeros((255,), dtype = np.int)
    for i in array:
        spectrum[i] += 1
    return spectrum

def sp3(array):
    spectrum = np.zeros((3, 255), dtype = np.int)
    for i in array:
        spectrum[0][i[0]] += 1
        spectrum[1][i[1]] += 1
        spectrum[2][i[2]] += 1
    return spectrum

"affichage du spectre"

def plotsp(fig, spectrum, c, n):
    sf = fig.addsubplot(111)

    N = 256
    ind = np.arange(N)
    width = 1

    sf.bar(ind, spectrum, width, color = c)

    sf.set_title(n)
    sf.set_xlabel("pixel value")
    sf.set_ylabel("nb of pixels")

def plotsp3(fig, spectrum3, n):
    sf = fig.addsubplot(111)

    N = 256 * 3
    ind1 = np.arange(0, N, 3)
    ind2 = np.arange(1, N, 3)
    ind3 = np.arange(2, N, 3)
    width = 1

    sf.bar(ind1, spectrum3[0], width, color = "red")
    sf.bar(ind1, spectrum3[1], width, color = "green")
    sf.bar(ind1, spectrum3[2], width, color = "blue")

    sf.set_title(n)
    sf.set_xlabel("pixel value")
    sf.set_ylabel("nb of pixels")

