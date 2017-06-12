import numpy

import wep

l = numpy.load("../analyse_image/res/roller_poly.npy")

l2 = [(a[0], a[1]) for a in l]

l3 = [(l2[i], l2[i+1]) for i in range(len(l2) - 1)] + [(l2[-1], l2[0])]

s = wep.silhouette((0, 0, 0), 0, (0, 0, 0), l3)
s.plot()
print(l3)