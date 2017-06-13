import math
import pdb
import time
import copy

import wep

c1 = wep.simpleCube()
# c1.translate(wep.vector(-1 / 2, -1 / 2, -1 / 2))
# c1.apply(wep.rotateFunction((0.5, 0.5, 0.5), wep.vector(1, 0.8, 3)))


# t0 = time.time()
# for i in range(100):
#     c1 = wep.simpleCube()
#     c1.translate(wep.vector(-1 / 2, -1 / 2, -1 / 2))
#     c1.apply(wep.rotateFunction((0.5, 0.5, 0.5), wep.vector(1, 0.8, 3)))
#     # c1.plot()
#     c2 = wep.simpleCube()
#     c3 = c2.intersection(c1)
#     # c3.plot()
# print(time.time() - t0)

lims = [0, 0, 0]
lims[0] = [min(v.x for v in c1.vertices),
           max(v.x for v in c1.vertices)]
lims[1] = [min(v.y for v in c1.vertices),
           max(v.y for v in c1.vertices)]
lims[2] = [min(v.z for v in c1.vertices),
           max(v.z for v in c1.vertices)]

assert c1.pnFacesInPoly()
un = copy.deepcopy(c1)
sils = []
m = 5

for i in range(m):
    sils.append(c1.silhouettePhoto((2 * math.cos(2 * i * math.pi / m),
                                    2 * math.sin(2 * i * math.pi / m),
                                    0),
                                   0.10,
                                   (math.pi + i * 2 * math.pi / m,
                                    0,
                                    0)))
    un.union(sils[-1].cone(2))
    # sils[-1].plot()
    # # sils[-1].cone(2).plot()
    # u = sils[-1].cone(1).union(c1)
    # lims = [0, 0, 0]
    # lims[0] = [min(v.x for v in u.vertices),
    #            max(v.x for v in u.vertices)]
    # lims[1] = [min(v.y for v in u.vertices),
    #            max(v.y for v in u.vertices)]
    # lims[2] = [min(v.z for v in u.vertices),
    #            max(v.z for v in u.vertices)]
    # u.plot(ort=True)
    # c1.plot(lims=lims, ort=True)
# un.plot()

assert sils[0].cone(15).pnFacesInPoly()
# pdb.set_trace()
r = wep.polyhedron.visualHull(sils, 10)
print('wait')
r.plot()
