import wep
import math

c1 = wep.simpleCube()
c1.translate(wep.vector(-1 / 2, -1 / 2, -1 / 2))
assert c1.pnFacesInPoly()

sils = []
m = 17

for i in range(m):
    sils.append(c1.silhouettePhoto((10 * math.cos(2 * i * math.pi / m),
                                    10 * math.sin(2 * i * math.pi / m),
                                    10 * math.cos(2 * i * math.pi / m)),
                                   0.30,
                                   (math.pi + i * 2 * math.pi / m,
                                    0,
                                    (-2) * i * math.pi / m)))
    # sils[-1].plot()
assert sils[0].cone(15).pnFacesInPoly()
r = wep.polyhedron.visualHull(sils, 10)
