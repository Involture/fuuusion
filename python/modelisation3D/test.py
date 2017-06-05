import wep
import math

c1 = wep.simpleCube()

# c2 = wep.simpleCube()
c1.translate(wep.vector(-1 / 2, -1 / 2, -1 / 2))
# c2.apply(wep.rotateFunction((0.5, 0.5, 0.5), wep.vector(1, 1, 1)))
#
# c3 = c1.intersection(c2)
# print(c3)
# c3.plot(True)

sil = c1.silhouettePhoto((-1, 0, 0), 0.1, (0.1, 0, 0))
yaw, roll, pitch = sil.angles
x = wep.vector(1, 0, 0)
y = wep.vector(0, 1, 0)
z = wep.vector(0, 0, 1)
N = math.cos(yaw) * x + math.sin(yaw) * y
A = (N * z).normalize()
Z = math.cos(roll) * z + math.sin(roll) * A
B = (Z * N).normalize()
X = math.cos(pitch) * N + math.sin(pitch) * B
Y = Z * X
# wep.simplex(X, Y, Z).plot()
# sil.plot()
sils = []
m = 8

for i in range(m):
    sils.append(c1.silhouettePhoto((3 * math.cos(2 * math.pi / m),
                                    3 * math.sin(2 * math.pi / m),
                                    0),
                                   0.2,
                                   (i * 2 * math.pi / m,
                                    0,
                                    0)))
    print(sils[i].angles)
    yaw, roll, pitch = sils[i].angles
    x = wep.vector(1, 0, 0)
    y = wep.vector(0, 1, 0)
    z = wep.vector(0, 0, 1)
    N = math.cos(yaw) * x + math.sin(yaw) * y
    A = (N * z).normalize()
    Z = math.cos(roll) * z + math.sin(roll) * A
    B = (Z * N).normalize()
    X = math.cos(pitch) * N + math.sin(pitch) * B
    Y = Z * X
    print(X, Y, Z)
    wep.simplex(X, Y, Z).plot()
# print(len(sil.segments))
# sil.plot()
cone = sil.cone(10)
# cone.plot()
r = wep.polyhedron.visualHull(sils, 10)
