import wep


c1 = wep.simpleCube()
c2 = wep.simpleCube()
c2.translate(wep.vector(1 / 2,  1 / 2, 1 / 2))
c2.apply(wep.rotateFunction((0.5, 0.5, 0.5), wep.vector(1, 1, 1)))

c3 = c1.intersection(c2)
# print(c3)
c1.union(c2).plot()
c3.plot(True)
