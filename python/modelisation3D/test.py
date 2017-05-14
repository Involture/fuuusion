import wep


c1 = wep.simpleCube()
c2 = wep.simpleCube()
c2.translate(wep.vector(1 / 2, 1 / 2, 1 / 2))
c3 = c1.intersection(c2)
print(c3)
c3.plot(True)
