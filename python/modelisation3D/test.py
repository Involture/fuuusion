import wep


c1 = wep.simpleCube()

sil = c1.silhouettePhoto((-1, 1/2, 1/2), 0.1, (0.5, 0.5, 0))
sil.plot()

sil.cone(5).plot()