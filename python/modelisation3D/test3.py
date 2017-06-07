import wep
import random

for i in range(1000000000000):
    equ = (random.random(),
           random.random(),
           random.random(),
           random.random() * 10)
    v1 = wep.vector(random.random() * 50,
                    random.random() * 50,
                    random.random() * 50)
    v2 = wep.vector(random.random() * 50,
                    random.random() * 50,
                    random.random() * 50)
    vi = wep.vector(wep.planeLineIntersect(v1.coords(), v2.coords(), equ))
    assert ((v1 - v2) * (vi - v2)).norm() < wep.COMPARISON_EPSILON
    print(i)
