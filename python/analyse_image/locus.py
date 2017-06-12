import numpy as np

def n(P):
    return np.sqrt(P[0] ** 2 + P[1] ** 2 + P[2] ** 2)

def v(P1, P2):
    return (P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2])

def s(U, V):
    return (U[0] * V[0], U[1] * V[1], U[2] * V[2])

def make_f(A, B, C, k, l, m):
    a = n(C - B)
    b = n(C - A)
    c = n(B - C)

    r = a / 2 * np.sin(k)
    s = r * cos(k)

    ACB_cos = (a ** 2 + b ** 2 + c ** 2)
    ACB_sin = np.sqrt(1 - ACB ** 2)
    h = b * ACB_sin
    d = b * ACB_cos - a / 2

    BB = (-s, a / 2, 0)
    CC = (-s, - a / 2, 0)

    def f(om):
        LL = (r * np.cos(om), r * np.sim(om), 0)
        LB = v(LL, BB)
        LC = v(LL, CC)
        nLB = n(LB)
        nLC = n(LC)



def cameraLocusSolving(A, B, C, k, l, m):

