from glob import *

#geometric tools

def scal(u, v):
    x1, y1 = u
    x2, y2 = v
    return x1 * x2 + y1 * y2

def norm(u):
    x, y = u
    return np.sqrt(x ** 2 + y ** 2)

def normScal(u, v):
    x, y = v
    x2, y2 = x / norm(v), y / norm(v)
    w = (x2, y2)
    a, b = u
    return scal(u, w)

def distFilt(t, v):
    center = (t - 1) / 2
    def vect(i, j):
        x, y = i - center, j - center
        return(x, y)
    return np.fromfunction(lambda i, j : normScal(vect(i, j), v), (t,t))

def circleDist(t):
    hdist = distFilt(t, (1,0))
    vdist = distFilt(t, (0,1))
    return norm((hdist, vdist))

def utheta(theta):
    return (np.cos(theta), np.sin(theta))

def dirGene(ndir):
    return (utheta(k * np.pi / ndir) for k in range(ndir))

#filters and filter generators

def circleCut(t, cut):
    cDist = circleDist(t)
    return cflt(cDist <= cut)

def gauss(t, u, sigma):
    d = distFilt(t, u)
    return np.exp(- (d / sigma) ** 2) / np.sqrt(sigma)

def gaussd(t, u, sigma) :
    d = distFilt(t, u)
    return (-2 * d / sigma ** 2) * np.exp(- (d / sigma) ** 2)

def fGaussGene(winSize, ndir, sigma):
    g = lambda u: gauss(winSize, u, sigma)
    cg = lambda u: g(u) * circleCut(winSize, winSize / 2)
    return (switchQuad(cg(u)) for u in dirGene(ndir))

def gaussdGene(winSize, ndir, sigma):
    g = lambda u: gaussd(winSize, u, sigma)
    cg = lambda u: g(u) * circleCut(winSize, winSize / 2)
    return (cg(u) for u in dirGene(ndir))
