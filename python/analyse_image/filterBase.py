from imclass import * 
im1 = open("image1.png")
im2 = open("image2.png")

def blur(size):
    return np.ones(((2 * size) + 1, (2 * size) + 1), dtype = lint)

def grad():
    grad1 = np.array(
            [[ 1, 1, 1],
             [ 0, 0, 0],
             [-1,-1,-1]], dtype = lint)
    grad2 = np.array(
            [[ 1, 0,-1],
             [ 1, 0,-1],
             [ 1, 0,-1]], dtype = lint)
    grad3 = np.array(
            [[ 1, 1, 0],
             [ 1, 0,-1],
             [ 0,-1,-1]], dtype = lint)
    return [grad1, grad2, grad3]

def lapl():
    lapl1 = np.array(
            [[ 1, 2, 1],
             [ 0, 0, 0],
             [-1,-2,-1]], dtype = lint)
    lapl2 = np.array(
            [[ 2, 1, 0],
             [ 1, 0,-1],
             [ 0,-1,-2]], dtype = lint)
    lapl3 = np.array(
            [[ 1, 0,-1],
             [ 2, 0,-2],
             [ 1, 0,-1]], dtype = lint)
    return [lapl1, lapl2, lapl3]

def squareSet(size):
    return np.full((size, size), True, dtype = np.bool)

def hbarSet(size):
    return np.array([[True] * size, [False] * size], dtype = np.bool)

def vbarSet(size):
    return np.array([[True,False]] * size, dtype = np.bool)
