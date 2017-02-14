from imclass import *

def equalize(im):
    arr = np.copy(im.array)
    ispRef = [[],[],[]] #list containing indexes of pixels sorted by intensity value
    if im._ispExist() :
        isp = np.copy(im.isp())
    else:
        isp = np.zeros((N,3))
    #initializing ispRef lists
    for i in range(N):
        for colour in range(3):
            ispRef[colour].append([])
    #initializing ispRef
    for i in range(im.shape[0]):
        for j in range(im.shape[1]):
            for colour in range(3):
                intensity = arr[i,j,colour]
                ispRef[colour][intensity].append((i,j))
                if not im._ispExist():
                    isp[intensity, colour] += 1

    nbPix = im.shape[0] * im.shape[1]
    q = nbPix // N
    r = nbPix % N
    equalRep = np.full((N), q, dtype = np.int64) #desired intensity spectrum
    #intializing equalRep
    for i in range(r):
        equalRep += 1

    newIsp = np.zeros((N,3)) #the isp of the futur equalized image
    for colour in range(3):
        remaining = equalRep[0] 
        newI = 0
        for intensity in range(N):
            currentIPixNb = isp[intensity, colour]
            #loop until desired newI is found
            while  currentIPixNb > remaining:
                currentIPixNb -= remaining
                newI += 1
                remaining = equalRep[newI]
            remaining -= currentIPixNb
            #updating newIsp and changing pixels value to match the desired repartition
            newIsp[newI, colour] += isp[intensity, colour]
            for pix in ispRef[colour][intensity]:
                arr[pix[0], pix[1], colour] = newI

    #setting the argument image isp in case it wasn't im._ispExist()
    if not im._ispExist():
        im._isp = isp
        
    return Image(im.name + '_equalized', arr, newIsp)
