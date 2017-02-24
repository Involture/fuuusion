from glob import *

def _charListToInt(clist):
    i = 0
    res = []
    while i < len(clist):
        if clist[i] == "-":
            res.append(-int(clist[i+1]))
            i += 2
        else:
            res.append(int(clist[i]))
            i += 1
    return res

def _toArr(string):
    lines = string.split("\n")
    return np.stack([np.fromiter(_charListToInt(list(line)), dtype = lint) for line in lines])

def _toArrList(l):
    arrList = []
    for subl in l:
        arrList.append(_toArr(subl))
    if len(arrList) == 1:
        return arrList[0]
    else:
        return arrList

def read(dbaseName):
    file = open(dbaseName + ".dbase", "r")
    string0 = file.read()
    funcs = string0.split("\n\n")
    for i in range(len(funcs)):
        split = funcs[i].split("\n!\n")
        result = (split[0], _toArrList(split[1:]))
        funcs[i] = result
    dic = {}
    for func in funcs:
        dic[func[0]] = func[1]
    return dic

            
        


