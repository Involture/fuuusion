from glob import *

def charListToVal(clist):
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

def toArr(string):
    lines = string.split("\n")
    rows = []
    for line in lines:
        valList = charListToVal(list(line))
        rows.append(np.fromiter(valList, dtype = lint))
    return np.row_stack(rows)

def toArrList(l):
    arrList = []
    for subl in l:
        arrList.append(toArr(subl))
    return arrList

def read(dbaseName):
    file = open(dbaseName + ".dbase", "r")
    string0 = file.read()[:-1]
    funcs = string0.split("\n\n")
    for i in range(len(funcs)):
        split = funcs[i].split("\n!\n")
        result = (split[0], toArrList(split[1:]))
        funcs[i] = result
    dic = {}
    for func in funcs:
        dic[func[0]] = func[1]
    return dic
