#! /usr/bin/env python
import argparse

def genElementaryElemsRec(remainingNumbers, current=[]):
    if remainingNumbers:
        for i in range(len(remainingNumbers)):
            c = current.copy()
            c.append(remainingNumbers[i])
            yield from genElementaryElemsRec(remainingNumbers[0:i]+remainingNumbers[i+1:len(remainingNumbers)], c)
    else:
        yield current

def genElementaryElems(n):
    remainingNumbers = list(range(1,n+1))
    yield from genElementaryElemsRec(remainingNumbers)

def valueX(elem):
    return elem[0]

def valueY(elem):
    r = 0
    for i in range(len(elem)):
        if i + 1 == elem[i]:
            r += 1
    return r

def valueXY(elem):
    return valueX(elem), valueY(elem)

def countRandomVariable(elems, var):
    r = {}
    for e in elems:
        val = var(e)
        if val not in r:
            r[val] = 0
        r[val] += 1
    return r

def countXY(elems, varXdist, varYdist):
    dist = []
    dist = {x: {y: 0 for y,v2 in varYdist.items()} for x,v1 in varXdist.items()}
    for e in elems:
        x, y = valueXY(e)
        dist[x][y] += 1
    return dist

def mapDist1DToPropapilities(dist, elems):
    return {i: v / len(elems) for i,v in dist.items()}

def mapDist2DToPropapilities(dist, elems):
    return {x: {y: v / len(elems) for y,v in line.items()} for x,line in dist.items()}

def expectedValue(propabilites):
    r = 0
    for k,v in propabilites.items():
        r+= k*v
    return r

def expectedValueXY(propabilites):
    r = 0
    for kx,line in propabilites.items():
        for ky,v in line.items():
            r+= kx * ky * v
    return r

def covar(propabilites, eX, eY):
    r = 0
    for kx,line in propabilites.items():
        for ky,v in line.items():
            r+= (kx - eX) * (ky - eY) * v
    return r

def printDist(dist):
    WIDTH = 10
    FORMAT_I = "{:10d}"
    FORMAT_F = "{:10.5f}"
    FORMAT_S = "{:10s}"
    first = True
    for x,line in dist.items():
        if first:
            first = False
            h = sorted([k for k,v in line.items()])
            print(" ".join([FORMAT_S.format("x / y")] + [FORMAT_I.format(i) for i in h]))
        print(" ".join([FORMAT_I.format(x)] + [FORMAT_F.format(line[k]) for k in h]))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('N', type=int, help='Number of elements in the tuple')
    args = parser.parse_args()

    elementaryElements = list(genElementaryElems(args.N))
    varXdist = countRandomVariable(elementaryElements, valueX)
    varYdist = countRandomVariable(elementaryElements, valueY)
    disttribution = countXY(elementaryElements, varXdist, varYdist)
    propVarXY = mapDist2DToPropapilities(disttribution, elementaryElements)

    propVarX = mapDist1DToPropapilities(varXdist, elementaryElements)
    propVarY = mapDist1DToPropapilities(varYdist, elementaryElements)
    printDist(propVarXY)

    expectedX = expectedValue(propVarX)
    expectedY = expectedValue(propVarY)
    expectedXY = expectedValueXY(propVarXY)

    covarXY = covar(propVarXY, expectedX, expectedY)

    print("Expected Value X : {:10.5f}".format(expectedX))
    print("Expected Value Y : {:10.5f}".format(expectedY))
    print("Expected Value XY: {:10.5f}".format(expectedXY))
    print("Covariance     XY: {:10.5f}".format(covarXY))


