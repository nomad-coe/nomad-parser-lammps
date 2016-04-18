import numpy as np
import math as m
import os, sys

####
every = 400
total = 80000
nofFrames = total/every
####

fName = sys.argv[1]
traj = open(fName).readlines()

# SPLIT AND CLEAN THE TRAJ FILE LINES
trajTotal = []
for line in traj:
    line = line.strip('\n' + '').split(' ')
    line = filter(None, line)

        # If line is just empty
    if line != []:
        pass
        trajTotal.append(line)


for i in range(len(trajTotal)):
    try:
        trajTotal[i][:2] = map(int, trajTotal[i][:2])
        trajTotal[i][-3:] = map(int, trajTotal[i][-3:])
        trajTotal[i][2:8] = map(float, trajTotal[i][2:8])
    except ValueError:
        pass


nofLinesPerFrame = len(trajTotal)/(nofFrames+1)

trajByFrame = [ trajTotal[i:i + nofLinesPerFrame] for i in xrange(0, len(trajTotal), nofLinesPerFrame) ]  # every frame is stored in a list


frameHeader = []
frameAtomInfo = []
for frame in trajByFrame:
    header = frame[:9]
    atomInfo = frame[9:]
    frameHeader.append(header)     # header info per frame
    frameAtomInfo.append(atomInfo) # atom info per frame

for frame in frameAtomInfo:        # sorting frame atom by index thoughout frameAtomInfo
    frame.sort(key=lambda x: x[0])

for header in frameHeader:         # box boundaries to floating
    for line in header:
        try:
            line[:2] = map(float, line[:2])
        except ValueError:
            pass

# Reading frame column header (stored quantity name)
variables = frameHeader[0][8][1:]
# print variables

# Boundary periodicity bool
def strToBool(x):
  return x.lower() in ('pp')

pbcBool = []
for header in frameHeader:
    for line in header:
        if 'BOX' and 'BOUNDS' in line:
            xbool = strToBool(line[3])
            ybool = strToBool(line[4])
            zbool = strToBool(line[5])
            store = [xbool, ybool, zbool]
            pbcBool.append(store)

# Calculating simulation cell frame by frame
simulationCell = []
for header in frameHeader:
    xx = [ header[5][0] - header[5][1], 0, 0 ]
    yy = [ 0, header[6][0] - header[6][1], 0 ]
    zz = [ 0, 0, header[7][0] - header[7][1] ]
    store = [xx, yy, zz]
    simulationCell.append(store)

# Atom position (unwrapped)
xInd = variables.index('x')-1
yInd = variables.index('y')-1
zInd = variables.index('z')-1
# print len(frameAtomInfo[0][0])
# print frameAtomInfo[0][0]
# print xInd, yInd, zInd

store = []
store_1 =[]
atomPosition = []
for frame in frameAtomInfo:
    for line in frame:

        store = [ line[xInd], line[yInd], line[zInd] ]
        store_1.append(store)
    atomPosition.append(store_1)

print len(atomPosition)
# for line in atomPosition:
#     print line

# Atom periodic image index
ixInd = variables.index('ix')-1
iyInd = variables.index('iy')-1
izInd = variables.index('iz')-1
# print ixInd, iyInd, izInd

store = []
store_1 =[]
imageFlagIndex = []
for frame in frameAtomInfo:
    for line in frame:

        store = [ line[ixInd], line[iyInd], line[izInd] ]
        store_1.append(store)
    imageFlagIndex.append(store_1)

# for line in imageFlagIndex:
#     print line


# Atom position wrapped
i = -1 # frame
for frame in atomPosition:

    i += 1
    j = -1
    for atom in frame:

        j += 1 # atom
        xw = atom[0] + imageFlagIndex[i][j][0] * simulationCell[i][0][0]
        xw = atom[1] + imageFlagIndex[i][j][1] * simulationCell[i][1][1]
        xw = atom[2] + imageFlagIndex[i][j][2] * simulationCell[i][2][2]

        print xw


# Atom velocities
vxInd = variables.index('vx')-1
vyInd = variables.index('vy')-1
vzInd = variables.index('vz')-1

store = []
store_1 =[]
atomVelocity = []
for frame in frameAtomInfo:
    for line in frame:

        store = [ line[vxInd], line[vyInd], line[vzInd] ]
        store_1.append(store)
    atomVelocity.append(store_1)

# for line in atomVelocity:
#     print line


