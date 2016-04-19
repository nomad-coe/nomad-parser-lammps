import numpy as np
import math as m
import os, sys

from LAMMPSParserInput import readDumpFileName, simulationTime

fNameTraj, stepsPrintFrame, trajDumpStyle = readDumpFileName()
frame_length, simulation_length, stepsPrintThermo, integrationSteps = simulationTime()

examplesPath = os.path.dirname(os.path.abspath(sys.argv[1]))  # address of the LAMMPS calculation's directory

########################################################################################################################
########################################################################################################################
# HERE THE TRAJECTORY FILE IS READ AND STORED *IF FOUND*

if fNameTraj:  # if the trajectory file name is defined within the input file, open that file
    n = fNameTraj
    try:
        traj = open(examplesPath + '/' + n).readlines()
    except IOError:
        pass
else:
    traj = []

def trajFileOpen():  # skip section_frame_sequence if an output log file is not found, i.e., storedOutput = None
    skipTraj = True
    if traj:
        skipTraj = False

    return skipTraj

skipTraj = trajFileOpen()
########################################################################################################################

if trajDumpStyle == 'custom' and skipTraj == False and stepsPrintFrame == stepsPrintFrame: # also checks if the thermo and trajectory output settings are consistent  <=== NOTA BENE

    nofFrames = integrationSteps/stepsPrintFrame  # number of frames in the trajectory file

    def readCustomTraj():

        # SPLIT AND CLEAN THE TRAJ FILE LINES
        trajTotal = []
        for line in traj:
            line = line.strip('\n' + '').split(' ')
            line = filter(None, line)

            # If line is just empty
            if line != []:
                pass
                trajTotal.append(line)


        for i in range(len(trajTotal)):  # atom index and type converted to int
            try:
                trajTotal[i][:2] = map(int, trajTotal[i][:2])
            except ValueError:
                pass


        nofLinesPerFrame = len(trajTotal)/(nofFrames+1)

        trajByFrame = [ trajTotal[i:i + nofLinesPerFrame] for i in xrange(0, len(trajTotal), nofLinesPerFrame) ]  # stepsPrintFrame frame is stored in a list


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


        # Calculating simulation cell vectors frame by frame
        simulationCell = []
        for header in frameHeader:
            xx = [ header[5][0] - header[5][1], 0, 0 ]
            yy = [ 0, header[6][0] - header[6][1], 0 ]
            zz = [ 0, 0, header[7][0] - header[7][1] ]
            store = [xx, yy, zz]
            simulationCell.append(store)


        # Atom position (unwrapped)
        try:
            xInd = variables.index('x')-1
        except ValueError:
            xInd = 0

        try:
            yInd = variables.index('y')-1
        except ValueError:
            yInd = 0

        try:
            zInd = variables.index('z')-1
        except ValueError:
            zInd = 0

        atomPositionBool = False
        if xInd and yInd and zInd:
            atomPositionBool = True

        store = []
        atomPosition = []
        for frame in frameAtomInfo:
            store_1 =[]
            for line in frame:

                if xInd and yInd and zInd:
                    store = [ float(line[xInd]), float(line[yInd]), float(line[zInd]) ]
                    store_1.append(store)
            atomPosition.append(store_1)


        # Atom periodic image index
        try:
            ixInd = variables.index('ix')-1
        except ValueError:
            ixInd = 0

        try:
            iyInd = variables.index('iy')-1
        except ValueError:
            iyInd = 0

        try:
            izInd = variables.index('iz')-1
        except ValueError:
            izInd = 0

        imageFlagIndexBool = False
        if ixInd and izInd and izInd:
            imageFlagIndexBool = True

        store = []
        imageFlagIndex = []
        for frame in frameAtomInfo:
            store_1 = []
            for line in frame:

                if ixInd and izInd and izInd:
                    store = [ int(line[ixInd]), int(line[iyInd]), int(line[izInd]) ]
                    store_1.append(store)

            imageFlagIndex.append(store_1)

        # for line in imageFlagIndex[0]:
        #     print line


        # Atom position wrapped
        atomPositionWrappedBool = False

        if len(atomPosition)==(nofFrames+1) and len(atomPosition)==len(imageFlagIndex) and atomPosition and imageFlagIndex:
            atomPositionWrappedBool = True

            atomPositionWrapped = []
            store = []
            i = -1 # frame index
            for frame in atomPosition:
                i += 1
                j = -1
                store_1 = []
                for atom in frame:

                    j += 1 # atom in frame index
                    xw = atom[0] + imageFlagIndex[i][j][0] * np.linalg.norm(simulationCell[i][0])
                    yw = atom[1] + imageFlagIndex[i][j][1] * np.linalg.norm(simulationCell[i][1])
                    zw = atom[2] + imageFlagIndex[i][j][2] * np.linalg.norm(simulationCell[i][2])
                    store = [xw, yw, zw]
                    store_1.append(store)
                atomPositionWrapped.append(store_1)

        # xyz_file = []
        # xyz_file.append(len(atomPositionWrapped[0]))
        # xyz_file.append([' '])
        # for line in atomPositionWrapped[0]:
        #     index = 1
        #     xyz_line = [index, line[0], line[1], line[2]]
        #     xyz_file.append(xyz_line)
        #     print index, line[0], line[1], line[2]

        # with open('generated_from_data_file.xyz', 'w') as xyz:
        #     xyz.writelines('  '.join(str(j) for j in i) + '\n' for i in xyz_file)    # WRITE XYZ ATOMIC NUMBER AND COORDINATES


        # Atom velocities
        try:
            vxInd = variables.index('vx')-1
        except ValueError:
            vxInd = 0

        try:
            vyInd = variables.index('vy')-1
        except ValueError:
            vyInd = 0

        try:
            vzInd = variables.index('vz')-1
        except ValueError:
            vzInd =0

        atomVelocityBool = False
        if vxInd and vyInd and vzInd:
            atomVelocityBool = True

        store = []
        atomVelocity = []
        for frame in frameAtomInfo:
            store_1 = []
            for line in frame:

                if vxInd and vyInd and vzInd:
                    store = [ float(line[vxInd]), float(line[vyInd]), float(line[vzInd]) ]
                    store_1.append(store)

            atomVelocity.append(store_1)

        # for line in atomVelocity:
        #     print line


        # Atom forces
        try:
            fxInd = variables.index('fx')-1
        except ValueError:
            fxInd = 0

        try:
            fyInd = variables.index('fy')-1
        except ValueError:
            fyInd = 0

        try:
            fzInd = variables.index('fz')-1
        except ValueError:
            fzInd = 0

        atomForceBool = False
        if fxInd and fyInd and fzInd:
            atomForceBool = True

        store = []
        atomForce = []
        for frame in frameAtomInfo:
            store_1 = []
            for line in frame:

                if fxInd and fyInd and fzInd:
                    store = [ float(line[fxInd]), float(line[fyInd]), float(line[fzInd]) ]
                    store_1.append(store)

            atomForce.append(store_1)

        # for line in atomVelocity:
        #     print line

        return (simulationCell, atomPosition, imageFlagIndex, atomPositionWrapped, atomVelocity, atomForce,
                atomPositionBool, atomPositionBool, imageFlagIndexBool, atomPositionWrappedBool, atomVelocityBool, atomForceBool)