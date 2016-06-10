import numpy as np
import math as m
import os, sys

from LAMMPSParserInput import readDumpFileName, simulationTime

stepsPrintFrame = 0
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
    if traj and stepsPrintFrame == stepsPrintThermo: # also checks if the thermo and trajectory output settings are consistent  <=== NOTA BENE
        skipTraj = False                             # FOR NOW, I PARSE THE TRAJ ONLY IF THERMO OUTPUTS AND MD FRAMES ARE PRINTED
                                                     # WITH THE SAME FREQUENCY
    return skipTraj

skipTraj = trajFileOpen()


########################################################################################################################
###### HERE custom LAMMPS TRAJECORY IS PARSED (*.lammpstrj)
########################################################################################################################
if trajDumpStyle == 'custom' and skipTraj == False:

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


        # for i in range(len(trajTotal)):  # atom index and type converted to int
        #     try:
        #         trajTotal[i][:2] = map(int, trajTotal[i][:2])
        #     except ValueError:
        #         pass


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

        pbcBool = []  # periodic boundary conditions (true or false) for x, y and z
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



        return (simulationCell, atomPosition, imageFlagIndex, atomPositionWrapped, atomVelocity, atomForce,
                atomPositionBool, atomPositionBool, imageFlagIndexBool, atomPositionWrappedBool, atomVelocityBool, atomForceBool)



########################################################################################################################
###### HERE atom LAMMPS TRAJECORY IS PARSED (*.atom)
########################################################################################################################
if trajDumpStyle == 'atom' and skipTraj == False:

    nofFrames = integrationSteps/stepsPrintFrame  # number of frames in the trajectory file

    def readAtomTraj():

        # SPLIT AND CLEAN THE TRAJ FILE LINES
        trajTotal = []
        for line in traj:
            line = line.strip('\n' + '').split(' ')
            line = filter(None, line)

            # If line is just empty
            if line != []:
                pass
                trajTotal.append(line)


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

        isAtomPositionScaled  = False
        isAtomPosition        = False
        isImageFlagIndex      = False

        if 'xs' in variables and 'ys' in variables and 'zs' in variables:  # if true, scaled coord are dumped
            isAtomPositionScaled = True

        if 'x' in variables and 'y' in variables and 'z' in variables:     # if true, unwrapped coord are dumped
            isAtomPosition = True

        if 'ix' in variables and 'iy' in variables and 'iz' in variables:  # if true, image flag indexes are dumped
            isImageFlagIndex = True


        # Boundary periodicity bool
        def strToBool(x):
            return x.lower() in ('pp')

        pbcBool = []  # periodic boundary conditions (true or false) for x, y and z
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


        # Atom position (scaled) [0, 1]
        if isAtomPositionScaled == True:

            try:
                xsInd = variables.index('xs')-1
            except ValueError:
                xsInd = 0

            try:
                ysInd = variables.index('ys')-1
            except ValueError:
                ysInd = 0

            try:
                zsInd = variables.index('zs')-1
            except ValueError:
                zsInd = 0

            atomPositionScaledBool = False
            if xsInd and ysInd and zsInd:
                atomPositionScaledBool = True

            store = []
            atomPositionScaled = []
            for frame in frameAtomInfo:
                store_1 =[]
                for line in frame:

                    if xsInd and ysInd and zsInd:
                        store = [ float(line[xsInd]), float(line[ysInd]), float(line[zsInd]) ]
                        store_1.append(store)
                atomPositionScaled.append(store_1)


            # Atom position (converted from scaled positions)
            atomPositionBool = False

            if len(atomPositionScaled)==(nofFrames+1) and atomPositionScaled:
                atomPositionBool = True

                atomPosition = []
                store = []
                i = -1 # frame index
                for frame in atomPositionScaled:
                    i += 1
                    j = -1
                    store_1 = []
                    for atom in frame:

                        j += 1 # atom in frame index
                        x = atom[0] * np.linalg.norm(simulationCell[i][0]) + np.min(simulationCell[i][0])
                        y = atom[1] * np.linalg.norm(simulationCell[i][1]) + np.min(simulationCell[i][1])
                        z = atom[2] * np.linalg.norm(simulationCell[i][2]) + np.min(simulationCell[i][2])
                        store = [x, y, z]
                        store_1.append(store)
                    atomPosition.append(store_1)

        else:
            atomPositionScaled = []
            atomPositionScaledBool = False
            atomPosition = []
            atomPositionBool = False



        # Atom position (unwrapped)
        if isAtomPosition == True:

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

        else:
            atomPosition = []
            atomPositionBool = False


        # Atom position wrapped
        if isAtomPosition == True and isImageFlagIndex == True:

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

            atomPosition = list(atomPosition)
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

        else:
            atomPositionWrapped = []
            atomPositionWrappedBool = False
            imageFlagIndex = []
            imageFlagIndexBool = False


        return (simulationCell, atomPositionScaled, atomPositionScaledBool, atomPosition, atomPositionBool,
                atomPositionWrapped, atomPositionWrappedBool , imageFlagIndex, imageFlagIndexBool )



########################################################################################################################
###### HERE xyz LAMMPS TRAJECORY IS PARSED (*.xyz)
########################################################################################################################
if trajDumpStyle == 'xyz' and skipTraj == False:

    nofFrames = integrationSteps/stepsPrintFrame  # number of frames in the trajectory file

    def readXyzTraj():

        # SPLIT AND CLEAN THE TRAJ FILE LINES
        trajTotal = []
        for line in traj:
            line = line.strip('\n' + '').split(' ')
            line = filter(None, line)

            # If line is just empty
            if line != []:
                pass
                trajTotal.append(line)


        atomPositionBool = False   # just a boolean flag
        if trajTotal:
            atomPositionBool = True


        nofLinesPerFrame = len(trajTotal)/(nofFrames+1)

        trajByFrame = [ trajTotal[i:i + nofLinesPerFrame] for i in xrange(0, len(trajTotal), nofLinesPerFrame) ]  # stepsPrintFrame frame is stored in a list


        frameHeader = []
        frameAtomInfo = []
        for frame in trajByFrame:
            header = frame[:2]
            atomInfo = frame[2:]
            frameHeader.append(header)     # header info per frame  ---> we have just the step number and number of atoms record
            frameAtomInfo.append(atomInfo) # atom info per frame

        for header in frameHeader:         # number of atoms and step number to integer
            for line in header:
                try:
                    line = map(int, line)
                except ValueError:
                    pass


        store = []
        atomPosition = []
        for frame in frameAtomInfo:
            store_1 =[]
            for line in frame:

                store = [ float(line[1]), float(line[2]), float(line[3]) ]
                store_1.append(store)
            atomPosition.append(store_1)


        return (atomPosition, atomPositionBool)
