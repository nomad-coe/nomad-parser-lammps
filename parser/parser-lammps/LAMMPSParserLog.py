import fnmatch
import os
import numpy as np

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))

# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
for file in os.listdir(examplesPath):
    if fnmatch.fnmatch(file, '*log.64*'):
       n = file

lines = open(examplesPath + '/' + n).readlines()


def readFrames():

    frame_filter = filter(lambda x: fnmatch.fnmatch(x, '---------------- Step*'), lines)

    frames = []
    for line in frame_filter:
        line_split = line.split()

        frame = line_split[2]
        frames.append(frame)

        frames_count = len(frames)

    return frames_count


def readPotEnergy():

    poten_filter = filter(lambda x: fnmatch.fnmatch(x, 'PotEng*'), lines)

    poten = []
    for line in poten_filter:
        line_split = line.split()

        poten_fr = float(line_split[2])
        poten.append(poten_fr)
    poten = np.asarray(poten)

    return poten


def readKinEnergy():

    kinen_filter = filter(lambda x: fnmatch.fnmatch(x, 'TotEng*'), lines)

    kinen = []
    temp  = []
    for line in kinen_filter:
        line_split = line.split()

        kinen_fr = float(line_split[5])
        temp_fr  = float(line_split[8])
        kinen.append(kinen_fr)
        temp.append(temp_fr)
    kinen = np.asarray(kinen)
    temp = np.asarray(temp)

    return (kinen, temp)






