import fnmatch
import os, sys
import numpy as np
from LAMMPSParserInput import readLoggedThermoOutput

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))

# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
def logFileStatus():

    n    = str
    skip = bool
    extFile = str

    if sys.argv[1].endswith("multi"):
        extFile = "-thermo_style_multi"

    if sys.argv[1].endswith("one"):
        extFile = "-thermo_style_one"

    if sys.argv[1].endswith("custom"):
        extFile = "-thermo_style_custom"

    for file in os.listdir(examplesPath):

        if file.endswith(extFile):
            n = file
            skip = False

    return (n, skip)


n, skip = logFileStatus()


if skip == False:
    lines = open(examplesPath + '/' + n).readlines()



var, thermo_style = readLoggedThermoOutput()


########################################################################################################################
########################################################################################################################
# READING THERMO OUTPUTS FOR thermo_style = multi

if thermo_style:

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


    def readPressure():

        press_filter = filter(lambda x: fnmatch.fnmatch(x, 'E_coul*'), lines)

        press = []
        for line in press_filter:
            line_split = line.split()

            press_fr = float(line_split[8])
            press.append(press_fr)
        press = np.asarray(press)

        return press


    def readVolume():

        vol_filter = filter(lambda x: fnmatch.fnmatch(x, 'Volume*'), lines)

        if vol_filter:

            vol = []
            for line in vol_filter:
                line_split = line.split()

                vol_fr = float(line_split[8])
                vol.append(vol_fr)
            vol = np.asarray(vol)

            return vol

else:
    pass


########################################################################################################################
########################################################################################################################







