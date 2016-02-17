import fnmatch
import os, sys
import numpy as np
from LAMMPSParserInput import readLoggedThermoOutput

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))

# FIRST I FIND THE LAMMPS OUTPUT LOG FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES

n    = str
extFile = str

if sys.argv[1].endswith("1_methyl_naphthalene"):
    extFile = "naph_298_396_20ns"

if sys.argv[1].endswith("multi"):
    extFile = "-thermo_style_multi"

if sys.argv[1].endswith("one"):
    extFile = "-thermo_style_one"

if sys.argv[1].endswith("custom"):
    extFile = "-thermo_style_custom"

lines = []
for file in os.listdir(examplesPath):
    if file.endswith(extFile):
        n = file
        lines = open(examplesPath + '/' + n).readlines()

def logFileOpen():
    skip = True
    if lines:
        skip = False

    return skip

skip = logFileOpen()

#print skip
var, thermo_style = readLoggedThermoOutput()


########################################################################################################################
########################################################################################################################
# READING THERMO OUTPUTS FOR thermo_style = multi

if thermo_style == 'multi' and skip == False:

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
# READING THERMO OUTPUTS FOR thermo_style = multi

if thermo_style == 'custom' and skip == False:

    # SPLIT AND CLEAN THE LOG FILE LINES
    data = []
    for line in lines:
        line = line.strip('\n' + '').split(' ')
        line = filter(None, line)

        # If line is just empty
        if line != []:
            pass
            data.append(line)

    #print data

    # PICK AND STORE LINES STARTING BY FLOAT (DATA LINES)
    logged_tmp = []
    for i in range(0, len(data)):
        x = data[i]
        try:
            float(x[0])
            logged_tmp.append(x)
        except ValueError:
            continue

    logged = []
    for i in range(len(logged_tmp)):
        x = logged_tmp[i]
        try:
            float(x[1])
            logged.append(x)
        except ValueError:
            continue


    for i in range(0, len(logged)):
        logged[i] = map(float, logged[i])

    logged = np.array(logged)

    n_vars = len(var)

    logged_vars = {}
    for i in range(0, n_vars):
        if var[i] != 'step':
            logged_var = { var[i] : [row[i] for row in logged] }  # MATCH VARIABLE NAMES AND VALUES IN A DICT
            logged_vars.update(logged_var)


    def pickNOMADVarsCustom():

        ke    = logged_vars.get('ke')
        pe    = logged_vars.get('pe')
        press = logged_vars.get('press')
        temp  = logged_vars.get('temp')

        return (ke, pe, press, temp)


    ke, pe, press, temp = pickNOMADVarsCustom()


########################################################################################################################
########################################################################################################################
# READING THERMO OUTPUTS FOR thermo_style = one

if thermo_style == 'one' and skip == False:

    # SPLIT AND CLEAN THE LOG FILE LINES
    data = []
    for line in lines:
        line = line.strip('\n' + '').split(' ')
        line = filter(None, line)

        # If line is just empty
        if line != []:
            pass
            data.append(line)

    #print data

    # PICK AND STORE LINES STARTING BY FLOAT (DATA LINES)
    logged_tmp = []
    for i in range(0, len(data)):
        x = data[i]
        try:
            float(x[0])
            logged_tmp.append(x)
        except ValueError:
            continue

    logged = []
    for i in range(len(logged_tmp)):
        x = logged_tmp[i]
        try:
            float(x[1])
            logged.append(x)
        except ValueError:
            continue


    for i in range(0, len(logged)):
        logged[i] = map(float, logged[i])

    logged = np.array(logged)

    n_vars = len(var)

    logged_vars = {}
    for i in range(0, n_vars):
        if var[i] != 'step':
            logged_var = { var[i] : [row[i] for row in logged] }  # MATCH VARIABLE NAMES AND VALUES IN A DICT
            logged_vars.update(logged_var)


    def pickNOMADVarsOne():

        press = logged_vars.get('press')
        temp  = logged_vars.get('temp')

        return (press, temp)


    press, temp = pickNOMADVarsOne()








