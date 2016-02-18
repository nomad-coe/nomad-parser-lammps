import fnmatch
import os, sys
import numpy as np
from LAMMPSParserInput import readLoggedThermoOutput

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))

########################################################################################################################
########################################################################################################################
# FIRST I FIND THE LAMMPS OUTPUT LOG FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES

n       = str
extFile = str

if sys.argv[1].endswith("1_methyl_naphthalene"):
    extFile = "naph_298_396_20ns"

if sys.argv[1].endswith("multi"):
    extFile = "-thermo_style_multi"

if sys.argv[1].endswith("one"):
    extFile = "-thermo_style_one"

if sys.argv[1].endswith("custom"):
    extFile = "-thermo_style_custom"

########################################################################################################################
########################################################################################################################
# HERE THE THERMO OUTPUT LOG FILE IS READ AND STORED *IF FOUND*

storedOutput = []
for file in os.listdir(examplesPath):
    if file.endswith(extFile):
        n = file
        storedOutput = open(examplesPath + '/' + n).readlines()    # storing the output log file in the list "storedOutput"

def logFileOpen():  # skip section_frame_sequence if an output log file is not found, i.e., storedOutput = None
    skip = True
    if storedOutput:
        skip = False

    return skip

skip = logFileOpen()

var, thermo_style = readLoggedThermoOutput()   # var = list of output varibles defined when thermo_style = custom
                                               # thermo_style =  output style read in the lammps input file

########################################################################################################################
########################################################################################################################
# NAME CONVENTION (AS IN LAMMPS)

# ke    = Kinetic Energy    (defined in section_frame_sequence)
# pe    = Potential Energy  (defined in section_frame_sequence)
# temp  = Temperature       (defined in section_frame_sequence)
# press = Pressure          (defined in section_frame_sequence)

# LAMMPS thermo_style SETTINGS

# thermo_style = multi   --> multiple-line listing of thermodynamic info that is the equivalent of
#                           'thermo_style custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press'.
#                            The listing contains numeric values and a string ID for each quantity.

# thermo_style = one     --> one-line summary of thermodynamic info that is the equivalent of
#                            'thermo_style custom step temp epair emol etotal press'. The line contains only numeric values.

# thermo_style = custom  --> Style custom is the most general setting and allows you to specify which of
#                            the keywords listed above you want printed on each thermodynamic timestep.
#

########################################################################################################################
########################################################################################################################
# READING THERMO OUTPUTS FOR thermo_style = multi

if thermo_style == 'multi' and skip == False:

    def readFrames():   # reading frames' time step

        frame_filter = filter(lambda x: fnmatch.fnmatch(x, '---------------- Step*'), storedOutput)

        frames = []
        for line in frame_filter:
            line_split = line.split()

            frame = line_split[2]
            frames.append(frame)

            frames_count = len(frames)

        return frames_count


    def readPotEnergy():    # reading frames' potential energy

        pe_filter = filter(lambda x: fnmatch.fnmatch(x, 'PotEng*'), storedOutput)

        pe = []
        for line in pe_filter:
            line_split = line.split()

            pe_fr = float(line_split[2])
            pe.append(pe_fr)
        pe = np.asarray(pe)

        return pe


    def readKinEnergy():    # reading frames' kinetic energy

        ke_filter = filter(lambda x: fnmatch.fnmatch(x, 'TotEng*'), storedOutput)

        ke = []
        temp  = []
        for line in ke_filter:
            line_split = line.split()

            ke_fr = float(line_split[5])
            temp_fr  = float(line_split[8])
            ke.append(ke_fr)
            temp.append(temp_fr)
        ke = np.asarray(ke)
        temp = np.asarray(temp)

        return (ke, temp)


    def readPressure():     # reading frames' pressure

        press_filter = filter(lambda x: fnmatch.fnmatch(x, 'E_coul*'), storedOutput)

        press = []
        for line in press_filter:
            line_split = line.split()

            press_fr = float(line_split[8])
            press.append(press_fr)
        press = np.asarray(press)

        return press


    def readVolume():       # reading frames' simulation box volume

        vol_filter = filter(lambda x: fnmatch.fnmatch(x, 'Volume*'), storedOutput)

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
    for line in storedOutput:
        line = line.strip('\n' + '').split(' ')
        line = filter(None, line)

        # If line is just empty
        if line != []:
            pass
            data.append(line)

    #print data

    # PICK AND STORE LINES STARTING BY FLOAT (THERMO OUTPUT LINES)
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

    # A DICTIONARY IS CREATED WHERE ALL LOGGED VARIABLES ARE STORED (KEYS ARE LAMMPS STANDARD VARIABLE'S NAMES)
    logged_vars = {}
    for i in range(0, n_vars):
        if var[i] != 'step':
            logged_var = { var[i] : [row[i] for row in logged] }  # MATCH VARIABLE NAMES AND VALUES IN A DICT
            logged_vars.update(logged_var)


    # FUNCTION RETURNING THE LISTS STORING THE VARIABLES DEFINED IN section_frame_sequence (NOT THE USER DEFINED ONES)
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
    for line in storedOutput:
        line = line.strip('\n' + '').split(' ')
        line = filter(None, line)

        # If line is just empty
        if line != []:
            pass
            data.append(line)

    #print data

    # PICK AND STORE LINES STARTING BY FLOAT (THERMO OUTPUT LINES)
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

    # A DICTIONARY IS CREATED WHERE ALL LOGGED VARIABLES ARE STORED (KEYS ARE LAMMPS STANDARD VARIABLE'S NAMES)
    logged_vars = {}
    for i in range(0, n_vars):
        if var[i] != 'step':
            logged_var = { var[i] : [row[i] for row in logged] }  # MATCH VARIABLE NAMES AND VALUES IN A DICT
            logged_vars.update(logged_var)


    # FUNCTION RETURNING THE LISTS STORING THE VARIABLES DEFINED IN section_frame_sequence (NOT THE USER DEFINED ONES)
    def pickNOMADVarsOne():

        press = logged_vars.get('press')
        temp  = logged_vars.get('temp')

        return (press, temp)


    press, temp = pickNOMADVarsOne()








