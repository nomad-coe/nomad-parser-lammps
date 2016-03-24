import fnmatch
import os, sys
import numpy as np


#examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))
########################################################################################################################
# FIRST I OPEN THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
########################################################################################################################

n = sys.argv[1]
examplesPath = os.path.dirname(os.path.abspath(n))

lines = open(n).readlines()  # read LAMMPS input file

########################################################################################################################
# HERE I FIND USER DEFINED VARIABLES AND SUBSTITUTE THEIR NUMERIC VALUE THROUGHOUT
# example: variable tempe equal 300.0
########################################################################################################################
storeInput = []
for line in lines:
    line = line.strip('\n' + '').split(' ')
    line = filter(None, line)

    # If line is just empty
    if line != []:
        pass
        storeInput.append(line)


var_name  = []
var_value = []

for line in storeInput:
    if "variable" and "equal" in line:
        name = '$'+line[1]
        try:
            numb = float(line[3])
        except ValueError:
            continue

    else:
        continue

    var_name.append(name)
    var_value.append(numb)

# print var_name, var_value

for i in range(0, len(var_name)):
    storeInput = [[w.replace(var_name[i],str(var_value[i])) for w in line] for line in storeInput]

storeInput = map(' '.join, storeInput)
storeInput = [line for line in storeInput if not line.startswith('#')]  # CLEAR COMMENT LINES

#print storeInput

################################################################################################################################

################################
##### LAMMPS INPUT PARSING #####
################################

################################################################################################################################

def readLogFileName():

    log_filter = filter(lambda x: x.startswith("log "), storeInput)

    if not log_filter:
        logFileName = 'log.lammps'  # thermo output is put in the log.lammps file
                                    # with input echo if a file name is not given from the input

    for line in log_filter:
        line_split = line.split()

        logFileName = line_split[1]

    return logFileName


def readDumpFileName():

    dump_filter = filter(lambda x: x.startswith("dump "), storeInput)

    for line in dump_filter:
        line_split = line.split()

        dumpFileName = line_split[5]

    return dumpFileName

################################################################################################################################

def readStyles():  # HERE WE COLLECT CALCULATIONS STYLES (ATOM, BONDS, ANGLES, DIHEDRALS, ELECTROSTATICS, PAIR STYLES)

    specs_filter = filter(lambda x: fnmatch.fnmatch(x, '*_style*'), lines)
    pm_filter    = filter(lambda x: x.startswith("pair_modify"), lines)
    sb_filter    = filter(lambda x: x.startswith("special_bonds"), lines)

    list_of_styles = {}
    styles_dict    = {}
    for line in specs_filter:
        line_split = line.split()
        if len(line_split)==2:   # the first 2 terms are strings

            index1 = str(line_split[0])
            index2 = str(line_split[1])

    # creat a list
            styles = [index1, index2]

    # create a dictionary
            styles_dict = { index1 : index2 }

        elif len(line_split)>2 and line_split[0] != "thermo_style":  # this reads lj/cut, lj/class2, lj/charmm, lj/gromacs, etc. with cutoff

            index1 = str(line_split[0])
            index2 = str(line_split[1])

            cut = []
            for i in range(2, len(line_split)):
                if "#" in line_split[i]:
                        break

                try:
                    index = float(line_split[i])
                    cut.append(index)
                except ValueError:
                    index = line_split[i]
                    cut.append(index)

            # creat a list
            styles = [index1, index2, cut]

        # create a dictionary
            styles_dict = { index1 : [index2, cut] }
        list_of_styles.update(styles_dict)

    for line in pm_filter:       # reading pair_modify specs
        line_split = line.split()

        modify_dict = { line_split[0] : [ line_split[1:] ] }

    for line in sb_filter:       # reading special_bonds specs
        line_split = line.split()

        special_dict = { line_split[0] : [ line_split[1:] ] }

    list_of_styles.update(modify_dict)
    list_of_styles.update(special_dict)

    return list_of_styles


################################################################################################################################

def readEnsemble():  # HERE I READ THE INTEGRATION TYPE AND POTENTIAL CONSTRAINT ALGORITHM

	ensemble_filter = filter(lambda x: fnmatch.fnmatch(x, 'fix*'), storeInput)

	for line in ensemble_filter:
		line_split = line.split()

		if line_split[3] == "langevin":
			sampling = "langevin_dynamics"

        if line_split[3] == "nve":
            sampling = "molecular_dynamics"
            ensemble = "NVE"

        if line_split[3] == "nvt":
            sampling = "molecular_dynamics"
            ensemble = "NVT"

        if line_split[3] == "npt":
            sampling = "molecular_dynamics"
            ensemble = "NPT"

        if line_split[3] == "nph":
            sampling = "molecular_dynamics"
            ensemble = "NPH"

	return (ensemble, sampling)


def readTPSettings():  # HERE THERMOSTAT/BAROSTAT TARGETS AND RELAXATION TIMES ARE READ

    target_t       = 0
    target_p       = 0
    thermo_tau     = 0
    baro_tau       = 0
    langevin_gamma = 0

    ensemble_filter = filter(lambda x: fnmatch.fnmatch(x, 'fix*'), storeInput)

    for line in ensemble_filter:
        line_split = line.split()

        if line_split[3] == "nvt":
            target_t = float(line_split[6])
            thermo_tau = float(line_split[7])

        if line_split[3] == "npt":
            target_t = float(line_split[6])
            thermo_tau = float(line_split[7])
            target_p = float(line_split[10])
            baro_tau = float(line_split[11])

        if line_split[3] == "nph":
            target_p = float(line_split[6])
            baro_tau = float(line_split[7])

        if line_split[3] == "langevin":
            target_t = float(line_split[5])
            langevin_gamma = float(line_split[6])

    return (target_t, thermo_tau, langevin_gamma, target_p, baro_tau)



################################################################################################################################

def readIntegratorSettings():  # HERE I READ INTEGRATOR SETTINGS (TYPE, TIME STEP, NUMBER OF STEPS, ...)

    int_filter = filter(lambda x: fnmatch.fnmatch(x, 'run_style*'), storeInput)

    if int_filter:
        for line in int_filter:
            line_split = line.split()
            int_type = line_split[1]
    else:
        int_type = "verlet"  # if no run_style command, the integrator is standard Verlet


    run_filt = filter(lambda x: x.startswith("run"), storeInput)  # OK FOR A SINGLE RUN INPUT SCRIPT

    for line in run_filt:
        line_split = line.split()
        steps = int(line_split[1])


    ts_filter = filter(lambda x: fnmatch.fnmatch(x, 'timestep*'), storeInput)

    for line in ts_filter:
        line_split = line.split()

        tstep = float(line_split[1])

    return (int_type, tstep, steps)


int_type, tstep, steps = readIntegratorSettings()

################################################################################################################################

#cacchio = updateAtomicTypes()

#print cacchio

################################################################################################################################

def readPairCoeff():  # HERE WE COLLECT PAIR COEFFICIENTS (LJ)

    lj_filt = filter(lambda x: x.startswith("pair_coeff"), storeInput)

    list_of_ljs = {}
    ljs_dict     = {}
    at_types    = []
    index       = 0
    for line in lj_filt:
        line_split = line.split()

        index += 1
        atom1 = int(line_split[1])
        atom2 = int(line_split[2])

        coeff = []
        for i in range(3, len(line_split)):  # this reads several pair styles coeffs with explicit cutoff (not for ReaxFF)
            if "#" in line_split[i]:
                        break

            try:
                param = float(line_split[i])
                coeff.append(param)
            except ValueError:
                    param = line_split[i]
                    coeff.append(param)

        #index3 = float(line_split[3])
        #index4 = float(line_split[4])

    # creat a list
        lj_coeff = [atom1, atom2, coeff]
        at_types.append(lj_coeff)

    # create dictionaries
        lj_pair = { index : [atom1, atom2] }
        ljs_dict.update(lj_pair)

        lj_param = {index : coeff}
        list_of_ljs.update(lj_param)

    #list_of_ljs = { "Pair coefficients" : list_of_ljs}
    return (list_of_ljs, ljs_dict)


########################################################################################################################

def readBonds():   # HERE WE COLLECT BONDS COEFFICIENTS

    bond_filt = filter(lambda x: x.startswith("bond_coeff"), storeInput)

    list_of_bonds={}
    for line in bond_filt:
        line_split = line.split()
        index1 = int(line_split[1])
        index2 = float(line_split[2])
        index3 = float(line_split[3])

    # creat a list
        bond_coeff = [index1, index2, index3]

    # create a dictionary

    #index1 = 'Bond ' + str(index1)

        #bond = { 'Force Constant' : index2, 'Equilibrium Length' : index3 }
        bond = [ index2, index3 ]
        bond_dict = {index1 : bond }
        list_of_bonds.update(bond_dict)
    #list_of_bonds = { "Covalent bonds [Force constant, Lenght]" : list_of_bonds }

    #print bond_list

    return list_of_bonds


########################################################################################################################

def readAngles():

    angle_filt = filter(lambda x: x.startswith("angle_coeff"), storeInput)

    list_of_angles={}
    for line in angle_filt:
        line_split = line.split()
        index1 = int(line_split[1])
        index2 = float(line_split[2])
        index3 = float(line_split[3])

    # creat a list
        angle_coeff = [index1, index2, index3]

    # create a dictionary

    #index1 = 'Angle ' + str(index1)

        #angle = { 'Force Constant' : index2, 'Equilibrium Angle' : index3 }
        angle = [ index2, index3 ]
        angle_dict = {index1 : angle }
        list_of_angles.update(angle_dict)
    #list_of_angles = { "Bond angles [Force constant, Rest angle]" : list_of_angles }
    return list_of_angles


########################################################################################################################

def readDihedrals():

    dihedral_filt = filter(lambda x: x.startswith("dihedral_coeff"), storeInput)

    list_of_dihedrals={}
    for line in dihedral_filt:
        line_split = line.split()
        index1 = int(line_split[1])
        index2 = float(line_split[2])
        index3 = float(line_split[3])
        index4 = float(line_split[4])
        index5 = float(line_split[5])

    # creat a list
        dihedral_coeff = [index1, index2, index3, index4, index5]

    # create a dictionary

    #index1 = 'Dihedral ' + str(index1)

        #dihedral = { 'Fourier Coefficients' : [ index2, index3, index4, index5 ] }
        dihedral = [ index2, index3, index4, index5 ]
        dihedral_dict = {index1 : dihedral }
        list_of_dihedrals.update(dihedral_dict)
    #list_of_dihedrals = { "Dihedral parameters" : list_of_dihedrals }
    return list_of_dihedrals


################################################################################################################################

def readLoggedThermoOutput():  # HERE I READ THE LIST OF THERMO VARIABLES THAT ARE LOGGED (excludind user-defined ones)

    logvars_filter = filter(lambda x: fnmatch.fnmatch(x, 'thermo_style*'), storeInput)

    var          = []
    thermo_style = str
    for line in logvars_filter:
        line_split = line.split()

        if line_split[1] == "one" or logvars_filter:
            var = ['step', 'temp', 'epair', 'emol', 'etotal', 'press']
            thermo_style = 'one'

        if line_split[1] == "custom":
            var = line_split[2:]
            thermo_style = 'custom'

        if line_split[1] == "multi":
            thermo_style = 'multi'

    var = [i for i in var if not i.startswith('f_') and not i.startswith('v_') and not i.startswith('c_')]

    return (var, thermo_style)


################################################################################################################################

def readUnits():  # HERE I READ THE UNITS STYLE

    units_filter = filter(lambda x: fnmatch.fnmatch(x, '*units*'), storeInput)

    units = []
    for line in units_filter:
            line_split = line.split()
            unitsKey  = line_split[0]
            unitsType = line_split[1]

    units = [unitsKey, unitsType]

    unitsDict = { unitsKey : unitsType }

    return (unitsDict, unitsType)

unitsDict, unitsType = readUnits()


################################################################################################################################

def simulationTime():

    run_filt   = filter(lambda x: x.startswith("run"), storeInput)
    frame_filt = filter(lambda  x: x.startswith("thermo "), storeInput)

    for line in run_filt:
        line_split = line.split()
        steps = float(line_split[1])

        if unitsType == "real":
            time_length = steps * tstep
            #time_length = str(time_length) + " ns"
            #time_length = { "Simulation time" : time_length }

    for line in frame_filt:
        line_split = line.split()
        frame = int(line_split[1])
        frame_length = frame * tstep


    return (frame_length, time_length)




