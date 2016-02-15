import fnmatch
import os
import numpy as np

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))
########################################################################################################################
# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
########################################################################################################################
for file in os.listdir(examplesPath):
    if fnmatch.fnmatch(file, '*input*'):
       n = file

lines = open(examplesPath + '/' + n).readlines()

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

lines = map(' '.join, storeInput)

################################################################################################################################

################################
##### LAMMPS INPUT PARSING #####
################################

################################################################################################################################

def readEnsemble():  # HERE I READ THE INTEGRATION TYPE AND POTENTIAL CONSTRAINT ALGORITHM

	ensemble_filter = filter(lambda x: fnmatch.fnmatch(x, 'fix*'), lines)

	for line in ensemble_filter:
		line_split = line.split()

		if line_split[3] == "nve":
			ensemble = "NVE"

		if line_split[3] == "langevin":
			sampling = "langevin_dynamics"

		if line_split[3] == "nvt":
			sampling = "molecular_dynamics"
			ensemble = "NVT"

		if line_split[3] == "npt":
			sampling = "molecular_dynamics"
			ensemble = "NPT"

	return (ensemble, sampling)


def readTPSettings():  # HERE THERMOSTAT/BAROSTAT TARGETS AND RELAXATION TIMES ARE READ

    target_t = 0
    target_p = 0
    thermo_tau = 0
    baro_tau = 0

    ensemble_filter = filter(lambda x: fnmatch.fnmatch(x, 'fix*'), lines)

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

    return (target_t, thermo_tau, target_p, baro_tau)



################################################################################################################################

def readIntegratorSettings():  # HERE I READ INTEGRATOR SETTINGS (TYPE, TIME STEP, NUMBER OF STEPS, ...)

    int_filter = filter(lambda x: fnmatch.fnmatch(x, 'run_style*'), lines)

    if int_filter:
        for line in int_filter:
            line_split = line.split()
            int_type = line_split[1]
    else:
        int_type = "verlet"  # if no run_style command, the integrator is standard Verlet


    run_filt = filter(lambda x: x.startswith("run"), lines)  # OK FOR A SINGLE RUN INPUT SCRIPT

    for line in run_filt:
        line_split = line.split()
        steps = float(line_split[1])


    ts_filter = filter(lambda x: fnmatch.fnmatch(x, 'timestep*'), lines)

    for line in ts_filter:
        line_split = line.split()

        tstep = float(line_split[1])

    return (int_type, tstep, steps)


################################################################################################################################

def readPairCoeff():  # HERE WE COLLECT PAIR COEFFICIENTS (LJ)

    lj_filt = filter(lambda x: x.startswith("pair_coeff"), lines)

    list_of_ljs = {}
    at_types = []
    for line in lj_filt:
        line_split = line.split()
        index1 = line_split[1]
        index2 = line_split[2]

        coeff = []
        for i in range(3, len(line_split)):  # this reads several pair styles coeffs with explicit cutoff (not for ReaxFF)
            if "#" in line_split[i]:
                        break

            try:
                index = float(line_split[i])
                coeff.append(index)
            except ValueError:
                    index = line_split[i]
                    coeff.append(index)

        #index3 = float(line_split[3])
        #index4 = float(line_split[4])

    # creat a list
        lj_coeff = [index1, index2, coeff]
        at_types.append(lj_coeff)

    # create a dictionary

        index1 = index1 + " <=> " + index2

        lj_dict = { index1 : coeff }
        list_of_ljs.update(lj_dict)
    list_of_ljs = { "Pair coefficients" : list_of_ljs}
    return (list_of_ljs, at_types)


########################################################################################################################

def readBonds():   # HERE WE COLLECT BONDS COEFFICIENTS

    bond_filt = filter(lambda x: x.startswith("bond_coeff"), lines)

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

        bond = { 'Force Constant' : index2, 'Equilibrium Length' : index3 }
        bond_dict = {index1 : bond }
        list_of_bonds.update(bond_dict)
    #list_of_bonds = { "Covalent bonds [Force constant, Lenght]" : list_of_bonds }

    return list_of_bonds


########################################################################################################################

def readAngles():

    angle_filt = filter(lambda x: x.startswith("angle_coeff"), lines)

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

        angle = { 'Force Constant' : index2, 'Equilibrium Angle' : index3 }
        angle_dict = {index1 : angle }
        list_of_angles.update(angle_dict)
    #list_of_angles = { "Bond angles [Force constant, Rest angle]" : list_of_angles }
    return list_of_angles


########################################################################################################################

def readDihedrals():

    dihedral_filt = filter(lambda x: x.startswith("dihedral_coeff"), lines)

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

        dihedral = { 'Fourier Coefficients' : [ index2, index3, index4, index5 ] }
        dihedral_dict = {index1 : dihedral }
        list_of_dihedrals.update(dihedral_dict)
    #list_of_dihedrals = { "Dihedral parameters" : list_of_dihedrals }
    return list_of_dihedrals