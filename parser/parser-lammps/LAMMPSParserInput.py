import fnmatch
import os
import numpy as np

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples"))

# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
for file in os.listdir(examplesPath):
    if fnmatch.fnmatch(file, '*input*'):
       n = file

lines = open(examplesPath + '/' + n).readlines()

# HERE I FIND USER DEFINED VARIABLES AND SUBSTITUTE THEIR NUMERIC VALUE THROUGHOUT
# example: variable tempe equal 300.0

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
