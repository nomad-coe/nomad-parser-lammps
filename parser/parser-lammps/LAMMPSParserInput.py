import fnmatch
import os
import numpy as np


n = '../../test/examples/in.eq_1'                # INPUT FILE NAME

lines = open(n).readlines()

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
