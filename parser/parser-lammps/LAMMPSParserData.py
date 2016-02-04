import fnmatch
import os
import numpy as np

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples"))

# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
for file in os.listdir(examplesPath):
    if fnmatch.fnmatch(file, '*data*'):
       n = file

lines = open(examplesPath + '/' + n).readlines()

# SPLIT AND CLEAN THE DATA FILE LINES
data = []
for line in lines:
    line = line.strip('\n' + '').split(' ')
    line = filter(None, line)

        # If line is just empty
    if line != []:
        pass
        data.append(line)


for line in data:
    if "atom" in line:
        if "types" in line:
            at_types = int(line[0])  # NUMBER OF ATOM TYPES


for line in data:
    if "atoms" in line:
        if len(line)==2:
            at_count = int(line[0])  # NUMBER OF ATOMS


topo_list = []    #  LIST STORING ATOMIC CHARGES AND COORDINATES
for i in range(0, len(data)):
    if "Atoms" in data[i]:

        for j in range(0, at_count):
            topo = data[i+j+1]
            topo_list.append(topo)


def readMass():  # READ ATOMIC MASSES AND CALCULATE ATOMIC NUMBER FOR XYZ RENDERING

    mass_dict = {}
    mass_list = []
    mass_xyz  = []
    for i in range(0, len(data)):
        if "Masses" in data[i]:

            for j in range(0, at_types):
                mass = data[i+j+1]
                index = int(mass[0])
                val   = float(mass[1])
                val1  = int(val/2)


                # create list
                store = [index, val]
                mass_xyz.append(val1)
                mass_list.append(store)

                # create dictionary
                masses = { index : val }
                mass_dict.update(masses)

    for i in range(0,len(mass_xyz)):
        if mass_xyz[i] == 0:
            mass_xyz[i] = 1

    #mass_dict = { "Atomic masses" : mass_dict}
    return (mass_dict, mass_list, mass_xyz, at_types)


def readCharge():  # READ ATOMIC CHARGES

    charge_dict = {}
    charge_list = []
    for line in topo_list:
        index = int(line[2])
        charge = float(line[3])
        seen = [index, charge]
        store = { index : charge }

        if seen not in charge_list:
            charge_list.append(seen)
            charge_dict.update(store)

    #charge_dict = { "Atomic charges" : charge_dict}
    return(charge_dict, charge_list)