import fnmatch
import os
import numpy as np

examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples"))

# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
for file in os.listdir(examplesPath):
    if fnmatch.fnmatch(file, '*data*'):
       n = file

lines = open(examplesPath + '/' + n).readlines()


########################################################################################################################
# SPLIT AND CLEAN THE DATA FILE LINES
data = []
for line in lines:
    line = line.strip('\n' + '').split(' ')
    line = filter(None, line)

        # If line is just empty
    if line != []:
        pass
        data.append(line)


########################################################################################################################
for line in data:
    if "atom" in line:
        if "types" in line:
            at_types = int(line[0])  # NUMBER OF ATOM TYPES

for line in data:
    if "atoms" in line:
        if len(line)==2:
            at_count = int(line[0])  # NUMBER OF ATOMS

for line in data:
    if "bond" in line:
        if "types" in line:
            bd_types = int(line[0])  # NUMBER OF BOND TYPES

for line in data:
    if "bonds" in line:
        if len(line)==2:
            bd_count = int(line[0])  # NUMBER OF BONDS

for line in data:
    if "angles" in line:
        if len(line)==2:
            ag_count = int(line[0])  # NUMBER OF ANGLES


########################################################################################################################
topo_list = []    #  LIST STORING ATOMIC CHARGES AND COORDINATES
for i in range(0, len(data)):
    if "Atoms" in data[i]:

        for j in range(0, at_count):
            topo = data[i+j+1]
            topo_list.append(topo)


bond_list = []    #  LIST STORING ALL BONDS
for i in range(0, len(data)):
    if "Bonds" in data[i]:

        for j in range(0, bd_count):
            bd = data[i+j+1]
            bond_list.append(bd)

angle_list = []    #  LIST STORING ALL ANGLES
for i in range(0, len(data)):
    if "Angles" in data[i]:

        for j in range(0, ag_count):
            ag = data[i+j+1]
            angle_list.append(ag)

########################################################################################################################
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
    return (mass_dict, mass_list, mass_xyz)


########################################################################################################################
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


########################################################################################################################
def assignBonds():  # ASSIGN COVALENT BOND TO ITS ATOM PAIR
    bond_ass_d = {}
    bond_ass = []
    for line in bond_list:
        nr    = line[0]
        index = line[1]
        at1   = line[2]
        at2   = line[3]

        store = { nr : [index, at1, at2] }

        if index not in bond_ass:
            bond_ass.append(index)
            bond_ass_d.update(store)

    bond_pairs = []
    for key, value in bond_ass_d.iteritems():
        temp = value
        bond_pairs.append(temp)

    bond_dict = {}
    for line in bond_pairs:
        ind = line[0]
        at1 = int(line[1])
        at2 = int(line[2])

        a = topo_list[at1-1][2]
        b = topo_list[at2-1][2]
        bd = { ind : [int(a), int(b)] }
        bond_dict.update(bd)
    #bond_dict = { "Bond assignement (index, at_type1, at_type1)" : bond_dict }

    return bond_dict


########################################################################################################################
def assignAngles():  # ASSIGN ANGLE TO ITS ATOM TRIPLET
    angle_ass_d = {}
    angle_ass = []
    for line in angle_list:
        nr    = line[0]
        index = line[1]
        at1   = line[2]
        at2   = line[3]
        at3   = line[4]

        store = { nr : [index, at1, at2, at3] }

        if index not in angle_ass:
            angle_ass.append(index)
            angle_ass_d.update(store)

    angle_triplets = []
    for key, value in angle_ass_d.iteritems():
        temp = value
        angle_triplets.append(temp)

    angle_dict = {}
    for line in angle_triplets:
        ind = int(line[0])
        at1 = int(line[1])
        at2 = int(line[2])
        at3 = int(line[3])

        a = topo_list[at1-1][2]
        b = topo_list[at2-1][2]
        c = topo_list[at3-1][2]
        ag = { ind : [int(a), int(b), int(c)] }
        angle_dict.update(ag)
    #angle_dict = { "Angle assignement (index, at_type1, at_type1, at_type3)" : angle_dict }

    return angle_dict
