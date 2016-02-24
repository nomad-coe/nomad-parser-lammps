import fnmatch
import os, sys
import numpy as np

examplesPath = os.path.dirname(os.path.abspath(sys.argv[1]))  # address of the LAMMPS calculation's directory
#examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))

# FIRST I FIND THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
for file in os.listdir(examplesPath):
    if fnmatch.fnmatch(file, '*data.*'):
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


def numberOfTopologyAtoms():

    at_count = 0
    for line in data:
        if "atoms" in line:
            if len(line)==2:
                at_count = int(line[0])  # NUMBER OF ATOMS
    return at_count

at_count = numberOfTopologyAtoms()

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

for line in data:
    if "dihedrals" in line:
        if len(line)==2:
            dh_count = int(line[0])  # NUMBER OF DIHEDRALS


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

dihedral_list = []    #  LIST STORING ALL DIHEDRALS
for i in range(0, len(data)):
    if "Dihedrals" in data[i]:

        for j in range(0, dh_count):
            dh = data[i+j+1]
            dihedral_list.append(dh)
            

########################################################################################################################
def readMass():  # READING ATOMIC MASSES AND CALCULATING ATOMIC NUMBER FOR XYZ RENDERING

    mass_dict = {}
    mass_list = []
    mass_xyz  = []
    for i in range(0, len(data)):
        if "Masses" in data[i]:

            for j in range(0, at_types):
                mass = data[i+j+1]
                index = int(mass[0])
                val   = float(mass[1])
                val1  = int(val/2)     # REQUIRES DOUBLE CHECK


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

mass_dict, mass_list, mass_xyz     = readMass()


########################################################################################################################
def readCharge():  # READING ATOMIC CHARGES

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
def assignBonds():  # ASSIGNING COVALENT BOND TO ITS ATOM PAIR

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
    bondTypeList = []
    for line in bond_pairs:
        ind = int(line[0])
        at1 = int(line[1])
        at2 = int(line[2])
        bondTypeList.append(ind)

        a = topo_list[at1-1][2]
        b = topo_list[at2-1][2]
        bd = { ind : [int(a), int(b)] }
        bond_dict.update(bd)
    #bond_dict = { "Bond assignement (index, at_type1, at_type1)" : bond_dict }

    bondTypeList = sorted(bondTypeList)

    bond_interaction_atoms = []
    for i in bondTypeList:
        for line in bond_list:

            if int(line[1]) == i:
                #count += 1
                hit = [ int(line[1]), int(line[2]), int(line[3]) ]
                bond_interaction_atoms.append(hit)

    return (bond_dict, bondTypeList, bond_interaction_atoms)

bond_dict, bondTypeList, bond_interaction_atoms = assignBonds()

print bondTypeList
#print store_interaction_atoms



########################################################################################################################
def assignAngles():  # ASSIGNING ANGLE TO ITS ATOM TRIPLET

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
    angleTypeList = []
    for line in angle_triplets:
        ind = int(line[0])
        at1 = int(line[1])
        at2 = int(line[2])
        at3 = int(line[3])
        angleTypeList.append(ind)

        a = topo_list[at1-1][2]
        b = topo_list[at2-1][2]
        c = topo_list[at3-1][2]
        ag = { ind : [int(a), int(b), int(c)] }
        angle_dict.update(ag)
    #angle_dict = { "Angle assignement (index, at_type1, at_type1, at_type3)" : angle_dict }

    angleTypeList = sorted(angleTypeList)

    angle_interaction_atoms = []
    for i in angleTypeList:
        for line in angle_list:

            if int(line[1]) == i:
                #count += 1
                hit = [ int(line[1]), int(line[2]), int(line[3]), int(line[4]) ]
                angle_interaction_atoms.append(hit)

    return (angle_dict, angleTypeList, angle_interaction_atoms)

angle_dict, angleTypeList, angle_interaction_atoms = assignAngles()

print angleTypeList



########################################################################################################################
def assignDihedrals():  # ASSIGNING DIHEDRAL TO ITS ATOM QUARTET

    dihedral_ass_d = {}
    dihedral_ass = []
    for line in dihedral_list:
        nr    = line[0]
        index = line[1]
        at1   = line[2]
        at2   = line[3]
        at3   = line[4]
        at4   = line[5]

        store = { nr : [index, at1, at2, at3, at4] }
    
        if index not in dihedral_ass:
            dihedral_ass.append(index)
            dihedral_ass_d.update(store)

    dihedral_quartets = []
    for key, value in dihedral_ass_d.iteritems():
        temp = value
        dihedral_quartets.append(temp)

    dihedral_dict = {}
    dihedralTypeList = []
    for line in dihedral_quartets:
        ind = int(line[0])
        at1 = int(line[1])
        at2 = int(line[2])
        at3 = int(line[3])
        at4 = int(line[4])
        dihedralTypeList.append(ind)
    
        a = topo_list[at1-1][2]
        b = topo_list[at2-1][2]
        c = topo_list[at3-1][2]
        d = topo_list[at4-1][2]
        dh = { ind : [int(a), int(b), int(c), int(d)] }
        dihedral_dict.update(dh)
    #dihedral_dict = { "Dihedral assignement (index, at_type1, at_type1, at_type3, at_type4)" : dihedral_dict }

    dihedralTypeList = sorted(dihedralTypeList)

    dihedral_interaction_atoms = []
    for i in dihedralTypeList:
        for line in dihedral_list:

            if int(line[1]) == i:
                #count += 1
                hit = [ int(line[1]), int(line[2]), int(line[3]), int(line[4]), int(line[5])  ]
                dihedral_interaction_atoms.append(hit)

    return (dihedral_dict, dihedralTypeList, dihedral_interaction_atoms)

dihedral_dict, dihedralTypeList, dihedral_interaction_atoms = assignDihedrals()

print dihedralTypeList
#print dihedral_interaction_atoms


########################################################################################################################
def assignMolecules():  # DEFINING MOLECULE FROM COORDINATION PATTERN

    moleculeId = 0
    store = []
    for i in range(len(bond_list)):
        at1   = int(bond_list[i][2])  # atom 1 index
        at2   = int(bond_list[i][3])  # atom 2 index
        store.append([at1, at2])

    atomIndexInMolecule = []
    for i in range(0,len(store)-1):
        atomIndexInMolecule.append(store[i][0])
        atomIndexInMolecule.append(store[i][1])

        if store[i+1][0] not in atomIndexInMolecule and store[i+1][1] not in atomIndexInMolecule:
            break

    atomIndexInMolecule = list(set(atomIndexInMolecule))  # clear doubles and return the list of atoms in the molecule

    atomTypeInMolecule = []
    for i in atomIndexInMolecule:
        temp = int(topo_list[i-1][2])
        atomTypeInMolecule.append(temp)

    numberOfAtomsInMolecule = len(atomIndexInMolecule)
    numberOfMolecules  = divmod(at_count, len(atomIndexInMolecule))[0]
    mixtureFlag        = divmod(at_count, len(atomIndexInMolecule))[1]

    return (atomIndexInMolecule, atomTypeInMolecule, numberOfAtomsInMolecule, numberOfMolecules)



atomIndexInMolecule, atomTypeInMolecule, numberOfAtomsInMolecule, numberOfMolecules = assignMolecules()
#print atomIndexInMolecule, atomTypeInMolecule, numberOfAtomsInMolecule, numberOfMolecules


########################################################################################################################
xyz_file = []     # WRITE AN XYZ FILE FROM LAMMPS TOPOLOGY DATA
xyz_file.append([at_count])
xyz_file.append([' '])
for line in topo_list:
    index = int(line[2])
    xyz_line = [mass_xyz[index-1], float(line[4]), float(line[5]),  float(line[6])]
    xyz_file.append(xyz_line)

with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), 'generated_from_data_file.xyz')), 'w') as xyz:
    xyz.writelines('  '.join(str(j) for j in i) + '\n' for i in xyz_file)    # WRITE XYZ ATOMIC NUMBER AND COORDINATES