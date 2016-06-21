from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import fnmatch
import os, sys, copy, tempfile
from LAMMPSParserInput import readDataFileName, readDumpFileName, readUnits
from LAMMPSParserUnitConversion import unitConversion


fNameData = readDataFileName()                                 ### reading topology data file name from LAMMPS input
fNameTraj, stepsPrintFrame, trajDumpStyle = readDumpFileName() ### reading topology data file name from LAMMPS input

examplesPath = os.path.dirname(os.path.abspath(sys.argv[1]))  # address of the LAMMPS calculation's directory

# FIRST I FIND THE LAMMPS TOPOLOGY DATA FILE
if fNameData:
    lines = open(examplesPath + '/' + fNameData).readlines()

else:
    for file in os.listdir(examplesPath):
        if fnmatch.fnmatch(file, '*data.*'):
            n = file

        lines = open(examplesPath + '/' + n).readlines()


########################################################################################################################
# SPLIT AND CLEAN THE DATA FILE LINES
data = []
for line in lines:
    line = line.strip('\n' + '').split(' ')
    line = [_f for _f in line if _f]

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
    if "angle" in line:
        if "types" in line:
            ag_types = int(line[0])  # NUMBER OF ANGLE TYPES

for line in data:
    if "dihedral" in line:
        if "types" in line:
            dh_types = int(line[0])  # NUMBER OF DIHEDRAL TYPES

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
topo_list.sort(key=lambda x: int(x[0]))  # ordering the atom list (important if the data file is generated from a binary restart file)
# print topo_list

bond_list = []    #  LIST STORING ALL BONDS
for i in range(0, len(data)):
    if "Bonds" in data[i]:

        for j in range(0, bd_count):
            bd = data[i+j+1]
            bond_list.append(bd)
bond_list.sort(key=lambda x: int(x[0]))

angle_list = []    #  LIST STORING ALL ANGLES
for i in range(0, len(data)):
    if "Angles" in data[i]:

        for j in range(0, ag_count):
            ag = data[i+j+1]
            angle_list.append(ag)
angle_list.sort(key=lambda x: int(x[0]))

dihedral_list = []    #  LIST STORING ALL DIHEDRALS
for i in range(0, len(data)):
    if "Dihedrals" in data[i]:

        for j in range(0, dh_count):
            dh = data[i+j+1]
            dihedral_list.append(dh)
dihedral_list.sort(key=lambda x: int(x[0]))
            
#########################################################################################################################################################
# LOAD UNIT CONVERSION

unitsDict, unitsType = readUnits()
toMass,toDistance,toTime,toEnergy,toVelocity,toForce,toTorque,toTemp,toPress,toDynVisc,toCharge,toDipole,toElField,toDensity = unitConversion(unitsType)
toRadians = 0.0174533  # multiply to convert deg to rad

#########################################################################################################################################################

########################################################################################################################
# CREATE A .pdb TOPOLOGY TO BE FED TO MDTraj  (a topology file is required to read binary trajectory files)

with tempfile.NamedTemporaryFile(dir=os.path.dirname('top.pdb')) as pdb:
    for line in topo_list:
        atID = line[0]
        atTy = line[2]
        atCh = line[3]
        atX  = float(line[4])
        atY  = float(line[5])
        atZ  = float(line[6])

        pdb.write('%-6s %4s %2s %5s %1s %3s %11s %7s %7s \n' % ('ATOM', atID, atTy, 'RES', 'X', '1', format(atX, '.3f'), format(atY, '.3f'), format(atZ, '.3f')))
    os.link(pdb.name, 'top.pdb')


########################################################################################################################
# def run_once(f):
#     def wrapper(*args, **kwargs):
#         if not wrapper.has_run:
#             wrapper.has_run = True
#             return f(*args, **kwargs)
#     wrapper.has_run = False
#     return wrapper
#
#
# @run_once

# topo_list_new = copy.deepcopy(topo_list)

def readChargeAndMass():  ### here we record atomic masses and partial charges

    charge_dict   = {}
    charge_list   = []
    mass_dict     = {}
    mass_list     = []
    mass_xyz      = []
    new_mass_list = []

    #topo_list_new = list(topo_list)

    for line in topo_list:
        index = int(line[2])
        charge = float(line[3])
        seen = [index, charge]
        store = { index : charge }

        if seen not in charge_list:
            charge_list.append(seen)
            charge_dict.update(store)

    # print charge_list

    switch = False
    if at_types == len(charge_list):

        switch = True

        for i in range(0, len(data)):
            if "Masses" in data[i]:

                for j in range(0, at_types):
                    mass = data[i+j+1]
                    index = int(mass[0])
                    val   = float(mass[1])
                    val1  = int(val/2)     # REQUIRES DOUBLE CHECK  (calculates atomic number for labelling purposes)

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


        xyz_file = []     # WRITE AN XYZ FILE FROM LAMMPS TOPOLOGY DATA
        atomLabelling = []
        xyz_file.append([at_count])
        xyz_file.append([' '])
        for line in topo_list:
            index = int(line[2])
            xyz_line = [mass_xyz[index-1], float(line[4]), float(line[5]),  float(line[6])]
            xyz_file.append(xyz_line)
            atomLabelling.append(xyz_line)

        with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), 'generated_from_data_file.xyz')), 'w') as xyz:
            xyz.writelines('  '.join(str(j) for j in i) + '\n' for i in xyz_file)    # WRITE XYZ ATOMIC NUMBER AND COORDINATES



        #### A SINGLE ATOM TYPE MIGHT HAVE MORE THAN ONE CHARGE (E.G. CARBON IN CH3, CH2, CH, ...)
        #### WITH THE CODE BELOW, WE CREATE A NEW ATOM TYPE WHENEVER A SINGLE ATOM TYPE HAS MORE THAN ONE ASSOCIATED PARTIAL CHARGE
    if switch == False:

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

        #print mass_list

        new_mass_list = []
        for type in charge_list:
            index = type[0]-1
            # print index
            new_mass_list.append([type[0], mass_list[index][1]])

        #print new_mass_list
        mass_list = list(new_mass_list)

        #########
        for type in mass_list:
            temp = int(type[1]/2)
            mass_xyz.append(temp)

        for i in range(0,len(mass_xyz)):
            if mass_xyz[i] == 0:
                mass_xyz[i] = 1
        #########

        for i in range(len(charge_list)):
            new_mass_list[i] = [ charge_list[i][0], mass_list[i][1] ]
            mass_list[i][0] = i + 1
            mass_dict.update({ mass_list[i][0] : mass_list[i][1] })
            charge_list[i][0] = i + 1
            charge_dict.update({ charge_list[i][0] : charge_list[i][1] })

            #print mass_list

        topo_list_new = []
        for line in topo_list:
            topo_list_new.append(line)

        for type in charge_list:
            for i in range(len(topo_list_new)):

                if type[1] == float(topo_list_new[i][3]):
                    topo_list_new[i][2] = str(type[0])

        topo_list_new.sort(key=lambda x: int(x[0]))

        #print topo_list_new

        xyz_file = []     # WRITE AN XYZ FILE FROM LAMMPS TOPOLOGY DATA
        atomLabelling = []
        xyz_file.append([at_count])
        xyz_file.append([' '])
        for line in topo_list_new:
            index = int(line[2])
            xyz_line = [mass_xyz[index-1], float(line[4]), float(line[5]),  float(line[6])]
            xyz_file.append(xyz_line)
            atomLabelling.append(xyz_line)

        with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), 'generated_from_data_file.xyz')), 'w') as xyz:
            xyz.writelines('  '.join(str(j) for j in i) + '\n' for i in xyz_file)    # WRITE XYZ ATOMIC NUMBER AND COORDINATES


        # for type in charge_list:
        #     for i in range(len(topo_list)):
        #
        #         if type[1] == float(topo_list[i][3]):
        #             topo_list[i][2] = str(type[0])
        #
        # topo_list.sort(key=lambda x: int(x[0]))


    return (charge_dict, charge_list, mass_dict, mass_list, mass_xyz, new_mass_list, atomLabelling)

#print topo_list_new



########################################################################################################################
def readBondsFromData(bondFunctional):

    bond_list = []    #  LIST STORING BOND PARAMETERS
    for i in range(0, len(data)):
        if "Bond" in data[i] and "Coeffs" in data[i]:

            for j in range(0, bd_types):
                bd = data[i+j+1]
                bond_list.append(bd)

    bond_list.sort(key=lambda x: int(x[0]))


    list_of_bonds={}
    for line in bond_list:

        if bondFunctional == "harmonic":
            index1 = int(line[0])
            index2 = float(line[1])*(toEnergy/(toDistance)**2)
            index3 = float(line[2])*toDistance

            bond = [ index2, index3 ]
            bond_dict = {index1 : bond }
            list_of_bonds.update(bond_dict)


        if bondFunctional == "class2":   # COMPASS
            index1 = int(line[0])
            index2 = float(line[1])*toDistance
            index3 = float(line[2])*(toEnergy/(toDistance)**2)
            index4 = float(line[3])*(toEnergy/(toDistance)**3)
            index5 = float(line[4])*(toEnergy/(toDistance)**4)

            bond = [ index2, index3, index4, index5 ]
            bond_dict = {index1 : bond}
            list_of_bonds.update(bond_dict)


        if bondFunctional == "nonlinear":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toDistance
            index4 = float(line[3])*toDistance

            bond = [ index2, index3, index4 ]
            bond_dict = {index1 : bond}
            list_of_bonds.update(bond_dict)


        if bondFunctional == "morse":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*(1/toDistance)
            index4 = float(line[3])*toDistance

            bond = [ index2, index3, index4 ]
            bond_dict = {index1 : bond}
            list_of_bonds.update(bond_dict)


    return list_of_bonds



########################################################################################################################
def readAnglesFromData(angleFunctional):

    angle_list = []    #  LIST STORING ANGLE PARAMETERS
    for i in range(0, len(data)):
        if "Angle" in data[i] and "Coeffs" in data[i]:

            for j in range(0, ag_types):
                ag = data[i+j+1]
                angle_list.append(ag)

    angle_list.sort(key=lambda x: int(x[0]))


    list_of_angles={}
    for line in angle_list:


        if angleFunctional == "harmonic":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toRadians

            angle = [ index2, index3 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "class2":   # COMPASS
            pass


        if angleFunctional == "charmm":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toRadians
            index4 = float(line[3])*(toEnergy/(toDistance)**2)
            index5 = float(line[4])*toDistance

            angle = [ index2, index3, index4, index5 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy

            angle = [ index2 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine/delta":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toRadians

            angle = [ index2, index3 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine/periodic":   # DREIDING
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = int(line[2])
            index4 = int(line[3])

            angle = [ index2, index3, index4 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine/squared":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toRadians

            angle = [ index2, index3 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


    return list_of_angles



########################################################################################################################
def readDihedralsFromData(dihedralFunctional):

    dihedral_list = []    #  LIST STORING DIHEDRAL PARAMETERS
    for i in range(0, len(data)):
        if "Dihedral" in data[i] and "Coeffs" in data[i]:

            for j in range(0, dh_types):
                dh = data[i+j+1]
                dihedral_list.append(dh)

    dihedral_list.sort(key=lambda x: int(x[0]))


    list_of_dihedrals={}
    for line in dihedral_list:


        if dihedralFunctional == "harmonic":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = int(line[2])
            index4 = int(line[3])

            dihedral = [ index2, index3, index4 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "class2":   # COMPASS
            pass


        if dihedralFunctional == "multi/harmonic":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toEnergy
            index4 = float(line[3])*toEnergy
            index5 = float(line[4])*toEnergy
            index6 = float(line[5])*toEnergy

            dihedral = [ index2, index3, index4, index5, index6 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "charmm":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = int(line[2])
            index4 = float(line[3])*toRadians
            index5 = float(line[4])

            dihedral = [ index2, index3, index4, index5 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "opls":   # OPLS aa
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toEnergy
            index4 = float(line[3])*toEnergy
            index5 = float(line[4])*toEnergy

            dihedral = [ index2, index3, index4, index5 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "helix":
            index1 = int(line[0])
            index2 = float(line[1])*toEnergy
            index3 = float(line[2])*toEnergy
            index4 = float(line[3])*toEnergy

            dihedral = [ index2, index3, index4 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


    return list_of_dihedrals



########################################################################################################################
def readPairCoeffFromData(updateAtomTypes, pairFunctional):

    lj_list = []    #  LIST STORING ALL PAIR INTERACTION PARAMETERS
    for i in range(0, len(data)):
        if "Pair" in data[i] and "Coeffs" in data[i]:

            for j in range(0, at_types):
                lj = data[i+j+1]
                lj_list.append(lj)


    lj_list.sort(key=lambda x: int(x[0]))


    # this currently gather pair coefficients for pairFunctional reported in the list of strings below

    supportedLJFunct      = ['lj/cut', 'lj/cut/coul/cut', 'lj/cut/coul/long', 'lj/cut/coul/debye',
                             'lj/cut/coul/dsf', 'lj/cut/coul/msm', 'lj/cut/tip4p/cut', 'lj/cut/tip4p/long']

    supportedCHARMMFunct  = ['lj/charmm/coul/charmm', 'lj/charmm/coul/charmm/implicit', 'lj/charmm/coul/long',
                             'lj/charmm/coul/msm']

    supportedGROMACSFunct = ['lj/gromacs', 'lj/gromacs/coul/gromacs']


    list_of_ljs  = {}
    ljs_dict     = {}
    at_types_lj  = []
    index       = 0
    for line in lj_list:


        if pairFunctional in supportedLJFunct:
            index += 1
            atom1 = int(line[0])   # pair atom type 1
            atom2 = int(line[0])   # pair atom type 2
            eps   = float(line[1])*toEnergy     # epsilon
            sigma = float(line[2])*toDistance   # sigma

            coeff = [eps, sigma]

            for i in range(4, len(line)):  # if another float is present, it is the pair style cutoff(s) for this pair interaction
                if "#" in line[i]:
                            break

                try:
                    rad = float(line[i])*toDistance
                    coeff.append(rad)
                except ValueError:
                        pass

        # creat a list
            lj_coeff = [atom1, atom2, coeff]
            at_types_lj.append(lj_coeff)

        # create dictionaries
            lj_pair = { index : [atom1, atom2] }
            ljs_dict.update(lj_pair)

            lj_param = {index : coeff}
            list_of_ljs.update(lj_param)

        else:
            pass


        if pairFunctional in supportedCHARMMFunct:
            index += 1
            atom1 = int(line[0])   # pair atom type 1
            atom2 = int(line[0])   # pair atom type 2
            eps   = float(line[1])*toEnergy     # epsilon
            sigma = float(line[2])*toDistance   # sigma
            eps14   = float(line[3])*toEnergy     # epsilon 1-4
            sigma14 = float(line[4])*toDistance   # sigma   1-4

            coeff = [eps, sigma, eps14, sigma14]

        # creat a list
            lj_coeff = [atom1, atom2, coeff]
            at_types_lj.append(lj_coeff)

        # create dictionaries
            lj_pair = { index : [atom1, atom2] }
            ljs_dict.update(lj_pair)

            lj_param = {index : coeff}
            list_of_ljs.update(lj_param)

        else:
            pass


        if pairFunctional in supportedGROMACSFunct:
            index += 1
            atom1 = int(line[0])   # pair atom type 1
            atom2 = int(line[0])   # pair atom type 2
            eps   = float(line[1])*toEnergy     # epsilon
            sigma = float(line[2])*toDistance   # sigma
            inner = float(line[3])*toDistance   # inner sigma
            outer = float(line[4])*toDistance   # outer sigma

            coeff = [eps, sigma, inner, outer]

        # creat a list
            lj_coeff = [atom1, atom2, coeff]
            at_types_lj.append(lj_coeff)

        # create dictionaries
            lj_pair = { index : [atom1, atom2] }
            ljs_dict.update(lj_pair)

            lj_param = {index : coeff}
            list_of_ljs.update(lj_param)

        else:
            pass



    if updateAtomTypes:  # here I create pair styles including the new atom types (to account for atoms of the same type, but with different partial charges)

        for line in updateAtomTypes:
            if line[0] != line[1]:

                list_of_ljs.setdefault(line[1], [])
                try:
                    list_of_ljs[line[1]].append(list_of_ljs[line[0]][0])
                except KeyError:
                    pass

                try:
                    list_of_ljs[line[1]].append(list_of_ljs[line[0]][1])
                except KeyError:
                    pass

                ljs_dict.setdefault(line[1], [])
                ljs_dict[line[1]].append(line[1])
                ljs_dict[line[1]].append(line[1])



    return (list_of_ljs, ljs_dict)



########################################################################################################################
def assignBonds(updateAtomTypes):  # ASSIGNING COVALENT BOND TO ITS ATOM PAIR

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
    for key, value in bond_ass_d.items():
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


    ##################
    if updateAtomTypes:

        bond_ass   = []
        for line in bond_list:

            at1   = int(line[2])
            at2   = int(line[3])

            type_at1 = int(topo_list[at1-1][2])
            type_at2 = int(topo_list[at2-1][2])

            store = [ type_at1, type_at2 ]
            store = sorted(store)

            if store not in bond_ass:
                bond_ass.append(store)

            bond_ass = sorted(bond_ass)

        bond_ass_eq = copy.deepcopy(bond_ass)

        for line in updateAtomTypes:
            if line[0] != line[1]:

                for i in range(len(bond_ass_eq)):
                    for j in range(2):
                        if bond_ass_eq[i][j] == line[1]:
                            bond_ass_eq[i][j] = line[0]

        # for line in bond_list:
        #     ind1 = int(line[2])
        #     ind2 = int(line[3])
        #     to_atom_type = sorted( [ int(topo_list[ind1-1][2]), int(topo_list[ind2-1][2]) ] )
        #
        #     new_bond_index = bond_ass.index(to_atom_type) + 1
        #
        #     new_bond_list_line = [ int(line[0]), new_bond_index, ind1, ind2 ]

        new_bond_list =[]
        for line in bond_list:
            bd   = int(line[1])
            ind1 = int(line[2])
            ind2 = int(line[3])
            # to_atom_type = sorted( [ int(topo_list[ind1-1][2]), int(topo_list[ind2-1][2]) ] )

            new_bond_list_line = [ int(line[0]), bd, int(topo_list[ind1-1][2]), int(topo_list[ind2-1][2]) ]
            # print new_bond_list_line
            new_bond_list.append(new_bond_list_line)

        bond_list_unique = []
        for bond in bond_ass:
            for line in new_bond_list:

                temp = [ line[1], [ line[2], line[3] ] ]
                if temp not in bond_list_unique and [temp[0], temp[1][::-1]] not in bond_list_unique:   # VIP
                    bond_list_unique.append(temp)

        updated_bond_list = list(sorted(bond_list_unique))

        bond_dict = {}
        for bond in updated_bond_list:
            bond_dict.setdefault(bond[0], [])
            bond_dict[bond[0]].append(bond[1])

        # for bond in updated_bond_list:
        #     temp = { bond[0] : bond[1] }
        #     new_bond_dict.update(temp)

    return (bond_dict, bondTypeList, bond_interaction_atoms)

# bond_dict, bondTypeList, bond_interaction_atoms = assignBonds(updateAtomTypes)

#print bond_dict, bondTypeList
#print bond_interaction_atoms
#print store_interaction_atoms



########################################################################################################################
def assignAngles(updateAtomTypes):  # ASSIGNING ANGLE TO ITS ATOM TRIPLET

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
    for key, value in angle_ass_d.items():
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


    ##################
    if updateAtomTypes:

        angle_ass   = []
        for line in angle_list:

            at1   = int(line[2])
            at2   = int(line[3])
            at3   = int(line[4])

            type_at1 = int(topo_list[at1-1][2])
            type_at2 = int(topo_list[at2-1][2])
            type_at3 = int(topo_list[at3-1][2])

            store = [ type_at1, type_at2, type_at3 ]
            #store = sorted(store)

            if store not in angle_ass:
                angle_ass.append(store)

        angle_ass = sorted(angle_ass)

        angle_ass_eq = copy.deepcopy(angle_ass)

        for line in updateAtomTypes:
            if line[0] != line[1]:

                for i in range(len(angle_ass_eq)):
                    for j in range(3):
                        if angle_ass_eq[i][j] == line[1]:
                            angle_ass_eq[i][j] = line[0]

        new_angle_list =[]
        for line in angle_list:
            bd   = int(line[1])
            ind1 = int(line[2])
            ind2 = int(line[3])
            ind3 = int(line[4])

            new_angle_list_line = [ int(line[0]), bd, int(topo_list[ind1-1][2]), int(topo_list[ind2-1][2]), int(topo_list[ind3-1][2]) ]
            new_angle_list.append(new_angle_list_line)

        angle_list_unique = []
        for angle in angle_ass:
            for line in new_angle_list:

                temp = [ line[1], [ line[2], line[3], line[4] ] ]
                if temp not in angle_list_unique and [temp[0], temp[1][::-1]] not in angle_list_unique:   # VIP
                    angle_list_unique.append(temp)

        updated_angle_list = list(sorted(angle_list_unique))

        angle_dict = {}
        for angle in updated_angle_list:
            angle_dict.setdefault(angle[0], [])
            angle_dict[angle[0]].append(angle[1])

        # print angle_ass, '########'
        # print angle_dict, '#######'

    return (angle_dict, angleTypeList, angle_interaction_atoms)

#angle_dict, angleTypeList, angle_interaction_atoms = assignAngles()

#print angleTypeList



########################################################################################################################
def assignDihedrals(updateAtomTypes):  # ASSIGNING DIHEDRAL TO ITS ATOM QUARTET

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
    for key, value in dihedral_ass_d.items():
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


    ##################
    if updateAtomTypes:

        dihedral_ass   = []
        for line in dihedral_list:

            at1   = int(line[2])
            at2   = int(line[3])
            at3   = int(line[4])
            at4   = int(line[5])


            type_at1 = int(topo_list[at1-1][2])
            type_at2 = int(topo_list[at2-1][2])
            type_at3 = int(topo_list[at3-1][2])
            type_at4 = int(topo_list[at4-1][2])

            store = [ type_at1, type_at2, type_at3, type_at4 ]
            #store = sorted(store)

            if store not in dihedral_ass:
                dihedral_ass.append(store)

        dihedral_ass = sorted(dihedral_ass)

        dihedral_ass_eq = copy.deepcopy(dihedral_ass)

        for line in updateAtomTypes:
            if line[0] != line[1]:

                for i in range(len(dihedral_ass_eq)):
                    for j in range(4):
                        if dihedral_ass_eq[i][j] == line[1]:
                            dihedral_ass_eq[i][j] = line[0]

        new_dihedral_list =[]
        for line in dihedral_list:
            bd   = int(line[1])
            ind1 = int(line[2])
            ind2 = int(line[3])
            ind3 = int(line[4])
            ind4 = int(line[5])

            new_dihedral_list_line = [ int(line[0]), bd, int(topo_list[ind1-1][2]), int(topo_list[ind2-1][2]), int(topo_list[ind3-1][2]), int(topo_list[ind4-1][2]) ]
            new_dihedral_list.append(new_dihedral_list_line)

        dihedral_list_unique = []
        for dihedral in dihedral_ass:
            for line in new_dihedral_list:

                temp = [ line[1], [ line[2], line[3], line[4], line[5] ] ]
                if temp not in dihedral_list_unique and [temp[0], temp[1][::-1]] not in dihedral_list_unique:   # VIP
                    dihedral_list_unique.append(temp)

        updated_dihedral_list = list(sorted(dihedral_list_unique))

        # print updated_dihedral_list, "#######"

        dihedral_dict = {}
        for dihedral in updated_dihedral_list:
            dihedral_dict.setdefault(dihedral[0], [])
            dihedral_dict[dihedral[0]].append(dihedral[1])

        # print dihedral_ass, '########'
        # print dihedral_dict, '#######'

    return (dihedral_dict, dihedralTypeList, dihedral_interaction_atoms)

#dihedral_dict, dihedralTypeList, dihedral_interaction_atoms = assignDihedrals()

#print dihedralTypeList
#print dihedral_interaction_atoms


########################################################################################################################
def assignMolecules():  # FINDING INDIVIDUAL MOLECULES FROM BONDING PATTERN

    store = []
    for i in range(len(bond_list)):
        at1   = int(bond_list[i][2])  # atom 1 index
        at2   = int(bond_list[i][3])  # atom 2 index
        store.append([at1, at2])


    storeLine = []
    for i in range(len(bond_list)):
        at1   = int(bond_list[i][2])  # atom 1 index
        at2   = int(bond_list[i][3])  # atom 2 index
        storeLine.append(at1)
        storeLine.append(at2)
    #print len(store)

### TOPOLOGY'S ATOM TYPE PATTERN ###
    atomTypeInTopology = []
    for i in range(len(topo_list)):
        temp = int(topo_list[i][2])
        atomTypeInTopology.append(temp) # Atom type pattern throughout topology

    topologyPattern = " ".join([ str(x) for x in atomTypeInTopology])
    #print topologyPattern
############

#####################################################################
###### FIND ALL MOLECULE TYPES AND ALL INDIVIDUAL MOLECULES #########
#####################################################################

    goon             = True # looping over all covalent bonds
    nextMolecule     = 0
    moleculeId       = 0    # molecular index
    moleculePattern  = str  # atom type pattern as string (useful for matching/counting)
    moleculeInfo     = []   # list containing molecule info


    while goon == True:
        atomIndexInMolecule = [] # list storing the index of atoms within that molecule
        atomTypeInMolecule  = [] # list storing the type of atoms within that molecule

        atomPositionInMolecule = [] # list storing the molecular atom index, example: OHH -> 123
        atomCount = 0

        for i in range(nextMolecule, len(store)):
            atomIndexInMolecule.append(store[i][0])
            atomIndexInMolecule.append(store[i][1])

            try:
                if store[i+1][0] not in atomIndexInMolecule and store[i+1][1] not in atomIndexInMolecule:
                    moleculeId += 1
                    #print moleculeId
                    break
            except IndexError:
                moleculeId += 1
                #print moleculeId
                break


        atomIndexToTopology = atomIndexInMolecule
        nextMolecule += len(atomIndexToTopology)/2 ######## NOTA BENE ## I will start from this  to find the next molecule in the list "store"
        #print nextMolecule

        atomIndexInMolecule = sorted(list(set(atomIndexInMolecule)))  # clear duplicates and return the list of atom indexes in the molecule

        for i in atomIndexInMolecule:
            temp = int(topo_list[i-1][2])
            atomTypeInMolecule.append(temp)

            atomCount += 1
            atomPositionInMolecule.append(atomCount)

        moleculePattern = " ".join([ str(x) for x in atomTypeInMolecule])

        newMolecule = [ moleculeId, atomIndexInMolecule, atomTypeInMolecule, atomPositionInMolecule ] # storing molecule information
        moleculeInfo.append(newMolecule)

        #print atomIndexInMolecule
        #print atomTypeInMolecule

        if nextMolecule == len(store):
            goon = False

#### Storing molecule info for each molecule type (the size of moleculeTypeInfo is the number of molecule types)
    moleculeTypeInfo = []
    ghost = []
    for line in moleculeInfo:
        seen = line[2]
        temp = [line[0], line[2], line[1]]

        if seen not in ghost:
            ghost.append(seen)
            moleculeTypeInfo.append(temp)
#################

    for i in range(len(moleculeTypeInfo)):
        for j in range(len(moleculeInfo)):

            if moleculeTypeInfo[i][1] == moleculeInfo[j][2]:
                moleculeInfo[j].insert(1, i+1)   ### moleculeInfo contains: [ moleculeId, moleculeType, atomIndexInMolecule,
                                                 ###                          atomTypeInMolecule, atomPositionInMolecule ]

    atomPositionInMoleculeList = []
    for i in range(len(moleculeInfo)):
        atomPositionInMoleculeList = atomPositionInMoleculeList + moleculeInfo[i][4] # complete list storing the molecular atom index, example: OHHOHHOHH -> 123123123


    moleculeInfoResolved = []
    for i in range(0, at_count):
        for j in range(len(moleculeInfo)):

            if i+1 in moleculeInfo[j][2]:
                moleculeInfoResolved.append([ i+1, moleculeInfo[j][0], moleculeInfo[j][1], atomPositionInMoleculeList[i] ])


    # print moleculeTypeInfo

    return (moleculeTypeInfo, moleculeInfo, moleculeInfoResolved)




########################################################################################################################
