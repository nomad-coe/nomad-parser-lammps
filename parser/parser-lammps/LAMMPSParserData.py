import fnmatch
import os, sys, copy, tempfile
from LAMMPSParserInput import readDumpFileName

fNameTraj, stepsPrintFrame, trajDumpStyle = readDumpFileName()

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
            

########################################################################################################################
# CREATE A .pdb TOPOLOGY TO BE FED TO MDTraj

# noNumb =  re.compile(r'[#A-Za-z]+')
# topology = []    #  LIST STORING ATOMIC CHARGES AND COORDINATES
# for i in range(at_count):
#
#     topo = topo_list[i]
#     topo_up = []
#     for el in topo:
#
#         if re.match(noNumb, el) is None:
#             exc = el
#             topo_up.append(exc)
#     topology.append(topo_up)


# os.remove('top.pdb')
with tempfile.NamedTemporaryFile(dir=os.path.dirname('top.pdb')) as pdb:
    for line in topo_list:
        atID = line[0]
        atTy = line[2]
        atCh = line[3]
        atX  = float(line[4])
        atY  = float(line[5])
        atZ  = float(line[6])

        pdb.write('%-6s %4s %2s %5s %1s %3s %11s %7s %7s \n' % ('ATOM', atID, atTy, 'RES', 'X', '1', format(atX, '.4f'), format(atY, '.4f'), format(atZ, '.4f')))
    os.link(pdb.name, 'top.pdb')
os.remove('top.pdb')

# mdTrajectory =  md.load(os.path.dirname(os.path.abspath(sys.argv[1])) + '/' + fNameTraj, top='top.pdb')

# mdTopology = md.load_topology('top.pdb')


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

def readChargeAndMass():

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
    #if len(mass_list) != len(charge_list):  #
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

#moleculeInfo, moleculeInfoResolved = assignMolecules()



########################################################################################################################
