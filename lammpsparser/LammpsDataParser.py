# Copyright 2015-2018 Massimo Riello, Adam Fekete, Fawzi Mohamed, Ankit Kariryaa
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import logging
import os
import sys
import operator
import copy

# Python 2 compability
from io import open

from .LammpsCommon import get_metaInfo
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.simple_parser import mainFunction
import nomadcore.ActivateLogging

from re import escape as esc



class sformat(object):
    float = r"[-+]?\d*\.?[0-9]*([Ee][+-]?[0-9]+)?"
    int = r"[-+]?\d+"




############################################################
# This is the parser for the DATA file of LAMMPS.
############################################################

logger = logging.getLogger(name="nomad.LammpsDataParser")





class LammpsDataParserContext(object):

    """Context for parsing LAMMPS DATA file.

    Attributes:
        # dos_energies: Stores parsed energies.
        # dos_values: Stores parsed DOS values.

    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of writing and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """

    def __init__(self, converter, writeMetaData=True):
        """Args:
            writeMetaData: Determines if metadata is written or stored in class attributes.
        """
        self.parser = None
        self.converter = converter
        self.writeMetaData = writeMetaData

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts and the compiled parser is obtained.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser

        # allows to reset values if the same superContext is used to parse different files
        # self.dos_energies = None
        # self.dos_values = None

        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()


    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        self.masses = {}


    # def onClose_section_dos(self, backend, gIndex, section):
    #     """Trigger called when section_dos is closed.
    #
    #     Store the parsed values and write them if writeMetaData is True.
    #     """
    #     # extract energies
    #     dos_energies = section['fhi_aims_dos_energy']
    #     if dos_energies is not None:
    #         pass
    #         # self.dos_energies = np.asarray(dos_energies)
    #     # extract dos values
    #     dos_values = []
    #     if section['fhi_aims_dos_value_string'] is not None:
    #         for string in section['fhi_aims_dos_value_string']:
    #             strings = string.split()
    #             dos_values.append(map(float, strings))
    #     if dos_values:
    #         # need to transpose array since its shape is [number_of_spin_channels,n_dos_values] in the metadata
    #         pass
    #         # self.dos_values = np.transpose(np.asarray(dos_values))
    #     # write metadata only if values were found for both quantities
    #     if self.writeMetaData:
    #         if dos_energies is not None and dos_values:
    #             pass
    #             # backend.addArrayValues('dos_energies', self.dos_energies)
    #             # backend.addArrayValues('dos_values', self.dos_values)

    def onClose_section_run(self, backend, gIndex, section):
        pass

    def onClose_section_topology(self, backend, gIndex, section):

        ####################################################################################################################################################################################################################################
        ####### RETRIEVING TOPOLOGY AND FORCE FIELD INFO FROM ANCILLARY PARSERS ############################################################################################################################################################
        ####################################################################################################################################################################################################################################

        topo_list = []
        if  section['x_lammps_data_topo_list_store']:
            for line in section['x_lammps_data_topo_list_store']:
                topo_list.append(line.split())

            topo_list.sort(key=lambda x: int(x[0]))  # ordering the atom list (important if the data file is generated from a binary restart file)

        # collecting atomic masses
        charge_dict, charge_list, mass_dict, mass_list, mass_xyz, new_mass_list, atomLabelling  = self.readChargeAndMass(topo_list, section)
        mass_dict = sorted(list(mass_dict.items()), key=operator.itemgetter(0))
        charge_dict = sorted(list(charge_dict.items()), key=operator.itemgetter(0)) # ordering atomic partial charges

        # return values of the original function
        self.charge_dict = charge_dict
        self.charge_list = charge_list
        self.mass_dict = mass_dict
        self.mass_list = mass_list
        self.mass_xyz = mass_xyz
        self.new_mass_list = new_mass_list
        self.atomLabelling = atomLabelling


        updateAtomTypes = []  ### to account for atoms with the same type in the topology, but with different partial charges
        ### see function readChargeAndMass in LAMMPSParserData
        if new_mass_list:
            for i in range(len(mass_list)):
                updateAtomTypes.append([ new_mass_list[i][0], mass_list[i][0] ])
        else:
            pass
            ###

        self.updateAtomTypes = updateAtomTypes


        # collecting information defining different molecules
        # LEGEND:
        # moleculeTypeInfo     --> for each molecule type we store [ molecule type index, atom type each atom in the molecule,
        #                                                            absolute atom index for the first occurence of this molecule type in the topology ]
        # moleculeInfo         --> for each molecule we store      [ absoulute molecule index, molecule type, absolute molecular atom index within the topology,
        #                                                            molecular atom index within the molecule ]
        # moleculeInfoResolved --> for each topology atom we store [ absolute atom index within the topology, molecule index,
        #                                                            molecule type index, molecular atom index within the molecule ]


        bond_list = []    #  LIST STORING ALL BONDS
        if section['x_lammps_data_bond_list_store']:
            for line in section['x_lammps_data_bond_list_store']:
                bond_list.append(line.split())

            bond_list.sort(key=lambda x: int(x[0]))


        updateAtomTypes = []  ### to account for atoms with the same type in the topology, but with different partial charges
        ### see function readChargeAndMass in LAMMPSParserData
        if new_mass_list:
            for i in range(len(mass_list)):
                updateAtomTypes.append([ new_mass_list[i][0], mass_list[i][0] ])
        else:
            pass
            ###


        # collecting information defining different molecules
        # LEGEND:
        # moleculeTypeInfo     --> for each molecule type we store [ molecule type index, atom type each atom in the molecule,
        #                                                            absolute atom index for the first occurence of this molecule type in the topology ]
        # moleculeInfo         --> for each molecule we store      [ absoulute molecule index, molecule type, absolute molecular atom index within the topology,
        #                                                            molecular atom index within the molecule ]
        # moleculeInfoResolved --> for each topology atom we store [ absolute atom index within the topology, molecule index,
        #                                                            molecule type index, molecular atom index within the molecule ]

        moleculeTypeInfo, moleculeInfo, moleculeInfoResolved = self.assignMolecules(topo_list, bond_list, section)

        # return values of assignMolecules
        self.moleculeTypeInfo = moleculeTypeInfo
        self.moleculeInfo = moleculeInfo
        self.moleculeInfoResolved = moleculeInfoResolved
        ###

        # # collecting list of force field functional styles
        # list_of_styles = readStyles()
        # pairFunctionalAndCutOff = list_of_styles.get('pair_style')  ### pair interactions style LAMMPS name + cutoff
        # pairFunctional = pairFunctionalAndCutOff[0]                 ### pair interactions style LAMMPS name only
        # bondFunctional = list_of_styles.get('bond_style')
        # angleFunctional = list_of_styles.get('angle_style')
        # dihedralFunctional = list_of_styles.get('dihedral_style')
        ###

        # collecting covalent bond definitions
        bond_dict, bondTypeList, bond_interaction_atoms  = self.assignBonds(topo_list, bond_list, updateAtomTypes, section)
        bond_dict = sorted(list(bond_dict.items()), key=operator.itemgetter(0))
        bd_types = len(bond_dict)

        self.bond_list = bond_list
        self.bond_dict = bond_dict
        self.bondTypeList = bondTypeList
        self.bond_interaction_atoms = bond_interaction_atoms
        self.bd_types = bd_types
        ###

        angle_list = []    #  LIST STORING ALL ANGLES
        if section['x_lammps_data_angle_list_store']:
            for line in section['x_lammps_data_angle_list_store']:
                angle_list.append(line.split())

            angle_list.sort(key=lambda x: int(x[0]))

        # collecting bond angles definitions
        angle_dict, angleTypeList, angle_interaction_atoms  = self.assignAngles(angle_list, topo_list, updateAtomTypes)
        angle_dict = sorted(list(angle_dict.items()), key=operator.itemgetter(0))
        ag_types = len(angle_dict)
        ###

        self.angle_list = angle_list
        self.angle_dict = angle_dict
        self.angleTypeList = angleTypeList
        self.angle_interaction_atoms = angle_interaction_atoms
        self.ag_types = ag_types

        dihedral_list = []    #  LIST STORING ALL DIHEDRALS
        if section['x_lammps_data_dihedral_list_store']:
            for line in section['x_lammps_data_dihedral_list_store']:
                dihedral_list.append(line.split())

            dihedral_list.sort(key=lambda x: int(x[0]))


        # collecting dihedral angles definitions
        dihedral_dict, dihedralTypeList, dihedral_interaction_atoms  = self.assignDihedrals(topo_list, dihedral_list, updateAtomTypes)
        dihedral_dict = sorted(list(dihedral_dict.items()), key=operator.itemgetter(0))
        dh_types = len(dihedral_dict)
        ###


        self.dihedral_list = dihedral_list
        self.dihedral_dict = dihedral_dict
        self.dihedralTypeList = dihedralTypeList
        self.dihedral_interaction_atoms = dihedral_interaction_atoms
        self.dh_types = dh_types

        dihedral_coeff_list = []
        if section['x_lammps_data_dihedral_coeff_list_store']:
            for line in section['x_lammps_data_dihedral_coeff_list_store']:
                dihedral_list.append(line.split())

            dihedral_list.sort(key=lambda x: int(x[0]))


        # # NUMBER OF BONDS (bd_count)
        # SM(r"\s(?P<x_lammps_data_bd_count_store>{0})\sbond".format(sformat.int)),
        #
        # # NUMBER OF ANGLES (ag_count)
        # SM(r"\s(?P<x_lammps_data_ag_count_store>{0})\sangles".format(sformat.int)),
        #
        # # NUMBER OF DIHEDRALS (dh_count)
        # SM(r"\s(?P<x_lammps_data_dh_count_store>{0})\sdihedrals".format(sformat.int)),
        # SM(r"\s(?P<x_lammps_dummy>{0})\simpropers".format(sformat.int)),
        #
        #
        #
        # # NUMBER OF BOND TYPES
        # SM(r"\s(?P<x_lammps_data_bd_types_store>{0})\sbond types".format(sformat.int)),
        #
        # # NUMBER OF ANGLE TYPES (ag_types)
        # SM(r"\s(?P<x_lammps_dummy>{0})\sangle types".format(sformat.int)),
        #
        # # NUMBER OF DIHEDRAL TYPES (dh_types)
        pass

    ########################################################################################################################
    # moleculeTypeInfo, moleculeInfo, moleculeInfoResolved = self.assignMolecules(backend, gIndex, section)
    ########################################################################################################################
    def assignMolecules(self,topo_list, bond_list, section):  # FINDING INDIVIDUAL MOLECULES FROM BONDING PATTERN

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
            nextMolecule += len(atomIndexToTopology) // 2 ######## NOTA BENE ## I will start from this  to find the next molecule in the list "store"
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

        # NUMBER OF ATOMS (at_count)
        at_count = int(section['number_of_topology_atoms'][0])
        self.at_count = at_count

        moleculeInfoResolved = []
        for i in range(0, at_count):
            for j in range(len(moleculeInfo)):

                if i+1 in moleculeInfo[j][2]:
                    moleculeInfoResolved.append([ i+1, moleculeInfo[j][0], moleculeInfo[j][1], atomPositionInMoleculeList[i] ])


        # print moleculeTypeInfo
        return (moleculeTypeInfo, moleculeInfo, moleculeInfoResolved)



        # collecting covalent bond definitions
        # bond_dict, bondTypeList, bond_interaction_atoms  = assignBonds(updateAtomTypes)
        # def assignBonds(updateAtomTypes):  # ASSIGNING COVALENT BOND TO ITS ATOM PAIR

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

        # return (bond_dict, bondTypeList, bond_interaction_atoms)

        # bond_dict = sorted(list(bond_dict.items()), key=operator.itemgetter(0))
        #
        # list_of_bonds = readBonds(bondFunctional)
        # list_of_bonds = sorted(list(list_of_bonds.items()), key=operator.itemgetter(0))
        # bd_types = len(bond_dict)
        ###

    def readChargeAndMass(self, topo_list, section):  ### here we record atomic masses and partial charges
        charge_dict   = {}
        charge_list   = []
        mass_dict     = {}
        mass_list     = []
        mass_xyz      = []
        new_mass_list = []


        for line in topo_list:
            index = int(line[2])
            charge = float(line[3])
            seen = [index, charge]
            store = { index : charge }

            if seen not in charge_list:
                charge_list.append(seen)
                charge_dict.update(store)


        # NUMBER OF ATOMS (at_count)
        at_count = int(section['number_of_topology_atoms'][0])

        # NUMBER OF ATOM TYPES
        at_types = section['x_lammps_data_at_types_store'][0]


        switch = False
        if at_types == len(charge_list):

            switch = True

            masses = section['x_lammps_masses_store']
            for data in masses:
                mass = data.split()
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

            # mode = 'w'
            # if sys.version_info.major < 3:
            #     mode += 'b'
            #
            # with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), 'generated_from_data_file.xyz')), mode) as xyz:
            #     xyz.writelines('  '.join(str(j) for j in i) + '\n' for i in xyz_file)    # WRITE XYZ ATOMIC NUMBER AND COORDINATES


        #### A SINGLE ATOM TYPE MIGHT HAVE MORE THAN ONE CHARGE (E.G. CARBON IN CH3, CH2, CH, ...)
        #### WITH THE CODE BELOW, WE CREATE A NEW ATOM TYPE WHENEVER A SINGLE ATOM TYPE HAS MORE THAN ONE ASSOCIATED PARTIAL CHARGE
        if switch == False:


            masses = section['x_lammps_masses_store']
            for data in masses:
                mass = data.split()

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

            # mode = 'w'
            # if sys.version_info.major < 3:
            #     mode += 'b'
            #
            # with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), 'generated_from_data_file.xyz')), mode) as xyz:
            #     xyz.writelines('  '.join(str(j) for j in i) + '\n' for i in xyz_file)    # WRITE XYZ ATOMIC NUMBER AND COORDINATES

        return charge_dict, charge_list, mass_dict, mass_list, mass_xyz, new_mass_list, atomLabelling


    ########################################################################################################################
    def assignBonds(self, topo_list, bond_list, updateAtomTypes, section):  # ASSIGNING COVALENT BOND TO ITS ATOM PAIR

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

    ########################################################################################################################
    def assignAngles(self, angle_list, topo_list, updateAtomTypes):  # ASSIGNING ANGLE TO ITS ATOM TRIPLET

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

        ########################################################################################################################
    def assignDihedrals(self, topo_list, dihedral_list, updateAtomTypes):  # ASSIGNING DIHEDRAL TO ITS ATOM QUARTET

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
    def readBondsFromData(self, bond_list, bondFunctional):

        toMass     = self.converter.ratioMass
        toDistance = self.converter.ratioDistance
        toTime     = self.converter.ratioTime
        toEnergy   = self.converter.ratioEnergy
        toVelocity = self.converter.ratioVelocity
        toForce    = self.converter.ratioForce
        toTorque   = self.converter.ratioTorque
        toTemp     = self.converter.ratioTemp
        toPress    = self.converter.ratioPress
        toDVisc    = self.converter.ratioDVisc
        toCharge   = self.converter.ratioCharge
        toDipole   = self.converter.ratioDipole
        toEField   = self.converter.ratioEField
        toDensity  = self.converter.ratioDensity


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
    def readAnglesFromData(self, angle_list, angleFunctional):

        toMass     = self.converter.ratioMass
        toDistance = self.converter.ratioDistance
        toTime     = self.converter.ratioTime
        toEnergy   = self.converter.ratioEnergy
        toVelocity = self.converter.ratioVelocity
        toForce    = self.converter.ratioForce
        toTorque   = self.converter.ratioTorque
        toTemp     = self.converter.ratioTemp
        toPress    = self.converter.ratioPress
        toDVisc    = self.converter.ratioDVisc
        toCharge   = self.converter.ratioCharge
        toDipole   = self.converter.ratioDipole
        toEField   = self.converter.ratioEField
        toDensity  = self.converter.ratioDensity

        toRadians = 0.0174533  # multiply to convert deg to rad


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
    def readDihedralsFromData(self, dihedral_list, dihedralFunctional):

        toMass     = self.converter.ratioMass
        toDistance = self.converter.ratioDistance
        toTime     = self.converter.ratioTime
        toEnergy   = self.converter.ratioEnergy
        toVelocity = self.converter.ratioVelocity
        toForce    = self.converter.ratioForce
        toTorque   = self.converter.ratioTorque
        toTemp     = self.converter.ratioTemp
        toPress    = self.converter.ratioPress
        toDVisc    = self.converter.ratioDVisc
        toCharge   = self.converter.ratioCharge
        toDipole   = self.converter.ratioDipole
        toEField   = self.converter.ratioEField
        toDensity  = self.converter.ratioDensity

        toRadians = 0.0174533  # multiply to convert deg to rad

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
    def readPairCoeffFromData(self, lj_list, updateAtomTypes, pairFunctional):

        toMass     = self.converter.ratioMass
        toDistance = self.converter.ratioDistance
        toTime     = self.converter.ratioTime
        toEnergy   = self.converter.ratioEnergy
        toVelocity = self.converter.ratioVelocity
        toForce    = self.converter.ratioForce
        toTorque   = self.converter.ratioTorque
        toTemp     = self.converter.ratioTemp
        toPress    = self.converter.ratioPress
        toDVisc    = self.converter.ratioDVisc
        toCharge   = self.converter.ratioCharge
        toDipole   = self.converter.ratioDipole
        toEField   = self.converter.ratioEField
        toDensity  = self.converter.ratioDensity


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





def build_LammpsDataFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the DATA file of LAMMPS.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    HEADERS = set([
        //'atoms',
        //'bonds',
        //'angles',
        //'dihedrals',
        //'impropers',
        //'atom types',
        //'bond types',
        //'angle types',
        //'dihedral types',
        //'improper types',
        //'extra bond per atom',
        //'extra angle per atom',
        //'extra dihedral per atom',
        //'extra improper per atom',
        //'extra special per atom',
        //'ellipsoids',
        //'lines',
        //'triangles',
        //'bodies',
        //'xlo xhi',
        //'ylo yhi',
        //'zlo zhi',
        //'xy xz yz',
    ])


    # Sections will all start with one of these words
    # and run until the next section title
    SECTIONS = set([
        //'Atoms',  # Molecular topology sections
        //'Velocities',
        //'Masses',
        //'Ellipsoids',
        //'Lines',
        //'Triangles',
        //'Bodies',
        //'Bonds',  # Forcefield sections
        //'Angles',
        //'Dihedrals',
        //'Impropers',
        //'Pair',
        //'Pair LJCoeffs',
        //'Bond Coeffs',
        //'Angle Coeffs',
        //'Dihedral Coeffs',
        //'Improper Coeffs',
        //'BondBond Coeffs',  # Class 2 FF sections
        //'BondAngle Coeffs',
        //'MiddleBondTorsion Coeffs',
        //'EndBondTorsion Coeffs',
        //'AngleTorsion Coeffs',
        //'AngleAngleTorsion Coeffs',
        //'BondBond13 Coeffs',
        //'AngleAngle Coeffs',
    ])


    Returns:
       SimpleMatcher that parses DATA file of FHI-LAMMPS.
    """
    atoms = SM(
        name='lammps-data-list-of-atoms',
        startReStr=r"\s*Atoms",
        repeats=True,
        subMatchers=[
            SM(
                r"\s*(?P<x_lammps_data_topo_list_store>{0}\s{0}\s{0}\s{1}\s{1}\s{1}\s{1})".format(sformat.int, sformat.float),
                repeats=True
            )
        ])

    velocities = SM(
        name='lammps-data-list-of-velocities',
        startReStr=r"\s*Velocities",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_topo_list_store>{0}\s{0}\s{0}\s{1}\s{1}\s{1}\s{1})".format(sformat.int, sformat.float),
            #     repeats=True
            # )
        ])

    masses = SM(
        name='lammps-data-list-of-masses',
        startReStr=r"\s*Masses",
        repeats=True,
        subMatchers=[
            SM(
                r"\s*(?P<x_lammps_masses_store>{0}\s{1})".format(sformat.int, sformat.float),
                repeats=True,
            )
        ])

    ellipsoids = SM(
        name='lammps-data-list-of-ellipsoids',
        startReStr=r"\s*Ellipsoids",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_topo_list_store>{0}\s{0}\s{0}\s{1}\s{1}\s{1}\s{1})".format(sformat.int, sformat.float),
            #     repeats=True
            # )
        ])



    lines = SM(
        name='lammps-data-list-of-lines',
        startReStr=r"\s*Lines",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_topo_list_store>{0}\s{0}\s{0}\s{1}\s{1}\s{1}\s{1})".format(sformat.int, sformat.float),
            #     repeats=True
            # )
        ])

    triangles = SM(
        name='lammps-data-list-of-triangles',
        startReStr=r"\s*Triangles",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_topo_list_store>{0}\s{0}\s{0}\s{1}\s{1}\s{1}\s{1})".format(sformat.int, sformat.float),
            #     repeats=True
            # )
        ])

    bodies = SM(
        name='lammps-data-list-of-bodies',
        startReStr=r"\s*Bodies",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_topo_list_store>{0}\s{0}\s{0}\s{1}\s{1}\s{1}\s{1})".format(sformat.int, sformat.float),
            #     repeats=True
            # )
        ])

    bonds = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Bonds",
        repeats=True,
        subMatchers=[
           SM(
               r"\s*(?P<x_lammps_data_bond_list_store>{0}\s{0}\s{0}\s{0})".format(sformat.int),
               repeats=True
           )
        ])

    angles = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Angles",
        repeats=True,
        subMatchers=[
            SM(
                r"\s*(?P<x_lammps_data_angle_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
                repeats=True
            )
        ])

    dihedrals = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Dihedrals",
        repeats=True,
        subMatchers=[
            SM(
                r"\s*(?P<x_lammps_data_dihedral_list_store>{0}\s{0}\s{0}\s{0}(\s{0})+)".format(sformat.int),
                repeats=True
            )
        ])






    impropers = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Impropers",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    pair = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Pair",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    pair_ljcoeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Pair LJCoeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    bond_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Bond Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    angle_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Angle Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    dihedral_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Dihedral Coeffs",
        repeats=True,
        subMatchers=[
            SM(
                r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
                repeats=True
            )
        ])

    improper_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*Improper Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    bondbond_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*BondBond Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    bondangle_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*BondAngle Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    middlebondtorsion_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*MiddleBondTorsion Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    endbondtorsion_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*EndBondTorsion Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    angletorsion_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*AngleTorsion Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    angleangletorsion_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*AngleAngleTorsion Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    bondbond13_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*BondBond13 Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])
    angleangle_coeffs = SM(
        name='lammps-data-list-of-bonds',
        startReStr=r"\s*AngleAngle Coeffs",
        repeats=True,
        subMatchers=[
            # SM(
            #     r"\s*(?P<x_lammps_data_dihedral_coeff_list_store>{0}\s{0}\s{0}\s{0}\s{0})".format(sformat.int),
            #     repeats=True
            # )
        ])





    return SM(
        name='root1',
        startReStr="",
        sections=['section_run'],
        forwardMatch=True,
        weak=True,
        subMatchers=[
            SM( name='root2',
                startReStr="",
                sections=['section_topology'],
                forwardMatch=True,  # The first line of the header is always skipped
                weak=True,
                # repeats=True,
                # The search is done unordered since the keywords do not appear in a specific order.
                subFlags=SM.SubFlags.Unordered,
                subMatchers=[
                    # NUMBER OF ATOMS (at_count)
                    SM(r"\s*(?P<number_of_topology_atoms>{0})\satoms".format(sformat.int), repeats=True),
                    # NUMBER OF BONDS (bd_count)
                    SM(r"\s*(?P<x_lammps_data_bd_count_store>{0})\s*bonds".format(sformat.int), repeats=True),
                    # NUMBER OF ANGLES (ag_count)
                    SM(r"\s*(?P<x_lammps_data_ag_count_store>{0})\sangles".format(sformat.int), repeats=True),
                    # NUMBER OF DIHEDRALS (dh_count)
                    SM(r"\s*(?P<x_lammps_data_dh_count_store>{0})\sdihedrals".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\simpropers".format(sformat.int), repeats=True),
                    # NUMBER OF ATOM TYPES
                    SM(r"\s*(?P<x_lammps_data_at_types_store>{0})\satom types".format(sformat.int), repeats=True),
                    # NUMBER OF BOND TYPES
                    SM(r"\s*(?P<x_lammps_data_bd_types_store>{0})\sbond types".format(sformat.int), repeats=True),
                    # # NUMBER OF ANGLE TYPES (ag_types)
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sangle types".format(sformat.int), repeats=True),
                    # # NUMBER OF DIHEDRAL TYPES (dh_types)
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sdihedral types".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\simproper types ".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sextra bond per atom".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sextra angle per atom".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sextra dihedral per atom".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sextra improper per atom ".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sextra special per atom ".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sellipsoids".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\slines".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\striangles".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sbodies".format(sformat.int), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0}\s{0})\s*xlo\sxhi".format(sformat.float), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0}\s{0})\s*ylo yhi".format(sformat.float), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0}\s{0})\s*zlo zhi".format(sformat.float), repeats=True),
                    SM(r"\s*(?P<x_lammps_dummy>{0})\sxy xz yz".format(sformat.int), repeats=True),
                      atoms,
                    velocities,
                    ellipsoids,
                    lines,
                    triangles,
                    bodies,
                    masses,
                    bonds,
                    angles,
                    dihedrals,
                    impropers,
                    pair,
                    pair_ljcoeffs,
                    bond_coeffs,
                    angle_coeffs,
                    dihedral_coeffs,
                    improper_coeffs,
                    bondbond_coeffs,
                    bondangle_coeffs,
                    middlebondtorsion_coeffs,
                    endbondtorsion_coeffs,
                    angletorsion_coeffs,
                    angleangletorsion_coeffs,
                    bondbond13_coeffs,
                    angleangle_coeffs,

                ]
            )
        ]
    )




def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
        CachingLvl: Sets the CachingLevel for the sections dos, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
        'section_run': CachingLvl,
        'section_topology': CachingLvl,
    }

    # Set all temporarly metadata to Cache as they need post-processsing.
    for name in metaInfoEnv.infoKinds:
        if name.endswith('_store'):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main(CachingLvl):
    """Main function.

    Set up everything for the parsing of the FHI-aims DOS file and run the parsing.

    Args:
        CachingLvl: Sets the CachingLevel for the sections dos, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.
    """
    # get dos file description
    LammpsDataSimpleMatcher = build_LammpsDataFileSimpleMatcher()

    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)

    # set parser info
    parserInfo = {'name': 'lammps-data-parser', 'version': '1.0'}

    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)

    # start parsing
    mainFunction(
        mainFileDescription=LammpsDataSimpleMatcher,
        metaInfoEnv=metaInfoEnv,
        parserInfo=parserInfo,
        cachingLevelForMetaName=cachingLevelForMetaName,
        superContext=LammpsDataParserContext()
    )


if __name__ == "__main__":
    main(CachingLevel.Forward)
