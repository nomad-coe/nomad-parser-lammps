import setup_paths
import numpy as np
import math, copy
import operator
from contextlib import contextmanager

from LAMMPSParserInput import readEnsemble, readBonds, readAngles, readDihedrals, readPairCoeff, readStyles, \
                              readTPSettings, readIntegratorSettings, readLoggedThermoOutput, \
                              simulationTime

from LAMMPSParserData  import readChargeAndMass, assignBonds, assignAngles, assignDihedrals, assignMolecules, \
                              numberOfTopologyAtoms

from LAMMPSParserLog import logFileOpen

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.parser_backend import JsonParseEventsWriterBackend
import re, os, sys, json, logging




################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
@contextmanager # SECTIONS ARE CLOSED AUTOMATICALLY
def open_section(p, name):
	gid = p.openSection(name)
	yield
	p.closeSection(name, gid)

################################################################################################################################################################################################################################################

parser_info = {"name": "parser-lammps", "version": "1.0"}
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
               "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))

metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)


################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
####### PARSED VALUES ARE ORDERED IN SECTIONS ##################################################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

def parse(fName):

    p  = JsonParseEventsWriterBackend(metaInfoEnv)
    o  = open_section
    p.startedParsingSession(fName, parser_info)


    # opening section_run
    with o(p, 'section_run'):
        p.addValue('program_name', 'LAMMPS')


        # opening section_topology
        with o(p, 'section_topology'):

            ####### RETRIEVING TOPOLOGY AND FORCE FIELD INFO FROM ANCILLARY PARSERS ############################################################################################################################################################
            ####################################################################################################################################################################################################################################

            # collecting atomic masses
            charge_dict, charge_list, mass_dict, mass_list, mass_xyz, new_mass_list  = readChargeAndMass()
            mass_dict = sorted(mass_dict.items(), key=operator.itemgetter(0))
            charge_dict = sorted(charge_dict.items(), key=operator.itemgetter(0)) # ordering atomic partial charges

            updateAtomTypes = []
            if new_mass_list:
                for i in range(len(mass_list)):
                    updateAtomTypes.append([ new_mass_list[i][0], mass_list[i][0] ])
            else:
                pass
            ###

            # collecting information defining different molecules
            moleculeTypeInfo, moleculeInfo, moleculeInfoResolved = assignMolecules()
            ###

            # collecting list of force field functional styles
            list_of_styles = readStyles()
            pairFunctional = list_of_styles.get('pair_style')
            bondFunctional = list_of_styles.get('bond_style')
            angleFunctional = list_of_styles.get('angle_style')
            dihedralFunctional = list_of_styles.get('dihedral_style')
            ###

            # collecting covalent bond definitions
            bond_dict, bondTypeList, bond_interaction_atoms  = assignBonds(updateAtomTypes)
            bond_dict = sorted(bond_dict.items(), key=operator.itemgetter(0))
            list_of_bonds = readBonds()
            list_of_bonds = sorted(list_of_bonds.items(), key=operator.itemgetter(0))
            bd_types = len(bond_dict)
            ###

            # collecting bond angles definitions
            angle_dict, angleTypeList, angle_interaction_atoms  = assignAngles(updateAtomTypes)
            angle_dict = sorted(angle_dict.items(), key=operator.itemgetter(0))
            list_of_angles = readAngles()
            list_of_angles = sorted(list_of_angles.items(), key=operator.itemgetter(0))
            ag_types = len(angle_dict)
            ###

            # collecting dihedral angles definitions
            dihedral_dict, dihedralTypeList, dihedral_interaction_atoms  = assignDihedrals(updateAtomTypes)
            dihedral_dict = sorted(dihedral_dict.items(), key=operator.itemgetter(0))
            list_of_dihedrals = readDihedrals()
            list_of_dihedrals = sorted(list_of_dihedrals.items(), key=operator.itemgetter(0))
            dh_types = len(dihedral_dict)
            ###

            # collecting dispersion interactions ff terms
            list_of_ljs, ljs_dict  = readPairCoeff(updateAtomTypes)
            ljs_dict = sorted(ljs_dict.items(), key=operator.itemgetter(0))
            list_of_ljs = sorted(list_of_ljs.items(), key=operator.itemgetter(0))
            lj_types = len(ljs_dict)
            ###
            ####################################################################################################################################################################################################################################


            #### BASIC TOPOLOGY INFORMATION IN section_topology ################################################################################################################################################################################
            ####################################################################################################################################################################################################################################
            number_of_topology_atoms = numberOfTopologyAtoms()
            p.addValue('number_of_topology_atoms', number_of_topology_atoms)  # backend add number of topology atoms


            atom_to_molecule = []
            for i in range(number_of_topology_atoms):
                atom_to_molecule.append([moleculeInfoResolved[i][1], moleculeInfoResolved[i][3]])

            p.addValue('number_of_topology_molecules', len(moleculeInfo))  # backend add number of topology molecules
            p.addArrayValues('atom_to_molecule', np.asarray(atom_to_molecule))

            ####################################################################################################################################################################################################################################


            #### ATOM TYPE INFORMATION IN section_atom_type ####################################################################################################################################################################################
            ####################################################################################################################################################################################################################################
            at_types = len(mass_list)
            for i in range(at_types):

                with o(p, 'section_atom_type'):
                    p.addValue('atom_type_name', [mass_xyz[i], i+1])  # Here atom_type_name is atomic number plus an integer index identifying the atom type
                    p.addValue('atom_type_mass', mass_list[i][1])     # Atomic mass
                    p.addValue('atom_type_charge', charge_dict[i][1]) # Atomic charge, either partial or ionic
                    pass

            ####################################################################################################################################################################################################################################


            #### COVALENT BONDS INFORMATION IN section_interaction (number_of_atoms_per_interaction = 2) #######################################################################################################################################
            ####################################################################################################################################################################################################################################

            if bd_types:  #  COVALENT BONDS

                store = []
                interaction_atoms = []
                for i in bondTypeList:
                    for j in bondTypeList:

                        store = [ [x[1], x[2]] for x in bond_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(bondTypeList)):

                    with o(p, 'section_interaction'):
                        p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))     # atom indexes of bound pairs for a specific covalent bond
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))             # number of covalent bonds of this type
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0])) # number of atoms involved (2 for covalent bonds)

                        if bondFunctional:
                            p.addValue('interaction_kind', bondFunctional)  # functional form of the interaction

                        int_index_store = bond_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                        else:
                            for line in int_index_store:
                                temp = sorted(map(lambda x:x-1, line))
                                interaction_atom_to_atom_type_ref.append(temp)

                        p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', list_of_bonds[i][1])  # interaction parameters for the functional

            ####################################################################################################################################################################################################################################


            #### BOND ANGLES INFORMATION IN section_interaction (number_of_atoms_per_interaction = 3) ##########################################################################################################################################
            ####################################################################################################################################################################################################################################

            if ag_types:  # BOND ANGLES

                store = []
                interaction_atoms = []
                for i in angleTypeList:
                    for j in angleTypeList:

                        store = [ [x[1], x[2], x[3]] for x in angle_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(angleTypeList)):

                    with o(p, 'section_interaction'):
                        p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))     # atom indexes of triplets for a specific bond angle
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))             # number of bond angles of this type
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0])) # number of atoms involved (3 for bond angles)

                        if bondFunctional:
                            p.addValue('interaction_kind', angleFunctional)  # functional form of the interaction

                        int_index_store = angle_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1]

                        else:
                            for line in int_index_store:
                                temp = map(lambda x:x-1, line)
                                interaction_atom_to_atom_type_ref.append(temp)

                        p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', list_of_angles[i][1])  # interaction parameters for the functional

            ####################################################################################################################################################################################################################################


            #### DIHEDRAL ANGLES INFORMATION IN section_interaction (number_of_atoms_per_interaction = 4) ######################################################################################################################################
            ####################################################################################################################################################################################################################################

            if dh_types:  # DIHEDRAL ANGLES

                store = []
                interaction_atoms = []
                for i in dihedralTypeList:
                    for j in dihedralTypeList:

                        store = [ [x[1], x[2], x[3], x[4]] for x in dihedral_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(dihedralTypeList)):

                    with o(p, 'section_interaction'):
                        p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))     # atom indexes of quartets for a specific dihedral angle
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))             # number of dihedral angles of this type
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0])) # number of atoms involved (4 for dihedral angles)

                        if bondFunctional:
                            p.addValue('interaction_kind', dihedralFunctional)  # functional form of the interaction

                        int_index_store = dihedral_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1, int_index_store[3]-1]

                        else:
                            for line in int_index_store:
                                temp = map(lambda x:x-1, line)
                                interaction_atom_to_atom_type_ref.append(temp)

                        p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', list_of_dihedrals[i][1])  # interaction parameters for the functional
            ####################################################################################################################################################################################################################################


            #### DISPERSIVE FORCES INFORMATION IN section_interaction (number_of_atoms_per_interaction = 2) ####################################################################################################################################
            ####################################################################################################################################################################################################################################

            if lj_types:  # LJ-like interactions

                with o(p, 'section_interaction'):

                    p.addValue('number_of_defined_pair_interactions', lj_types)  # number of LJ interaction types
                    p.addValue('number_of_atoms_per_interaction', len(ljs_dict[0][1]))  # = 2 for pair interactions

                    if pairFunctional:
                        p.addValue('interaction_kind', pairFunctional)  # functional form of the interaction

                    int_index_store = []
                    int_param_store = []

                    for i in range(lj_types):
                        int_index_store.append(ljs_dict[i][1])
                        int_param_store.append(list_of_ljs[i][1])

                    interaction_atom_to_atom_type_ref = []
                    if all(isinstance(elem, list) for elem in int_index_store) == False:
                        interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                    else:
                        for line in int_index_store:
                            temp = map(lambda x:x-1, line)
                            interaction_atom_to_atom_type_ref.append(temp)

                    # p.addValue('interaction_atoms', int_index_store)
                    p.addArrayValues('pair_interaction_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                    p.addArrayValues('pair_interaction_parameters', np.asarray(int_param_store))  # interaction parameters for the functional
                    pass

            ####################################################################################################################################################################################################################################


            #### MOLECULE TYPE INFORMATION IN section_molecule_type ############################################################################################################################################################################
            ####################################################################################################################################################################################################################################

            for i in range(len(moleculeTypeInfo)):

                with o(p, 'section_molecule_type'):
                    # gindex = 0
                    p.addValue('molecule_type_name', 'molecule'+'_'+str(moleculeTypeInfo[i][0]))
                    p.addValue('number_of_atoms_in_molecule', len(moleculeTypeInfo[i][1]))

                    p.addArrayValues('atom_in_molecule_to_atom_type_ref', np.asarray([x-1 for x in moleculeTypeInfo[i][1]]))


                    atom_in_molecule_name = []
                    for j in moleculeTypeInfo[i][1]:
                        atom_in_molecule_name.append([ mass_xyz[j-1], j ] ) # Here atom_in_molecule_name is atomic number plus an integer index

                    p.addArrayValues('atom_in_molecule_name', np.asarray(atom_in_molecule_name))

                    atom_in_molecule_charge = []
                    for j in moleculeTypeInfo[i][1]:
                        atom_in_molecule_charge.append(charge_list[j-1][1])

                    p.addValue('atom_in_molecule_charge', atom_in_molecule_charge)

                    ############################################################################################################################################################################################################################


                    #### COVALENT BONDS INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 2) ######################################################################################################################

                    if bd_types:

                        toMoleculeAtomIndex  = min( moleculeTypeInfo[i][2] )

                        store = []
                        molecule_interaction_atoms = []
                        molecule_interaction_type  = []
                        for h in bondTypeList:
                            for k in bondTypeList:

                                store   = [ [x[1] - toMoleculeAtomIndex, x[2] - toMoleculeAtomIndex] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]
                                # store   = [ [x[1], x[2]] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]

                            molecule_interaction_atoms.append(store)

                        for l in bondTypeList:
                            store1  = [ x[0] for x in bond_interaction_atoms if x[0]==l and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]
                            molecule_interaction_type.append(store1)

                        molecule_interaction_type = [ x for sublist in molecule_interaction_type for x in sublist ] # from list of lists to list

                        for bond in bondTypeList:
                                if bond in molecule_interaction_type:

                                    with o(p, 'section_molecule_interaction'):

                                        p.addArrayValues('molecule_interaction_atoms', np.asarray(molecule_interaction_atoms[bond-1]))
                                        p.addValue('number_of_molecule_interactions', len(molecule_interaction_atoms[bond-1]))
                                        p.addValue('number_of_atoms_per_interaction', len(molecule_interaction_atoms[0][0]))

                                        if bondFunctional:
                                            p.addValue('molecule_interaction_kind', bondFunctional)


                                        int_index_store = bond_dict[bond-1][1]

                                        molecule_interaction_atom_to_atom_type_ref = []
                                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                                            molecule_interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                                        else:
                                            for line in int_index_store:
                                                temp = sorted(map(lambda x:x-1, line))
                                                molecule_interaction_atom_to_atom_type_ref.append(temp)

                                        p.addArrayValues('molecule_interaction_atom_to_atom_type_ref', np.asarray(molecule_interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                                        p.addValue('molecule_interaction_parameters', list_of_bonds[bond-1][1])

                    ############################################################################################################################################################################################################################


                    #### BOND ANGLE INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 3) #########################################################################################################################

                    if ag_types:

                        toMoleculeAtomIndex  = min( moleculeTypeInfo[i][2] )

                        store = []
                        molecule_interaction_atoms = []
                        molecule_interaction_type  = []
                        for h in angleTypeList:
                            for k in angleTypeList:

                                store   = [ [x[1] - toMoleculeAtomIndex, x[2] - toMoleculeAtomIndex, x[3] - toMoleculeAtomIndex] for x in angle_interaction_atoms
                                            if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] and x[3] in moleculeTypeInfo[i][2] ]
                                # store   = [ [x[1], x[2], x[3]] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]

                            molecule_interaction_atoms.append(store)

                        for l in angleTypeList:
                            store1  = [ x[0] for x in angle_interaction_atoms if x[0]==l and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]
                            molecule_interaction_type.append(store1)

                        molecule_interaction_type = [ x for sublist in molecule_interaction_type for x in sublist ] # from list of lists to list

                        for angle in angleTypeList:
                                if angle in molecule_interaction_type:

                                    with o(p, 'section_molecule_interaction'):

                                        p.addArrayValues('molecule_interaction_atoms', np.asarray(molecule_interaction_atoms[angle-1]))
                                        p.addValue('number_of_molecule_interactions', len(molecule_interaction_atoms[angle-1]))
                                        p.addValue('number_of_atoms_per_interaction', len(molecule_interaction_atoms[0][0]))

                                        if bondFunctional:
                                            p.addValue('molecule_interaction_kind', angleFunctional)


                                        int_index_store = angle_dict[angle-1][1]
                                        molecule_interaction_atom_to_atom_type_ref = []

                                        # print int_index_store, '######'

                                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                                            molecule_interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1]

                                        else:
                                            for line in int_index_store:
                                                temp = map(lambda x:x-1, line)
                                                molecule_interaction_atom_to_atom_type_ref.append(temp)

                                        p.addArrayValues('molecule_interaction_atom_to_atom_type_ref', np.asarray(molecule_interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                                        p.addValue('molecule_interaction_parameters', list_of_angles[angle-1][1])

                    ############################################################################################################################################################################################################################


                    #### DIHEDRAL ANGLE INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 4) ######################################################################################################################

                    if dh_types:

                        toMoleculeAtomIndex  = min( moleculeTypeInfo[i][2] )

                        store = []
                        molecule_interaction_atoms = []
                        molecule_interaction_type  = []
                        for h in dihedralTypeList:
                            for k in dihedralTypeList:

                                store   = [ [x[1] - toMoleculeAtomIndex, x[2] - toMoleculeAtomIndex, x[3] - toMoleculeAtomIndex, x[4] - toMoleculeAtomIndex]
                                            for x in dihedral_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2]
                                            and x[3] in moleculeTypeInfo[i][2] and x[4] in moleculeTypeInfo[i][2] ]
                                # store   = [ [x[1], x[2], x[3]] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]

                            molecule_interaction_atoms.append(store)

                        for l in dihedralTypeList:
                            store1  = [ x[0] for x in dihedral_interaction_atoms if x[0]==l
                                        and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] and x[3] in moleculeTypeInfo[i][2] and x[4] in moleculeTypeInfo[i][2] ]
                            molecule_interaction_type.append(store1)

                        molecule_interaction_type = [ x for sublist in molecule_interaction_type for x in sublist ] # from list of lists to list

                        # print molecule_interaction_type,  '#######'
                        # print molecule_interaction_atoms, '#######'

                        for dihedral in dihedralTypeList:
                                if dihedral in molecule_interaction_type:

                                    with o(p, 'section_molecule_interaction'):

                                        p.addArrayValues('molecule_interaction_atoms', np.asarray(molecule_interaction_atoms[dihedral-1]))
                                        p.addValue('number_of_molecule_interactions', len(molecule_interaction_atoms[dihedral-1]))
                                        p.addValue('number_of_atoms_per_interaction', len(molecule_interaction_atoms[0][0]))

                                        if bondFunctional:
                                            p.addValue('molecule_interaction_kind', dihedralFunctional)


                                        int_index_store = dihedral_dict[dihedral-1][1]
                                        molecule_interaction_atom_to_atom_type_ref = []

                                        # print int_index_store, '######'

                                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                                            molecule_interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1, int_index_store[3]-1]

                                        else:
                                            for line in int_index_store:
                                                temp = map(lambda x:x-1, line)
                                                molecule_interaction_atom_to_atom_type_ref.append(temp)

                                        p.addArrayValues('molecule_interaction_atom_to_atom_type_ref', np.asarray(molecule_interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                                        p.addValue('molecule_interaction_parameters', list_of_dihedrals[dihedral-1][1])

                    ############################################################################################################################################################################################################################


                    #### DISPERSIVE FORCES INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 2) ###################################################################################################################

                    if lj_types:  # LJ-like interactions

                        with o(p, 'section_molecule_interaction'):

                            if pairFunctional:
                                p.addValue('molecule_interaction_kind', pairFunctional)  # functional form of the interaction

                            int_index_store = []
                            int_param_store = []

                            for z in range(lj_types):

                                if ljs_dict[z][1][0] and ljs_dict[z][1][1] in moleculeTypeInfo[i][1]:
                                    int_index_store.append(ljs_dict[z][1])
                                    int_param_store.append(list_of_ljs[z][1])

                            interaction_atom_to_atom_type_ref = []
                            if all(isinstance(elem, list) for elem in int_index_store) == False:
                                interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                            else:
                                for line in int_index_store:
                                    temp = map(lambda x:x-1, line)
                                    interaction_atom_to_atom_type_ref.append(temp)

                            p.addValue('number_of_defined_molecule_pair_interactions', lj_types)  # number of LJ interaction types
                            p.addValue('number_of_atoms_per_molecule_interaction', len(ljs_dict[0][1]))  # = 2 for pair interactions
                            p.addArrayValues('pair_molecule_interaction_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                            p.addArrayValues('pair_molecule_interaction_parameters', np.asarray(int_param_store))  # interaction parameters for the functional

                        pass
                    ############################################################################################################################################################################################################################


            #### MAPPING THE TOPOLOGY MOLECULES TO THE RELATIVE section_molecule type ##########################################################################################################################################################

            molecule_to_molecule_type_map = []
            for i in range(len(moleculeInfo)):
                molecule_to_molecule_type_map.append(moleculeInfo[i][1]-1) # mapping molecules to the relative section_molecule_type

            p.addArrayValues('molecule_to_molecule_type_map', np.asarray(molecule_to_molecule_type_map))

            ####################################################################################################################################################################################################################################


        #### BASIC SAMPLING INFORMATION IN section_topology ################################################################################################################################################################################
        ####################################################################################################################################################################################################################################

        with o(p, 'section_sampling_method'):
            ensemble, sampling = readEnsemble()
            target_t, thermo_tau, langevin_gamma, target_p, baro_tau = readTPSettings()
            int_type, tstep, steps = readIntegratorSettings()

            p.addValue('integrator_type', int_type)
            p.addValue('integrator_dt', tstep)
            p.addValue('number_of_steps_requested', steps)

            p.addValue('sampling_method', sampling)
            p.addValue('ensemble_type', ensemble)

            if target_t:
                p.addValue('thermostat_target_temperature', target_t)
                p.addValue('thermostat_tau', thermo_tau)

            if target_p:
                p.addValue('barostat_target_pressure', target_p)
                p.addValue('barostat_tau', baro_tau)

            if langevin_gamma:
                p.addValue('thermostat_target_temperature', target_t)
                p.addValue('langevin_gamma', langevin_gamma)

        ####################################################################################################################################################################################################################################


        #### DYNAMICS FRAME INFORMATION IN section_frame_sequece ###########################################################################################################################################################################
        ####################################################################################################################################################################################################################################

        skip = logFileOpen()                                # bool var telling if an output log file is open
        var, thermo_style = readLoggedThermoOutput()        # collecting the read thermodynamic output
        frame_length, simulation_length = simulationTime()  # frame and simulation time length

        #### THERMO OUTPUTS FOR thermo_style = multi TO THE BACKEND

        if thermo_style == 'multi' and skip == False:   # Open section_frame_sequence only if an output log file is found
            from LAMMPSParserLog   import readFrames, readPotEnergy, readKinEnergy, readPressure, readVolume

            with o(p, 'section_frame_sequence'):
                frames_count = readFrames()
                pe = readPotEnergy()
                ke, temp = readKinEnergy()
                press = readPressure()
                vol = readVolume()

                p.addValue('number_of_frames_in_sequence', int(simulation_length / frame_length))
                p.addValue('frame_sequence_time', [frame_length, simulation_length])
                p.addValue('number_of_frames_in_sequence', frames_count-1)
                p.addValue('frame_sequence_potential_energy_stats', [pe.mean(), pe.std()])
                # p.addArrayValues('frame_sequence_potential_energy', pe)
                p.addValue('frame_sequence_kinetic_energy_stats', [ke.mean(), ke.std()])
                # p.addArrayValues('frame_sequence_kinetic_energy', ke)
                p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])
                # p.addArrayValues('frame_sequence_temperature', temp)
                p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])
                # p.addArrayValues('frame_sequence_pressure', press)

        else:
            pass


        #### THERMO OUTPUTS FOR thermo_style = custom TO THE BACKEND

        if thermo_style == 'custom' and skip == False:   # Open section_frame_sequence only if an output log file is found
            from LAMMPSParserLog import pickNOMADVarsCustom


            with o(p, 'section_frame_sequence'):

                ke, pe, press, temp = pickNOMADVarsCustom()

                p.addValue('number_of_frames_in_sequence', int(simulation_length / frame_length))
                p.addValue('frame_sequence_time', [frame_length, simulation_length])

                if pe:
                    pe = np.asarray(pe)
                    p.addValue('frame_sequence_potential_energy_stats', [pe.mean(), pe.std()])
                    # p.addArrayValues('frame_sequence_potential_energy', pe)

                if ke:
                    ke = np.asarray(ke)
                    p.addValue('frame_sequence_kinetic_energy_stats', [ke.mean(), ke.std()])
                    # p.addArrayValues('frame_sequence_kinetic_energy', ke)

                if temp:
                    temp = np.asarray(temp)
                    p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])
                    # p.addArrayValues('frame_sequence_temperature', temp)

                if press:
                    press = np.asarray(press)
                    p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])
                    # p.addArrayValues('frame_sequence_pressure', press)

        else:
            pass


        #### THERMO OUTPUTS FOR thermo_style = one TO THE BACKEND

        if thermo_style == 'one' and skip == False:   # Open section_frame_sequence only if an output log file is found
            from LAMMPSParserLog import pickNOMADVarsOne


            with o(p, 'section_frame_sequence'):

                press, temp = pickNOMADVarsOne()
                temp = np.asarray(temp)
                press = np.asarray(press)

                p.addValue('number_of_frames_in_sequence', int(simulation_length / frame_length))
                p.addValue('frame_sequence_time', [frame_length, simulation_length])

                p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])
                p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])
                # p.addArrayValues('frame_sequence_temperature', temp)
                # p.addArrayValues('frame_sequence_pressure', press)

        else:
            pass


########################################################################################################################
    p.finishedParsingSession("ParseSuccess", None)    # PARSING FINISHED


########################################################################################################################
# POIT TO A FILE TO PARSE FROM COMMAND LINE
if __name__ == '__main__':
    import sys
    fName = sys.argv[1]
    parse(fName)
