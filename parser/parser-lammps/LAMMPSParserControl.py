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




########################################################################################################################
########################################################################################################################
@contextmanager # SECTIONS ARE CLOSED AUTOMATICALLY
def open_section(p, name):
	gid = p.openSection(name)
	yield
	p.closeSection(name, gid)

########################################################################################################################
parser_info = {"name": "parser-lammps", "version": "1.0"}
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
               "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))

metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)



########################################################################################################################
########################################################################################################################
####### PARSED VALUES ARE ORDERED IN SECTIONS
########################################################################################################################
########################################################################################################################


# collecting atomic masses
charge_dict, charge_list, mass_dict, mass_list, mass_xyz, new_mass_list  = readChargeAndMass()
mass_dict = sorted(mass_dict.items(), key=operator.itemgetter(0))

updateAtomTypes = []
if new_mass_list:
    for i in range(len(mass_list)):
        updateAtomTypes.append([ new_mass_list[i][0], mass_list[i][0] ])
else:
    pass

print updateAtomTypes
# print topo_list_new

# ordering atomic partial charges
charge_dict = sorted(charge_dict.items(), key=operator.itemgetter(0))

at_types = len(mass_list)

def parse(fName):

    p  = JsonParseEventsWriterBackend(metaInfoEnv)
    o  = open_section
    p.startedParsingSession(fName, parser_info)


    # opening section_run
    with o(p, 'section_run'):
        p.addValue('program_name', 'LAMMPS')


        # opening section_topology
        with o(p, 'section_topology'):


            # collecting information defining different molecules
            moleculeTypeInfo, moleculeInfo, moleculeInfoResolved = assignMolecules()

            # collection list of force field functional styles
            list_of_styles = readStyles()
            pairFunctional = list_of_styles.get('pair_style')
            bondFunctional = list_of_styles.get('bond_style')
            angleFunctional = list_of_styles.get('angle_style')
            dihedralFunctional = list_of_styles.get('dihedral_style')

            # collecting covalent bond definitions
            bond_dict, bondTypeList, bond_interaction_atoms  = assignBonds(updateAtomTypes)
            bond_dict = sorted(bond_dict.items(), key=operator.itemgetter(0))
            list_of_bonds = readBonds()
            list_of_bonds = sorted(list_of_bonds.items(), key=operator.itemgetter(0))
            bd_types = len(bond_dict)

            # collecting bond angles definitions
            angle_dict, angleTypeList, angle_interaction_atoms  = assignAngles(updateAtomTypes)
            angle_dict = sorted(angle_dict.items(), key=operator.itemgetter(0))
            list_of_angles = readAngles()
            list_of_angles = sorted(list_of_angles.items(), key=operator.itemgetter(0))
            ag_types = len(angle_dict)

            # collecting dihedral angles definitions
            dihedral_dict, dihedralTypeList, dihedral_interaction_atoms  = assignDihedrals(updateAtomTypes)
            dihedral_dict = sorted(dihedral_dict.items(), key=operator.itemgetter(0))
            list_of_dihedrals = readDihedrals()
            list_of_dihedrals = sorted(list_of_dihedrals.items(), key=operator.itemgetter(0))
            dh_types = len(dihedral_dict)

            # collecting dispersion interactions ff terms
            list_of_ljs, ljs_dict  = readPairCoeff(updateAtomTypes)
            ljs_dict = sorted(ljs_dict.items(), key=operator.itemgetter(0))
            list_of_ljs = sorted(list_of_ljs.items(), key=operator.itemgetter(0))
            lj_types = len(ljs_dict)



            number_of_topology_atoms = numberOfTopologyAtoms()
            p.addValue('number_of_topology_atoms', number_of_topology_atoms)


            atom_to_molecule = []
            for i in range(number_of_topology_atoms):
                atom_to_molecule.append([moleculeInfoResolved[i][1], moleculeInfoResolved[i][3]])

            p.addValue('number_of_topology_molecules', len(moleculeInfo))
            #p.addArrayValues('atom_to_molecule', np.asarray(atom_to_molecule))


            # opening section_atom_types
            for i in range(at_types):
                #gid = p.openSection('section_atom_type')
                #p.closeSection('section_atom_type', gid)
                #print gid

                with o(p, 'section_atom_type'):
                    p.addValue('atom_type_name', [mass_xyz[i], i+1]) # Here atom_type_name is atomic number plus an integer index
                    p.addValue('atom_type_mass', mass_list[i][1])
                    p.addValue('atom_type_charge', charge_dict[i][1])
                    pass


            # opening section_molecule_type
            for i in range(len(moleculeTypeInfo)):

                with o(p, 'section_molecule_type'):
                    p.addValue('molecule_type_name', 'molecule'+'_'+str(moleculeTypeInfo[i][0]))
                    p.addValue('number_of_atoms_in_molecule', len(moleculeTypeInfo[i][1]))

                    p.addArrayValues('atom_in_molecule_to_atom_type_ref', np.asarray([x-1 for x in moleculeTypeInfo[i][1]]))


                    atom_in_molecule_name = []
                    for j in moleculeTypeInfo[i][1]:
                        atom_in_molecule_name.append([ mass_xyz[j-1], j ] ) # Here atom_in_molecule_name is atomic number plus an integer index

                    p.addValue('atom_in_molecule_name', atom_in_molecule_name)

                    atom_in_molecule_charge = []
                    for j in moleculeTypeInfo[i][1]:
                        atom_in_molecule_charge.append(charge_list[j-1][1])

                    p.addValue('atom_in_molecule_charge', atom_in_molecule_charge)
                    pass

            molecule_to_molecule_type_map = []
            for i in range(len(moleculeInfo)):
                molecule_to_molecule_type_map.append(moleculeInfo[i][1]-1) # mapping molecules to the relative section_molecule_type

            #p.addArrayValues('molecule_to_molecule_type_map', np.asarray(molecule_to_molecule_type_map))
            pass

            # opening section_interaction for covalent bonds (number_of_atoms_per_interaction = 2)
            if bd_types:

                store = []
                interaction_atoms = []
                for i in bondTypeList:
                    for j in bondTypeList:

                        store = [ [x[1], x[2]] for x in bond_interaction_atoms if x[0]==i ]
                    #print store
                    interaction_atoms.append(store)

                #print interaction_atoms

                for i in range(len(bondTypeList)):

                    with o(p, 'section_interaction'):
                        #p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0]))

                        if bondFunctional:
                            p.addValue('interaction_kind', bondFunctional)

                        int_index_store = bond_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                        else:
                            for line in int_index_store:
                                temp = sorted(map(lambda x:x-1, line))
                                interaction_atom_to_atom_type_ref.append(temp)

                        # interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                        p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', list_of_bonds[i][1])


            #if bd_types:
            #    with o(p, 'section_interaction'):
            #        p.addValue('number_of_interactions', bd_types)
            #        p.addValue('number_of_atoms_per_interaction', len(bond_dict[0][1]))

            #        if bondFunctional:
            #            p.addValue('interaction_kind', bondFunctional)

            #        int_index_store = []
            #        int_param_store = []

            #        for i in range(bd_types):
            #            int_index_store.append(bond_dict[i][1])
            #            int_param_store.append(list_of_bonds[i][1])

            #        p.addValue('interaction_atoms', int_index_store)
            #        p.addValue('interaction_parameters', int_param_store)
            #        pass


            #for i in range(bd_types):

            #    with o(p, 'section_interaction'):
            #        p.addValue('number_of_interactions', bd_types)
            #        p.addValue('number_of_atoms_per_interaction', len(bond_dict[i][1]))
            #        p.addValue('interaction_atoms', bond_dict[i][1])
            #        p.addValue('interaction_parameters', list_of_bonds[i][1])
            #        pass


            # opening section_interaction for bond angles (number_of_atoms_per_interaction = 3)
            if ag_types:

                store = []
                interaction_atoms = []
                for i in angleTypeList:
                    for j in angleTypeList:

                        store = [ [x[1], x[2], x[3]] for x in angle_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(angleTypeList)):

                    with o(p, 'section_interaction'):
                        #p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0]))

                        if bondFunctional:
                            p.addValue('interaction_kind', angleFunctional)

                        int_index_store = angle_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1]

                        else:
                            for line in int_index_store:
                                temp = map(lambda x:x-1, line)
                                interaction_atom_to_atom_type_ref.append(temp)

                        # int_index_store = angle_dict[i][1]
                        # interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1]

                        p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', list_of_angles[i][1])


            #if ag_types:
            #    with o(p, 'section_interaction'):
            #        p.addValue('number_of_interactions', ag_types)
            #        p.addValue('number_of_atoms_per_interaction', len(angle_dict[0][1]))

            #        if angleFunctional:
            #            p.addValue('interaction_kind', angleFunctional)

            #        int_index_store = []
            #        int_param_store = []

            #        for i in range(ag_types):
            #            int_index_store.append(angle_dict[i][1])
            #            int_param_store.append(list_of_angles[i][1])

            #        p.addValue('interaction_atoms', int_index_store)
            #        p.addValue('interaction_parameters', int_param_store)
            #        pass


            #for i in range(ag_types):

            #    with o(p, 'section_interaction'):
            #        p.addValue('number_of_interactions', ag_types)
            #        p.addValue('number_of_atoms_per_interaction', len(angle_dict[i][1]))
            #        p.addValue('interaction_atoms', angle_dict[i][1])
            #        p.addValue('interaction_parameters', list_of_angles[i][1])
            #        pass


            # opening section_interaction for dihedral angles (number_of_atoms_per_interaction = 4)
            if dh_types:

                store = []
                interaction_atoms = []
                for i in dihedralTypeList:
                    for j in dihedralTypeList:

                        store = [ [x[1], x[2], x[3], x[4]] for x in dihedral_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(dihedralTypeList)):

                    with o(p, 'section_interaction'):
                        #p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0]))

                        if bondFunctional:
                            p.addValue('interaction_kind', dihedralFunctional)

                        int_index_store = dihedral_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1, int_index_store[3]-1]

                        else:
                            for line in int_index_store:
                                temp = map(lambda x:x-1, line)
                                interaction_atom_to_atom_type_ref.append(temp)

                        # int_index_store = dihedral_dict[i][1]
                        # interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1, int_index_store[3]-1]

                        p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', list_of_dihedrals[i][1])


            #if dh_types:
            #    with o(p, 'section_interaction'):

            #        p.addValue('number_of_interactions', dh_types)
            #        p.addValue('number_of_atoms_per_interaction', len(dihedral_dict[0][1]))

            #        if dihedralFunctional:
            #            p.addValue('interaction_kind', dihedralFunctional)

            #        int_index_store = []
            #        int_param_store = []

            #        for i in range(dh_types):
            #            int_index_store.append(dihedral_dict[i][1])
            #            int_param_store.append(list_of_dihedrals[i][1])

            #        p.addValue('interaction_atoms', int_index_store)
            #        p.addValue('interaction_parameters', int_param_store)
            #        pass


            # opening section_interaction for dispersive interactions
            if lj_types:
                with o(p, 'section_interaction'):

                    p.addValue('number_of_interactions', lj_types)
                    p.addValue('number_of_atoms_per_interaction', len(ljs_dict[0][1]))

                    if pairFunctional:
                        p.addValue('interaction_kind', pairFunctional)

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

                    # interaction_atom_to_atom_type_ref = []
                    # for i in range(lj_types):
                    #     temp = [int_index_store[i][0]-1, int_index_store[i][1]-1]
                    #     interaction_atom_to_atom_type_ref.append(temp)

                    #p.addValue('interaction_atoms', int_index_store)
                    p.addArrayValues('interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                    p.addArrayValues('interaction_parameters', np.asarray(int_param_store))
                    pass


        # opening section_sampling_method
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
            pass



        # opening section_frame_sequence
        skip = logFileOpen()
        var, thermo_style = readLoggedThermoOutput()
        frame_length, simulation_length = simulationTime()

########################################################################################################################
# THERMO OUTPUTS FOR thermo_style = multi TO THE BACKEND

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
                #p.addValue('number_of_frames_in_sequence', frames_count-1)
                p.addValue('frame_sequence_potential_energy_stats', [pe.mean(), pe.std()])
                #p.addArrayValues('frame_sequence_potential_energy', pe)
                p.addValue('frame_sequence_kinetic_energy_stats', [ke.mean(), ke.std()])
                #p.addArrayValues('frame_sequence_kinetic_energy', ke)
                p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])
                #p.addArrayValues('frame_sequence_temperature', temp)
                p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])
                #p.addArrayValues('frame_sequence_pressure', press)

        else:
            pass

########################################################################################################################
# THERMO OUTPUTS FOR thermo_style = custom TO THE BACKEND

        if thermo_style == 'custom' and skip == False:   # Open section_frame_sequence only if an output log file is found
            from LAMMPSParserLog import pickNOMADVarsCustom


            with o(p, 'section_frame_sequence'):

                ke, pe, press, temp = pickNOMADVarsCustom()

                p.addValue('number_of_frames_in_sequence', int(simulation_length / frame_length))
                p.addValue('frame_sequence_time', [frame_length, simulation_length])

                if pe:
                    pe = np.asarray(pe)
                    p.addValue('frame_sequence_potential_energy_stats', [pe.mean(), pe.std()])
                    #p.addArrayValues('frame_sequence_potential_energy', pe)

                if ke:
                    ke = np.asarray(ke)
                    p.addValue('frame_sequence_kinetic_energy_stats', [ke.mean(), ke.std()])
                    #p.addArrayValues('frame_sequence_kinetic_energy', ke)

                if temp:
                    temp = np.asarray(temp)
                    p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])
                    #p.addArrayValues('frame_sequence_temperature', temp)

                if press:
                    press = np.asarray(press)
                    p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])
                    #p.addArrayValues('frame_sequence_pressure', press)

        else:
            pass

########################################################################################################################
# THERMO OUTPUTS FOR thermo_style = one TO THE BACKEND

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
                #p.addArrayValues('frame_sequence_temperature', temp)
                #p.addArrayValues('frame_sequence_pressure', press)

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
