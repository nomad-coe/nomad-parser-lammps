import setup_paths
import numpy as np
import math
import operator
from contextlib import contextmanager

from LAMMPSParserInput import readEnsemble, readBonds, readAngles, readDihedrals, \
                              readTPSettings, readIntegratorSettings, readLoggedThermoOutput, \
                              simulationTime

from LAMMPSParserData  import readMass, readCharge, assignBonds, assignAngles, assignDihedrals

from LAMMPSParserLog import logFileOpen

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.parser_backend import JsonParseEventsWriterBackend
import re, os, sys, json, logging


@contextmanager
def open_section(p, name):
	gid = p.openSection(name)
	yield
	p.closeSection(name, gid)




parser_info = {"name": "parser-lammps", "version": "1.0"}
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/forcefield.nomadmetainfo.json"))

metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)



def parse(filename):

    p  = JsonParseEventsWriterBackend(metaInfoEnv)
    o  = open_section
    p.startedParsingSession(filename, parser_info)


    # opening section_run
    with o(p, 'section_run'):
        p.addValue('program_name', 'LAMMPS')


        # opening section_topology
        with o(p, 'section_topology'):
            # collecting atomic masses
            mass_dict, mass_list, mass_xyz  = readMass()
            mass_dict = sorted(mass_dict.items(), key=operator.itemgetter(0))
            at_types = len(mass_dict)

            # collecting atomic partial charges
            charge_dict, charge_list        = readCharge()
            charge_dict = sorted(charge_dict.items(), key=operator.itemgetter(0))

            # collecting covalent bond definitions
            bond_dict  = assignBonds()
            bond_dict = sorted(bond_dict.items(), key=operator.itemgetter(0))
            list_of_bonds = readBonds()
            list_of_bonds = sorted(list_of_bonds.items(), key=operator.itemgetter(0))
            bd_types = len(bond_dict)

            # collecting bond angles definitions
            angle_dict  = assignAngles()
            angle_dict = sorted(angle_dict.items(), key=operator.itemgetter(0))
            list_of_angles = readAngles()
            list_of_angles = sorted(list_of_angles.items(), key=operator.itemgetter(0))
            ag_types = len(angle_dict)

            # collecting dihedral angles definitions
            dihedral_dict  = assignDihedrals()
            dihedral_dict = sorted(dihedral_dict.items(), key=operator.itemgetter(0))
            list_of_dihedrals = readDihedrals()
            list_of_dihedrals = sorted(list_of_dihedrals.items(), key=operator.itemgetter(0))
            dh_types = len(dihedral_dict)

            p.addValue('number_of_topology_atoms', at_types)
            pass


            # opening section_atom_types
            for i in range(at_types):

                with o(p, 'section_atom_type'):
                    p.addValue('atom_type_name', [mass_xyz[i], i+1]) # Here atom_type_name is atomic number plus a integer index
                    p.addValue('atom_type_mass', mass_dict[i][1])
                    p.addValue('atom_type_charge', charge_dict[i][1])
                    pass


            # opening section_interaction for covalent bonds (number_of_atoms_per_interaction = 2)
            if bd_types:
                with o(p, 'section_interaction'):
                    p.addValue('number_of_interactions', bd_types)
                    p.addValue('number_of_atoms_per_interaction', len(bond_dict[0][1]))

                    int_index_store = []
                    int_param_store = []

                    for i in range(bd_types):
                        int_index_store.append(bond_dict[i][1])
                        int_param_store.append(list_of_bonds[i][1])

                    p.addValue('interaction_atoms', int_index_store)
                    p.addValue('interaction_parameters', int_param_store)
                    pass


            #for i in range(bd_types):

            #    with o(p, 'section_interaction'):
            #        p.addValue('number_of_interactions', bd_types)
            #        p.addValue('number_of_atoms_per_interaction', len(bond_dict[i][1]))
            #        p.addValue('interaction_atoms', bond_dict[i][1])
            #        p.addValue('interaction_parameters', list_of_bonds[i][1])
            #        pass


            # opening section_interaction for bond angles (number_of_atoms_per_interaction = 3)
            if ag_types:
                with o(p, 'section_interaction'):
                    p.addValue('number_of_interactions', ag_types)
                    p.addValue('number_of_atoms_per_interaction', len(angle_dict[0][1]))

                    int_index_store = []
                    int_param_store = []

                    for i in range(ag_types):
                        int_index_store.append(angle_dict[i][1])
                        int_param_store.append(list_of_angles[i][1])

                    p.addValue('interaction_atoms', int_index_store)
                    p.addValue('interaction_parameters', int_param_store)
                    pass


            #for i in range(ag_types):

            #    with o(p, 'section_interaction'):
            #        p.addValue('number_of_interactions', ag_types)
            #        p.addValue('number_of_atoms_per_interaction', len(angle_dict[i][1]))
            #        p.addValue('interaction_atoms', angle_dict[i][1])
            #        p.addValue('interaction_parameters', list_of_angles[i][1])
            #        pass


            # opening section_interaction for dihedral angles (number_of_atoms_per_interaction = 4)
            if dh_types:
                with o(p, 'section_interaction'):

                    p.addValue('number_of_interactions', dh_types)
                    p.addValue('number_of_atoms_per_interaction', len(dihedral_dict[0][1]))

                    int_index_store = []
                    int_param_store = []

                    for i in range(dh_types):
                        int_index_store.append(dihedral_dict[i][1])
                        int_param_store.append(list_of_dihedrals[i][1])

                    p.addValue('interaction_atoms', int_index_store)
                    p.addValue('interaction_parameters', int_param_store)
                    pass



        # opening section_sampling_method
        with o(p, 'section_sampling_method'):
            ensemble, sampling = readEnsemble()
            target_t, thermo_tau, target_p, baro_tau = readTPSettings()
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
                p.addValue('frame_sequence_kinetic_energy_stats', [ke.mean(), ke.std()])
                p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])
                p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])

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

                if ke:
                    ke = np.asarray(ke)
                    p.addValue('frame_sequence_kinetic_energy_stats', [ke.mean(), ke.std()])

                if temp:
                    temp = np.asarray(temp)
                    p.addValue('frame_sequence_temperature_stats', [temp.mean(), temp.std()])

                if press:
                    press = np.asarray(press)
                    p.addValue('frame_sequence_pressure_stats', [press.mean(), press.std()])

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

        else:
            pass



########################################################################################################################

    p.finishedParsingSession("ParseSuccess", None)    # PARSING FINISHED




if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    parse(filename)
