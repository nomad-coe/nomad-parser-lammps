import setup_paths
import numpy as np
import math
import operator
from contextlib import contextmanager
from LAMMPSParserInput import readEnsemble, readPairCoeff, readBonds, readAngles
from LAMMPSParserData import readMass, readCharge, assignBonds, assignAngles
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

    p = JsonParseEventsWriterBackend(metaInfoEnv)
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
            for i in range(bd_types):

                with o(p, 'section_interaction'):
                    p.addValue('number_of_interactions', bd_types)
                    p.addValue('number_of_atoms_per_interaction', len(bond_dict[i][1]))
                    p.addValue('interaction_atoms', bond_dict[i][1])
                    p.addValue('interaction_parameters', list_of_bonds[i][1])
                    pass


            # opening section_interaction for bond angles (number_of_atoms_per_interaction = 3)
            for i in range(ag_types):

                with o(p, 'section_interaction'):
                    p.addValue('number_of_interactions', ag_types)
                    p.addValue('number_of_atoms_per_interaction', len(angle_dict[i][1]))
                    p.addValue('interaction_atoms', angle_dict[i][1])
                    p.addValue('interaction_parameters', list_of_angles[i][1])
                    pass

        # opening section_sampling_method
        with o(p, 'section_sampling_method'):
            ensemble, sampling = readEnsemble()
            p.addValue('sampling_method', sampling)
            p.addValue('ensemble_type', ensemble)
            pass



    p.finishedParsingSession("ParseSuccess", None)




if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    parse(filename)
