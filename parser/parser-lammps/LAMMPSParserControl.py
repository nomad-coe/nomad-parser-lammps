import setup_paths
import numpy as np
import math
import operator
from contextlib import contextmanager
from LAMMPSParserInput import readEnsemble, readPairCoeff, readBonds
from LAMMPSParserData import readMass, readCharge, assignBonds
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

    with o(p, 'section_run'):
        p.addValue('program_name', 'LAMMPS')

        with o(p, 'section_topology'):
            mass_dict, mass_list, mass_xyz, at_types  = readMass()
            charge_dict, charge_list        = readCharge()
            mass_dict = sorted(mass_dict.items(), key=operator.itemgetter(0))
            charge_dict = sorted(charge_dict.items(), key=operator.itemgetter(0))
            bond_dict  = assignBonds()
            bond_dict = sorted(bond_dict.items(), key=operator.itemgetter(0))
            bd_types = len(bond_dict)
            list_of_bonds = readBonds()
            list_of_bonds = sorted(list_of_bonds.items(), key=operator.itemgetter(0))



            p.addValue('number_of_topology_atoms', len(mass_dict))
            pass

            for i in range(at_types):

                with o(p, 'section_atom_type'):
                    p.addValue('atom_type_name', mass_xyz[i])           # Here I use the atomic number
                    p.addValue('atom_type_mass', mass_dict[i][1])
                    p.addValue('atom_type_charge', charge_dict[i][1])
                    pass


            #with o(p, 'section_atom_type'):
            #    mass_dict, mass_list, mass_xyz  = readMass()
            #    charge_dict, charge_list        = readCharge()
            #    p.addValue('number_of_topology_atoms', len(mass_list))
            #    p.addValue('atom_type_name', mass_xyz)
            #    p.addValue('atom_type_mass', sorted(mass_dict.items(), key=operator.itemgetter(0)))
            #    p.addValue('atom_type_charge', sorted(charge_dict.items(), key=operator.itemgetter(0)))
            #    pass



            for i in range(bd_types):

                with o(p, 'section_interaction'):
                    p.addValue('interaction_atoms', bond_dict[i][1])
                    p.addValue('interaction_parameters', list_of_bonds[i])
                    pass

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
