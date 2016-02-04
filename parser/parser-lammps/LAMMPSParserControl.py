import setup_paths
import numpy as np
import math
from contextlib import contextmanager
from LAMMPSParserInput import readEnsemble, readPairCoeff
from LAMMPSParserData import readMass
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
            pass

            with o(p, 'section_atom_type'):
                mass_dict, mass_list, mass_xyz  = readMass()
                p.addValue('atom_type_name', mass_list)
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
