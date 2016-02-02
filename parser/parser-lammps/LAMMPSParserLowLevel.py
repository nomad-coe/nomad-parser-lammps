import setup_paths
import numpy as np
import math
from contextlib import contextmanager
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.parser_backend import JsonParseEventsWriterBackend
import re, os, sys, json, logging


@contextmanager
def open_section(p, name):
	gid = p.openSection(name)
	yield
	p.closeSection(name, gid)


parser_info = {"name": "parser-lammps", "version": "1.0"}
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))

metaInfoEnv, warns = loadJsonFile(filePath=metaInfoPath,
                                  dependencyLoader=None,
                                  extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                  uri=None)



def parse(filename):
	p = JsonParseEventsWriterBackend(metaInfoEnv)
	o = open_section
	#r = Reader(filename)
	p.startedParsingSession(filename, parser_info)

	lines = open(filename).readlines()
	store_input = []
	for line in lines:
		line = line.strip('\n' + '').split(' ')
		line = filter(None, line)
		# If line is just empty
		if line != []:
			pass
			store_input.append(line)

	with o(p, 'section_run'):
		p.addValue('program_name', 'LAMMPS')


		for line in store_input:
			if 'variable' and 'dt' in line:
				dt_list = str(line).split()
				dt = dt_list.pop()

		p.addValue('program_version', dt)

	#print lines


	p.finishedParsingSession("ParseSuccess", None)




if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    parse(filename)
