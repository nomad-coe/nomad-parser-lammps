import setup_paths
import numpy as np
import math
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import re, os, sys, json, logging


class LAMMPSParserContext(object):

	def __init__(self):
		self.cell	= []
		self.parser = None

	def startedParsing(self, path, parser):
    		self.parser = parser
		return

	def onClose_section_run(self, backend, gIndex, section):
		print "<onClose_section_run>"
		return

	def onClose_section_method(self, backend, gIndex, section):
		print "<onClose_section_method>"
		return

	def onClose_section_system_description(self, backend, gIndex, section):
		print "<onClose_section_system_description>"
		return

	def onClose_lammps_section_md_molecule_type(self, backend, gIndex, section):
		print "<onClose_molecule_type>"
		return

def test_adhoc(parser):
	#return None
	for i in range(3):
		ln = parser.fIn.readline()
		print ln
		parser.fIn.pushbackLine(ln)
	#print len(parser.backend.openSections)



mainFileDescription = SM(name = 'root',
	 weak = True,
	 startReStr = "",
	 subMatchers = [
	 
		SM(name = 'newRun',
		   startReStr = r"units\sreal\s*",
		   repeats = False,
		   required = True,
		   forwardMatch = True,
		   adHoc = test_adhoc,
		   sections   = ['section_run'],
		   subMatchers = [

                          ])  # close section_run
])  # close root


# THIS PARSER
parserInfo = {'name':'lammps-parser', 'version': '1.0'}

# LOAD METADATA
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)

# CUSTOMIZE CACHING
cachingLevelForMetaName = { }

if __name__ == "__main__":

	mainFunction(mainFileDescription, 
    	metaInfoEnv, 
    	parserInfo, 
    	superContext = LAMMPSParserContext(),
    	cachingLevelForMetaName = cachingLevelForMetaName, 
    	onClose = {},
		defaultSectionCachingLevel = False)

