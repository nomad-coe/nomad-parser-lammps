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

import setup_paths

import logging
import os


from LammpsCommon import get_metaInfo
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.simple_parser import mainFunction, AncillaryParser
import LammpsDataParser

nformat = {
    'float': r"[-+]?\d*\.\d+",
    'int': r"[-+]?\d+",
    'file': r"[0-9a-zA-Z_.]+"
}



# import numpy as np
# import setup_paths
# import logging, os, re, sys
# import nomadcore.ActivateLogging

############################################################
# This is the parser for the input file of LAMMPS.
############################################################

logger = logging.getLogger(name="nomad.LammpsInputParser")

# LAMMPS input script
# (1) LAMMPS does not read your entire input script and then perform a simulation with all the settings. Rather, the input script is read one line at a time and each command takes effect when it is read.
# default values ????

# Parsing rules: (source: http://lammps.sandia.gov/doc/Section_commands.html#parsing-rules)
# (1) If the last printable character on the line is a "&" character, the command is assumed to continue on the next line.
# (2) All characters from the first "#" character onward are treated as comment and discarded.
# (3) The line is searched repeatedly for $ characters, which indicate variables that are replaced with a text string.
# (4) The line is broken into "words" separated by whitespace (tabs, spaces). Note that words can thus contain letters, digits, underscores, or punctuation characters.
# (5) The first word is the command name. All successive words in the line are arguments.
# (6) If you want text with spaces to be treated as a single argument, it can be enclosed in either single or double or triple quotes.
#
# A LAMMPS input script typically has 4 parts:
#
# (1) Initialization
# (2) Atom definition
# (3) Settings
# (4) Run a simulation
#


def assignMolecules(bond_list, topo_list, at_count):  # FINDING INDIVIDUAL MOLECULES FROM BONDING PATTERN

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


    moleculeInfoResolved = []
    for i in range(0, at_count):
        for j in range(len(moleculeInfo)):

            if i+1 in moleculeInfo[j][2]:
                moleculeInfoResolved.append([ i+1, moleculeInfo[j][0], moleculeInfo[j][1], atomPositionInMoleculeList[i] ])


    # print moleculeTypeInfo

    return (moleculeTypeInfo, moleculeInfo, moleculeInfoResolved)







class LammpsInputParserContext(object):
    """Context for parsing LAMMPS input file.

    Attributes:
        # dos_energies: Stores parsed energies.
        # dos_values: Stores parsed DOS values.

    This class keeps tracks of several lammps settings to adjust the parsing to them.
    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of writing and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """

    def __init__(self):
        pass



    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName

        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv

        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()


    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        self.units='lj'



    def compile_data_parser(self):
        """Instantiate superContext and construct parser for external data file.
        """
        self.dataSuperContext = LammpsDataParser.LammpsDataParserContext(False)
        self.dataParser = AncillaryParser(
            fileDescription = LammpsDataParser.build_LammpsDataFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = LammpsDataParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = self.dataSuperContext
        )

    def onClose_section_run(self, backend, gIndex, section):
        pass

    def onClose_section_topology(self, backend, gIndex, section):
        """Trigger called when section_topology is closed.

        Input data is parsed from external file but the result is only stored and written later in lammps_topology
        in section_run.
        """

        tmp = section['x_lammps_dummy_text']


        if section['x_lammps_input_units_store'] is not None:
            self.units = section['x_lammps_input_units_store']


        # construct file name
        data_file = section['x_lammps_data_file_store']

        if data_file is not None:
            self.compile_data_parser()

            dirName = os.path.dirname(os.path.abspath(self.fName))
            fName = os.path.normpath(os.path.join(dirName, data_file[0]))

            try:
                with open(fName) as fIn:
                    # construct parser for DOS file if not present
                    if self.dataSuperContext is None or self.dataParser is None:
                        self.compile_data_parser()
                    # parse DOS file
                    self.dataParser.parseFile(fIn)
                    pass
                    # set flag that DOS was parsed successfully and store values
                    # if self.dosSuperContext.dos_energies is not None and self.dosSuperContext.dos_values is not None:
                    #     self.dosFound = True
                    #     self.dos_energies = self.dosSuperContext.dos_energies
                    #     self.dos_values = self.dosSuperContext.dos_values
                    # else:
                    #     logger.error("DOS parsing unsuccessful. Parsing of file %s in directory '%s' did not yield energies or values for DOS." % (data_file, dirName))
            except IOError:
                logger.error("DOS parsing unsuccessful. Could not find %s file in directory '%s'." % (data_file, dirName))

    def onClose_section_dos(self, backend, gIndex, section):
        """Trigger called when section_dos is closed.

        Store the parsed values and write them if writeMetaData is True.
        """
        # extract energies
        dos_energies = section['fhi_aims_dos_energy']
        if dos_energies is not None:
            pass
            # self.dos_energies = np.asarray(dos_energies)
        # extract dos values
        dos_values = []
        if section['fhi_aims_dos_value_string'] is not None:
            for string in section['fhi_aims_dos_value_string']:
                strings = string.split()
                dos_values.append(map(float, strings))
        if dos_values:
            # need to transpose array since its shape is [number_of_spin_channels,n_dos_values] in the metadata
            pass
            # self.dos_values = np.transpose(np.asarray(dos_values))
        # write metadata only if values were found for both quantities
        if self.writeMetaData:
            if dos_energies is not None and dos_values:
                pass
                # backend.addArrayValues('dos_energies', self.dos_energies)
                # backend.addArrayValues('dos_values', self.dos_values)


def build_LammpsInputFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the DOS file of FHI-aims.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses DOS file of FHI-aims.
    """

    # A LAMMPS input script typically has 4 parts:
    #
    # (1) Initialization
    # (2) Atom definition
    # (3) Settings
    # (4) Run a simulation


    # Initialization
    # Set parameters that need to be defined before atoms are created or read-in from a file.

    # initialization = [
    #     SM(r"(?P<x_lammps_dummy_text>\s*dimension)"),
    #     SM(r"(?P<x_lammps_dummy_text>\s*atom_modify)"),
    #     SM(r"(?P<x_lammps_dummy_text>\s*atom_style)"),
    #     SM(r"(?P<x_lammps_dummy_text>\s*boundary)"),
    #     SM(r"(?P<x_lammps_dummy_text>\s*newton)"),
    #     SM(r"(?P<x_lammps_dummy_text>\s*processors)"),
    #     SM(r"(?P<x_lammps_dummy_text>\s*units)")
    # ]


    # If force-field parameters appear in the files that will be read, these commands tell LAMMPS what kinds of force fields are being used:
    # pair_style, bond_style, angle_style, dihedral_style, improper_style.

    # initialization_force_field = [
    #     SM(r"\s*(?P<x_lammps_dummy_text>pair_style)"),
    #     SM(r"\s*(?P<x_lammps_dummy_text>bond_style)"),
    #     SM(r"\s*(?P<x_lammps_dummy_text>angle_style)"),
    #     SM(r"\s*(?P<x_lammps_dummy_text>dihedral_style)"),
    #     SM(r"\s*(?P<x_lammps_dummy_text>improper_style)")
    # ]


    # Atom definition
    # There are 3 ways to define atoms in LAMMPS. Read them in from a data or restart file via the read_data or read_restart commands.
    # These files can contain molecular topology information. Or create atoms on a lattice (with no molecular topology), using these commands:
    # lattice, region, create_box, create_atoms.
    # The entire set of atoms can be duplicated to make a larger simulation using the replicate command.

    # atom_definitation = SM(
    #     name='atom_definitation',
    #     startReStr="",
    #     # forwardMatch=True,
    #     weak=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>create_atoms)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>create_box)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>lattice)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>read_data)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>read_dump)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>read_restart)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>region)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>replicate)")
    #     ]
    #
    # )



    # Settings
    # Once atoms and molecular topology are defined, a variety of settings can be specified: force field coefficients, simulation parameters, output options, etc.
    # Settings:

    # settings_other = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>comm_style)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>mass)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>set)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>velocity)")
    #     ]
    # )
    #

    # settings_force_field_coefficients = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>pair_coeff)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>bond_coeff)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>angle_coeff)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>dihedral_coeff)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>improper_coeff)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>kspace_style)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>kspace_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>dielectric)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>special_bonds)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>pair_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>pair_write)")
    #     ]
    # )

    # settings_simulation_parameters = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>neighbor)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>neigh_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>group)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>timestep)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>reset_timestep)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>run_style)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>min_style)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>min_modify)"),
    #     ]
    # )

    # settings_fixes = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>fix)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>fix_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>unfix)"),
    #     ]
    # )


    # settings_output_options = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>dump)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>dump image)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>dump_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>dump movie)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>restart)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>thermo)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>thermo_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>thermo_style)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>undump)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>write_data)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>write_dump)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>write_restart)"),
    #     ]
    # )

    # Run a simulation (actions)
    # A molecular dynamics simulation is run using the run command. Energy minimization (molecular statics) is performed using the minimize command. A parallel tempering (replica-exchange) simulation can be run using the temper command.

    # actions = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>delete_atoms)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>delete_bonds)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>displace_atoms)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>change_box)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>minimize)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>neb prd)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>rerun)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>run)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>temper)"),
    #     ]
    # )

    # computes = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>compute)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>compute_modify)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>uncompute)"),
    #     ]
    # )

    # miscellaneous = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     # forwardMatch=True,
    #     subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         SM(r"\s*(?P<x_lammps_dummy_text>clear)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>echo)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>if)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>include)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>jump)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>label)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>log)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>next)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>print)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>shell)"),
    #         SM(r"\s*(?P<x_lammps_dummy_text>variable)"),
    #     ]
    # )

    #
    # header = SM(
    #     name='lammps-input-header',
    #     startReStr="",
    #     forwardMatch=True,
    #     sections=['section_topology'],
    #     subMatchers=[
    #         SM(r"\s*read_data\s(?P<x_lammps_data_file>{[file]})".format(nformat))
    #     ]
    # )

    ########################################
    # return main Parser


    initialization = SM(
        name="initialization",
        startReStr="",
        # endReStr = r"\s*Have a nice day\.",
        # repeats = True,
        # required=True,
        # forwardMatch=True,
        sections=['section_topology'],
        # forwardMatch=True,
        # weak=True,
        repeats=True,
        # subFlags=SM.SubFlags.Unordered,
        subMatchers=[
            SM(r"(?P<x_lammps_dummy_text>\s*dimension)"),
            SM(r"(?P<x_lammps_dummy_text>\s*atom_modify)"),
            SM(r"(?P<x_lammps_dummy_text>\s*atom_style)"),
            SM(r"(?P<x_lammps_dummy_text>\s*boundary)"),
            SM(r"(?P<x_lammps_dummy_text>\s*newton)"),
            SM(r"(?P<x_lammps_dummy_text>\s*processors)"),
            SM(r"(?P<x_lammps_dummy_text>\s*units)"),

            SM(r"\s*(?P<x_lammps_dummy_text>pair_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>bond_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>angle_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>dihedral_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>improper_style)"),

            SM(r"\s*(?P<x_lammps_dummy_text>create_atoms)"),
            SM(r"\s*(?P<x_lammps_dummy_text>create_box)"),
            SM(r"\s*(?P<x_lammps_dummy_text>lattice)"),
            SM(r"\s*(?P<x_lammps_dummy_text>read_data)"),
            SM(r"\s*(?P<x_lammps_dummy_text>read_dump)"),
            SM(r"\s*(?P<x_lammps_dummy_text>read_restart)"),
            SM(r"\s*(?P<x_lammps_dummy_text>region)"),
            SM(r"\s*(?P<x_lammps_dummy_text>replicate)"),

            SM(r"\s*(?P<x_lammps_dummy_text>comm_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>mass)"),
            SM(r"\s*(?P<x_lammps_dummy_text>set)"),
            SM(r"\s*(?P<x_lammps_dummy_text>velocity)"),

            SM(r"\s*(?P<x_lammps_dummy_text>pair_coeff)"),
            SM(r"\s*(?P<x_lammps_dummy_text>bond_coeff)"),
            SM(r"\s*(?P<x_lammps_dummy_text>angle_coeff)"),
            SM(r"\s*(?P<x_lammps_dummy_text>dihedral_coeff)"),
            SM(r"\s*(?P<x_lammps_dummy_text>improper_coeff)"),
            SM(r"\s*(?P<x_lammps_dummy_text>kspace_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>kspace_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>dielectric)"),
            SM(r"\s*(?P<x_lammps_dummy_text>special_bonds)"),
            SM(r"\s*(?P<x_lammps_dummy_text>pair_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>pair_write)"),

            SM(r"\s*(?P<x_lammps_dummy_text>neighbor)"),
            SM(r"\s*(?P<x_lammps_dummy_text>neigh_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>group)"),
            SM(r"\s*(?P<x_lammps_dummy_text>timestep)"),
            SM(r"\s*(?P<x_lammps_dummy_text>reset_timestep)"),
            SM(r"\s*(?P<x_lammps_dummy_text>run_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>min_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>min_modify)"),

            SM(r"\s*(?P<x_lammps_dummy_text>fix)"),
            SM(r"\s*(?P<x_lammps_dummy_text>fix_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>unfix)"),

            SM(r"\s*(?P<x_lammps_dummy_text>dump)"),
            SM(r"\s*(?P<x_lammps_dummy_text>dump image)"),
            SM(r"\s*(?P<x_lammps_dummy_text>dump_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>dump movie)"),
            SM(r"\s*(?P<x_lammps_dummy_text>restart)"),
            SM(r"\s*(?P<x_lammps_dummy_text>thermo)"),
            SM(r"\s*(?P<x_lammps_dummy_text>thermo_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>thermo_style)"),
            SM(r"\s*(?P<x_lammps_dummy_text>undump)"),
            SM(r"\s*(?P<x_lammps_dummy_text>write_data)"),
            SM(r"\s*(?P<x_lammps_dummy_text>write_dump)"),
            SM(r"\s*(?P<x_lammps_dummy_text>write_restart)"),
        ]
    )






    actions = SM(
        name='lammps-input-header',
        startReStr="",
        forwardMatch=True,
        sections=['section_topology'],
        subFlags=SM.SubFlags.Unordered,
        subMatchers=[
            SM(r"\s*(?P<x_lammps_dummy_text>delete_atoms)"),
            SM(r"\s*(?P<x_lammps_dummy_text>delete_bonds)"),
            SM(r"\s*(?P<x_lammps_dummy_text>displace_atoms)"),
            SM(r"\s*(?P<x_lammps_dummy_text>change_box)"),
            SM(r"\s*(?P<x_lammps_dummy_text>minimize)"),
            SM(r"\s*(?P<x_lammps_dummy_text>neb prd)"),
            SM(r"\s*(?P<x_lammps_dummy_text>rerun)"),
            SM(r"\s*(?P<x_lammps_dummy_text>run)"),
            SM(r"\s*(?P<x_lammps_dummy_text>temper)"),

            SM(r"\s*(?P<x_lammps_dummy_text>compute)"),
            SM(r"\s*(?P<x_lammps_dummy_text>compute_modify)"),
            SM(r"\s*(?P<x_lammps_dummy_text>uncompute)"),

            SM(r"\s*(?P<x_lammps_dummy_text>clear)"),
            SM(r"\s*(?P<x_lammps_dummy_text>echo)"),
            SM(r"\s*(?P<x_lammps_dummy_text>if)"),
            SM(r"\s*(?P<x_lammps_dummy_text>include)"),
            SM(r"\s*(?P<x_lammps_dummy_text>jump)"),
            SM(r"\s*(?P<x_lammps_dummy_text>label)"),
            SM(r"\s*(?P<x_lammps_dummy_text>log)"),
            SM(r"\s*(?P<x_lammps_dummy_text>next)"),
            SM(r"\s*(?P<x_lammps_dummy_text>print)"),
            SM(r"\s*(?P<x_lammps_dummy_text>shell)"),
            SM(r"\s*(?P<x_lammps_dummy_text>variable)"),
        ]
    )

    # test case: try to separate each basic steps in the input file
    # not yet working:
    # - if I try to close a section using endStr its eat my line
    # - if I try to catch the line with the next SM
    return SM(
        name="root",
        startReStr="",
        sections=['section_run'],
        fixedStartValues={'program_name': 'LAMMPS'},
        forwardMatch=True,
        weak=True,
        # subFlags=SM.SubFlags.Unordered,
        subMatchers=[
            SM(
                startReStr="",
                # endReStr=r"pair_style.*",
                # sections=['test'],
                forwardMatch=True,
                subFlags=SM.SubFlags.Unordered,
                # weak=True,
                repeats=True,
                subMatchers=[
                    SM(
                        r"(?P<x_lammps_dummy_text>a.*)",
                        repeats=True,
                    ),
                    SM(
                        r"(?P<x_lammps_dummy_text>b.*)",
                        repeats=True,
                    )
                ]
            ),
            # SM(r"(?P<x_lammps_dummy_text>pair_style.*)"),
            SM(
                startReStr="",
                # endReStr=r"pair_style.*",
                forwardMatch=True,
                subFlags=SM.SubFlags.Unordered,
                # weak=True,
                repeats=True,
                subMatchers=[

                    SM(r"(?P<x_lammps_dummy_text>pair_style.*)")
                ]
            )
        ]
    )

    # return SM(
    #     name='Root',
    #     startReStr="",
    #     forwardMatch=True,
    #     weak=True,
    #     sections=['section_run'],
    #     fixedStartValues={'program_name': 'LAMMPS'},
    #     # subFlags=SM.SubFlags.Unordered,
    #     subMatchers=[
    #         initialization,
    #         actions
    #         # SM(
    #         #     name="atom description",
    #         #     startReStr="",
    #         #     # forwardMatch=True,
    #         #     # weak=True,
    #         #     sections=['section_topology'],
    #         #     subMatchers=[
    #         #         atom_definitation
    #         #     ]
    #         # ),
    #         # SM(
    #         #     name="repetable section",
    #         #     startReStr="",
    #         #     forwardMatch=True,
    #         #     subMatchers=[
    #         #         settings_simulation_parameters,
    #         #         settings_fixes,
    #         #         settings_force_field_coefficients,
    #         #         settings_output_options,
    #         #         settings_other,
    #         #         actions
    #         #     ]
    #         # )
    #     ]
    # )
    #
    # # return SM(
    # #     name='root',
    # #     # repeats=True,
    # #     # required=True,
    # #     forwardMatch=True,
    # #     startReStr="",
    # #     fixedStartValues={'program_name': 'LAMMPS'},
    # #     sections=['section_run'],
    # #     subMatchers=[
    # #         header
    # #     ]
    # # )
    # # return SM(
    # #     name='root',
    # #     repeats=True,
    # #     required=True,
    # #     forwardMatch=True,
    # #     startReStr="",
    # #     subMatchers=[
    # #         SM(
    # #             name='newRun',
    # #             startReStr=r"",
    # #             repeats=True,
    #             # required=True,
    #             # forwardMatch=True,
    #             # sections=['section_run'],
    #             # subFlags=SM.SubFlags.Unordered,
    #             # fixedStartValues={'program_name': 'LAMMPS'},
    #             # subMatchers=[
    #             # ]),
    #         # SM(r"read_data\s(?P<x_lammps_dummy>{[file]})".format(nformat),
    #         #    sections=['section_topology'])
    #     # ])


def build_build_LammpsInputFileKeywordsSimpleMatcher():
    """Builds the list of SimpleMatchers to parse the control.in keywords of FHI-aims.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       List of SimpleMatchers that parses control.in keywords of FHI-aims.
    """
    # Now follows the list to match the keywords from the control.in.
    # Explicitly add ^ to ensure that the keyword is not within a comment.
    # Repating occurrences of the same keywords are captured.
    # List the matchers in alphabetical order according to keyword name.
    #
    return [
        # Initialization
        # Set parameters that need to be defined before atoms are created or read-in from a file.
        SM(r"(?P<x_lammps_dummy_text>\s*dimension)", repeats=True),
        SM(r"(?P<x_lammps_dummy_text>\s*atom_modify)", repeats=True),
        SM(r"(?P<x_lammps_dummy_text>\s*atom_style)", repeats=True),
        SM(r"(?P<x_lammps_dummy_text>\s*boundary)", repeats=True),
        SM(r"(?P<x_lammps_dummy_text>\s*newton)", repeats=True),
        SM(r"(?P<x_lammps_dummy_text>\s*processors)", repeats=True),

        # Units: lj, real, metal, si, cgs, electron, micro,nano
        SM(r"\s*units\s(?P<x_lammps_input_units_store>\w+)", repeats=True),

        # If force-field parameters appear in the files that will be read, these commands tell LAMMPS what kinds of force fields are being used:
        # pair_style, bond_style, angle_style, dihedral_style, improper_style.
        SM(r"\s*(?P<x_lammps_dummy_text>pair_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>bond_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>angle_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>dihedral_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>improper_style)", repeats=True),

        # Atom definition
        # There are 3 ways to define atoms in LAMMPS. Read them in from a data or restart file via the read_data or read_restart commands.
        # These files can contain molecular topology information. Or create atoms on a lattice (with no molecular topology), using these commands:
        # lattice, region, create_box, create_atoms.
        # The entire set of atoms can be duplicated to make a larger simulation using the replicate command.

        SM(r"\s*(?P<x_lammps_dummy_text>create_atoms)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>create_box)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>lattice)", repeats=True),
        SM(r"\s*read_data\s(?P<x_lammps_data_file_store>{[file]})".format(nformat), repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>read_dump)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>read_restart)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>region)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>replicate)", repeats=True),

        # Settings
        # Once atoms and molecular topology are defined, a variety of settings can be specified: force field coefficients, simulation parameters, output options, etc.
        # settings_other:

        SM(r"\s*(?P<x_lammps_dummy_text>comm_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>mass)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>set)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>velocity)", repeats=True),

        # settings_force_field_coefficients
        SM(r"\s*(?P<x_lammps_dummy_text>pair_coeff)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>bond_coeff)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>angle_coeff)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>dihedral_coeff)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>improper_coeff)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>kspace_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>kspace_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>dielectric)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>special_bonds)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>pair_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>pair_write)", repeats=True),

        # settings_simulation_parameters
        SM(r"\s*(?P<x_lammps_dummy_text>neighbor)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>neigh_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>group)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>timestep)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>reset_timestep)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>run_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>min_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>min_modify)", repeats=True),

        # settings_fixes
        SM(r"\s*(?P<x_lammps_dummy_text>fix)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>fix_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>unfix)", repeats=True),

        # settings_output_options
        SM(r"\s*(?P<x_lammps_dummy_text>dump)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>dump image)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>dump_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>dump movie)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>restart)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>thermo)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>thermo_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>thermo_style)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>undump)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>write_data)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>write_dump)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>write_restart)", repeats=True),

        # computes
        SM(r"\s*(?P<x_lammps_dummy_text>compute)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>compute_modify)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>uncompute)", repeats=True),


        # Run a simulation (actions)
        # A molecular dynamics simulation is run using the run command. Energy minimization (molecular statics) is performed using the minimize command. A parallel tempering (replica-exchange) simulation can be run using the temper command.
        SM(r"\s*(?P<x_lammps_dummy_text>delete_atoms)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>delete_bonds)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>displace_atoms)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>change_box)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>minimize)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>neb prd)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>rerun)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>run)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>temper)", repeats=True),

        # miscellaneous
        SM(r"\s*(?P<x_lammps_dummy_text>clear)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>echo)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>if)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>include)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>jump)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>label)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>log)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>next)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>print)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>shell)", repeats=True),
        SM(r"\s*(?P<x_lammps_dummy_text>variable)", repeats=True),

        #
        # SM (r"^\s*charge\s+(?P<fhi_aims_controlIn_charge>[-+0-9.eEdD]+)", repeats = True),
        # # need to distinguish different cases
        # SM (r"^\s*occupation_type\s+",
        #     forwardMatch = True,
        #     repeats = True,
        #     subMatchers = [
        #         SM (r"^\s*occupation_type\s+(?P<fhi_aims_controlIn_occupation_type>[-_a-zA-Z]+)\s+(?P<fhi_aims_controlIn_occupation_width>[-+0-9.eEdD]+)\s+(?P<fhi_aims_controlIn_occupation_order>[0-9]+)"),
        #         SM (r"^\s*occupation_type\s+(?P<fhi_aims_controlIn_occupation_type>[-_a-zA-Z]+)\s+(?P<fhi_aims_controlIn_occupation_width>[-+0-9.eEdD]+)")
        #     ]),
        # # using separate sections
        # SM (r"^\s*species\s*(?P<fhi_aims_controlIn_species_name>[a-zA-Z]+)",
        #     sections = ["fhi_aims_section_controlIn_basis_set"],
        #     repeats=True,
        #     subFlags = SM.SubFlags.Unordered,
        #     subMatchers = [
        #         SM (r"\s*\#\s*Definition of \"minimal\" basis",
        #             repeats=True,
        #             subFlags = SM.SubFlags.Unordered,
        #             subMatchers = [
        #                 SM (r"^\s*(?P<fhi_aims_controlIn_basis_func_type>[-_a-zA-Z0-9]+)"
        #                     "\s*(?P<fhi_aims_controlIn_basis_func_n>[0-9]+)"
        #                     "\s+(?P<fhi_aims_controlIn_basis_func_l>[a-zA-Z])"
        #                     "\s+(?P<fhi_aims_controlIn_basis_func_radius>[.0-9]+)",
        #                     repeats = True),
        #                 SM (r"^\s*(?P<fhi_aims_controlIn_basis_func_type>[-_a-zA-Z0-9]+)"
        #                     "\s*(?P<fhi_aims_controlIn_basis_func_n>[0-9]+)"
        #                     "\s+(?P<fhi_aims_controlIn_basis_func_l>[a-zA-Z])\s*auto",
        #                     repeats = True)
        #             ])
        #     ])
    ]

def build_LammpsInputFileSimpleMatcher2():
    """Builds the SimpleMatcher to parse the control.in file of FHI-aims.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses control.in file of FHI-aims.
    """
    return SM(
        name = 'Root1',
        startReStr = "",
        sections = ['section_run'],
        forwardMatch = True,
        weak = True,
        subMatchers = [
            SM(
                name = 'Root2',
                startReStr = "",
                sections = ['section_topology'],
                forwardMatch = True,
                weak = True,
                # The search is done unordered since the keywords do not appear in a specific order.
                subFlags = SM.SubFlags.Unordered,
                subMatchers = build_build_LammpsInputFileKeywordsSimpleMatcher()
            )
        ])



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
        'section_dos': CachingLvl,
        'section_run': CachingLvl,
        'section_single_configuration_calculation': CachingLvl,
    }
    # Set all dos metadata to Cache as they need post-processsing.
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
    LammpsInputSimpleMatcher = build_LammpsInputFileSimpleMatcher2()

    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)

    # set parser info
    parserInfo = {'name': 'lammps-input-parser', 'version': '1.0'}

    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)

    # start parsing
    mainFunction(
        mainFileDescription=LammpsInputSimpleMatcher,
        metaInfoEnv=metaInfoEnv,
        parserInfo=parserInfo,
        cachingLevelForMetaName=cachingLevelForMetaName,
        superContext=LammpsInputParserContext()
    )


if __name__ == "__main__":
    main(CachingLevel.Forward)
