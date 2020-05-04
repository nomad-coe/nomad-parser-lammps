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

from past.utils import old_div
from contextlib import contextmanager
import operator
import numpy as np


from .LammpsCommon import get_metaInfo
from .LammpsCommon import converter
from nomadcore.baseclasses import MainHierarchicalParser
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import AncillaryParser
from nomadcore.simple_parser import PushbackLineFile

from nomadcore.local_meta_info import InfoKindEnv, InfoKindEl

from . import LammpsDataParser
from . import LammpsTrjParser

nformat = {
    'float': r"[-+]?\d*\.\d+",
    'int': r"[-+]?\d+",
    'file': r"[0-9a-zA-Z_.]+"
}


################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
@contextmanager # SECTIONS ARE CLOSED AUTOMATICALLY
def open_section(p, name):
    gid = p.openSection(name)
    yield gid
    p.closeSection(name, gid)


############################################################
# This is the parser for the main LOG file of LAMMPS.
############################################################

logger = logging.getLogger(name="nomad.LammpsLogParser")


class LammpsMainParserContext(object):
    pass


class LammpsMainParser(MainHierarchicalParser):
    """The main parser class for crystal. This main parser will take the
    LAMMPS output file path as an argument.

    LAMMPS input script
    (1) LAMMPS does not read your entire input script and then perform a simulation with all the settings. Rather,
    the input script is read one line at a time and each command takes effect when it is read.
    default values ????

    Parsing rules: (source: http://lammps.sandia.gov/doc/Section_commands.html#parsing-rules)
    (1) If the last printable character on the line is a "&" character, the command is assumed to continue on the next
        line.
    (2) All characters from the first "#" character onward are treated as comment and discarded.
    (3) The line is searched repeatedly for $ characters, which indicate variables that are replaced with a text string.
    (4) The line is broken into "words" separated by whitespace (tabs, spaces). Note that words can thus contain
        letters, digits, underscores, or punctuation characters.
    (5) The first word is the command name. All successive words in the line are arguments.
    (6) If you want text with spaces to be treated as a single argument, it can be enclosed in either single or double
        or triple quotes.

    A LAMMPS input script typically has 4 parts:
    (1) Initialization
    (2) Atom definition
    (3) Settings
    (4) Run a simulation
    """

    # Define the output parsing tree for this version
    r_float = "[-+]?\d*\.\d+"  # Regex for a floating point value
    r_int = "[-+]?\d+"         # Regex for an integer
    r_file = "[0-9a-zA-Z_.]+"

    # class SubFlags:
    #     Sequenced = 0    # the subMatchers should be executed in sequence
    #     Unordered = 1    # the subMatchers can be in any order

    def __init__(self, file_path, parser_context):
        """Initialize an output parser.
        """
        super(LammpsMainParser, self).__init__(parser_context)

        # manually adjust caching of metadata
        self.caching_level_for_metaname = self.get_cachingLevelForMetaName()

        self.root_matcher = SM(
            name = 'root',
            startReStr = "",
            forwardMatch = True,
            weak = True,
            subMatchers = [
                SM(
                    name = 'new-run',
                    # startReStr="",
                    startReStr = r"(?P<program_name>\w*)\s*\((?P<program_version>.+)\)",
                    endReStr = r"clear",
                    sections = ['section_run'],
                    # fixedStartValues={'program_name': 'LAMMPS', },
                    # forwardMatch = True,
                    weak = True,
                    subMatchers = [
                        SM(
                            name = 'root',
                            startReStr = "",
                            sections = ['section_topology'],
                            forwardMatch = True,
                            weak = True,
                            # The search is done unordered since the keywords do not appear in a specific order.
                            subFlags = SM.SubFlags.Unordered,
                            subMatchers = self.build_LammpsLogFileSimpleMatcher()
                        )]
                ),
                SM(
                    name = 'new-run',
                    startReStr="",
                    # startReStr = r"(?P<program_name>\w*)\s*\((?P<program_version>.+)\)",
                    # endReStr = r"clear",
                    # sections = ['section_run'],
                    # fixedStartValues={'program_name': 'LAMMPS', },
                    # forwardMatch = True,
                    weak = True,
                    subFlags = SM.SubFlags.Unordered,
                    subMatchers = self.build_LammpsLogFileSimpleMatcher(),
                )
            ]
        )

    def get_cachingLevelForMetaName(self):
        """Sets the caching level for the metadata.

        Args:
            metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
            CachingLvl: Sets the CachingLevel for the sections dos, run, and single_configuration_calculation.
                This allows to run the parser without opening new sections.

        Returns:
            Dictionary with metaname as key and caching level as value.
        """

        caching_level_for_meta_name = {
            'x_lammps_dummy_text': CachingLevel.Cache,
            # 'section_run': CachingLevel.Cache,
            # 'section_single_configuration_calculation': CachingLevel.Cache,
        }

        # Set all _store metadata to Cache as they need post-processsing.
        for name in self.parser_context.metainfo_env.infoKinds:
            if name.endswith('_store'):
                caching_level_for_meta_name[name] = CachingLevel.Cache

        return caching_level_for_meta_name

    def build_LammpsLogFileSimpleMatcher(self):
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

            # Set the dimensionality of the simulation. By default LAMMPS runs 3d simulations.
            SM(r"(?P<x_lammps_dummy_text>\s*dimension)", ),

            SM(r"\s*boundary", forwardMatch=True, repeats=True,
               adHoc=lambda parser: parser.superContext.adHoc_boundary(parser)),

            SM(r"(?P<x_lammps_dummy_text>\s*newton)", ),
            SM(r"(?P<x_lammps_dummy_text>\s*processors)", ),

            # Units: lj, real, metal, si, cgs, electron, micro,nano
            SM(r"\s*units\s\w+", forwardMatch=True, repeats=True,
               adHoc=lambda parser: parser.superContext.adHoc_input_units(parser)),

            SM(r"(?P<x_lammps_dummy_text>\s*atom_modify)", ),
            SM(r"(?P<x_lammps_dummy_text>\s*atom_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_atom_style(parser)
               ),

            # If force-field parameters appear in the files that will be read, these commands tell LAMMPS what kinds of
            # force fields are being used:
            # pair_style, bond_style, angle_style, dihedral_style, improper_style.
            SM(r"\s*(?P<x_lammps_dummy_text>pair_style)", forwardMatch=True,  repeats=True,
               adHoc=lambda parser: self.adHoc_pair_style(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>bond_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_bond_style(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>angle_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_angle_style(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>dihedral_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_dihedral_style(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>improper_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_improper_style(parser)),

            # Atom definition
            # There are 3 ways to define atoms in LAMMPS. Read them in from a data or restart file via the read_data
            # or read_restart commands.
            # These files can contain molecular topology information. Or create atoms on a lattice (with no molecular
            # topology), using these commands:
            # lattice, region, create_box, create_atoms.
            # The entire set of atoms can be duplicated to make a larger simulation using the replicate command.

            SM(r"\s*read_data\s", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_read_data(parser),
               subFlags=SM.SubFlags.Unordered,
               subMatchers=[
                   SM(r"\s+(?P<x_lammps_dummy_text>orthogonal box = \({0} {0} {0}\) to \({0} {0} {0})".format(self.r_float)),
                   SM(r"\s+reading atoms ...",
                      subMatchers=[
                          SM(r"\s+(?P<x_lammps_dummy_text>{} atoms)".format(self.r_int)),
                      ]),
                   SM(r"\s+scanning bonds ...",
                      subMatchers=[
                          SM(r"\s+(?P<x_lammps_dummy_text>{} = max bonds/atom)".format(self.r_int)),
                      ]),
                   SM(r"\s+scanning angles ...",
                      subMatchers=[
                          SM(r"\s+(?P<x_lammps_dummy_text>{} = max angles/atom)".format(self.r_int)),
                      ]),
                   SM(r"\s+reading bonds ...",
                      subMatchers=[
                          SM(r"\s+(?P<x_lammps_dummy_text>{} bonds)".format(self.r_int)),
                      ]),
                   SM(r"\s+reading angles ...",
                       subMatchers=[
                           SM(r"\s+(?P<x_lammps_dummy_text>{} angles)".format(self.r_int)),
                       ]),
               ]),

            SM(r"\s*(?P<x_lammps_dummy_text>read_dump)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>read_restart)", ),

            SM(r"\s*(?P<x_lammps_dummy_text>lattice)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>region)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>create_atoms)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>create_box)", ),


            SM(r"\s*(?P<x_lammps_dummy_text>replicate)", ),

            # Settings
            # Once atoms and molecular topology are defined, a variety of settings can be specified: force field
            # coefficients, simulation parameters, output options, etc.
            # settings_other:

            SM(r"\s*(?P<x_lammps_dummy_text>comm_style)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>mass)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>set)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>velocity)", ),

            # settings_force_field_coefficients
            SM(r"\s*(?P<x_lammps_dummy_text>pair_coeff)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_pair_coeff(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>bond_coeff)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_bond_coeff(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>angle_coeff)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_angle_coeff(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>dihedral_coeff)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_dihedral_coeff(parser) ),
            SM(r"\s*(?P<x_lammps_dummy_text>improper_coeff)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>kspace_style)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>kspace_modify)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>dielectric)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>special_bonds)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_special_bonds(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>pair_modify)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_pair_modify(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>pair_write)", ),

            # settings_simulation_parameters
            SM(r"\s*(?P<x_lammps_dummy_text>neighbor)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>neigh_modify)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>group)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>timestep)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_timestep(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>reset_timestep)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>run_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_run_style(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>min_style)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>min_modify)", ),

            # settings_fixes
            SM(r"\s*(?P<x_lammps_dummy_text>fix)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_fix(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>fix_modify)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>unfix)", ),

            # settings_output_options
            SM(r"\s*(?P<x_lammps_dummy_text>dump)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_dump(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>dump image)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>dump_modify)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>dump movie)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>restart)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>thermo)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_thermo(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>thermo_modify)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>thermo_style)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_thermo_style(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>undump)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>write_data)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>write_dump)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>write_restart)", ),

            # computes
            SM(r"\s*(?P<x_lammps_dummy_text>compute)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>compute_modify)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>uncompute)", ),


            # Run a simulation (actions)
            # A molecular dynamics simulation is run using the run command. Energy minimization (molecular statics)
            # is performed using the minimize command. A parallel tempering (replica-exchange) simulation can be run
            # using the temper command.
            SM(r"\s*(?P<x_lammps_dummy_text>delete_atoms)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>delete_bonds)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>displace_atoms)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>change_box)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>minimize)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>neb prd)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>rerun)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>run)", forwardMatch=True,
               adHoc=lambda parser: self.adHoc_run(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>temper)", ),

            # miscellaneous
            SM(r"\s*(?P<x_lammps_dummy_text>clear)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>echo)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>if)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>include)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>jump)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>label)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>log)", forwardMatch=True, repeats=True,
               adHoc=lambda parser: self.adHoc_log(parser)),
            SM(r"\s*(?P<x_lammps_dummy_text>next)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>print)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>shell)", ),
            SM(r"\s*(?P<x_lammps_dummy_text>variable)", ),

        ]

    def parse(self, mainfile):
        """Starts the parsing. By default uses the SimpleParser scheme, if you
        want to use something else or customize the process just override this
        method in the subclass.
        """
        mainFunction(
            mainFileDescription=self.root_matcher,
            metaInfoEnv=self.parser_context.metainfo_env,
            parserInfo=self.parser_context.parser_info,
            # outF=self.parser_context.super_backend.fileOut,
            cachingLevelForMetaName=self.caching_level_for_metaname,
            superContext=self,
            # onClose=self.onClose,
            default_units=self.parser_context.default_units,
            metainfo_units=self.parser_context.metainfo_units,
            superBackend=self.parser_context.super_backend,
            metaInfoToKeep=self.parser_context.metainfo_to_keep,
            mainFile=mainfile
        )
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
        self.converter = converter()
        self.boundary = 3 * [True]
        self.log_files = []
        self.specs_filter = []
        self.pm_filter = []
        self.sb_filter = []
        self.bond_coeff = []
        self.angle_coeff = []
        self.dihedral_coeff = []
        self.pair_coeff = []
        self.run = []
        self.thermo = []
        self.thermo_style = []
        self.fix = []
        self.run_style = []
        self.run = []
        self.timestep = []
        self.dump = []
        self.read_data = []
        self.data_file = None

    def compile_log_parser(self):
        """Instantiate superContext and construct parser for external data file.
        """
        self.logParser = AncillaryParser(
            fileDescription=self.build_LammpsLogFileSimpleMatcher(),
            parser=self.parser,
            cachingLevelForMetaName=self.get_cachingLevelForMetaName(),
            superContext=self
        )


    def compile_data_parser(self):
        """Instantiate superContext and construct parser for external data file.
        """
        self.dataSuperContext = LammpsDataParser.LammpsDataParserContext(converter=self.converter, writeMetaData=False)
        self.dataParser = AncillaryParser(
            fileDescription=LammpsDataParser.build_LammpsDataFileSimpleMatcher(),
            parser=self.parser,
            cachingLevelForMetaName=LammpsDataParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext=self.dataSuperContext
        )

    def compile_trj_parser(self):
        """Instantiate superContext and construct parser for external data file.
        """
        self.trjSuperContext = LammpsTrjParser.LammpsTrjParserContext(converter=self.converter, dump_style="", writeMetaData=False)
        self.trjParser = AncillaryParser(
            fileDescription=LammpsTrjParser.build_LammpsTrjFileSimpleMatcher(),
            parser=self.parser,
            cachingLevelForMetaName=LammpsTrjParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext=self.trjSuperContext
        )



    #===========================================================================
    # The functions that trigger when sections are closed


    def onClose_section_topology(self, backend, gIndex, section):
        """Trigger called when section_topology is closed.

        Input data is parsed from external file but the result is only stored and written later in lammps_topology
        in section_run.
        """

        # if section['x_lammps_input_units_store'] is not None:
        #     self.units = section['x_lammps_input_units_store']

        p = self.super_backend
        o = open_section


        # # collecting list of force field functional styles
        list_of_styles = self.readStyles()
        pairFunctionalAndCutOff = list_of_styles.get('pair_style')  ### pair interactions style LAMMPS name + cutoff
        pairFunctional = pairFunctionalAndCutOff[0]                 ### pair interactions style LAMMPS name only
        bondFunctional = list_of_styles.get('bond_style')
        angleFunctional = list_of_styles.get('angle_style')
        dihedralFunctional = list_of_styles.get('dihedral_style')
        # ###

        # read the trajectory from the data file if it is presented
        if self.data_file is not None:
            self.compile_data_parser()

            dir_name = os.path.dirname(os.path.abspath(self.fName))
            f_name = os.path.normpath(os.path.join(dir_name, self.data_file))

            try:
                with open(f_name) as fIn:
                    # construct parser for DATA file if not present
                    if self.dataSuperContext is None or self.dataParser is None:
                        self.compile_data_parser()
                    # parse Data file
                    self.dataParser.parseFile(fIn)
            except IOError:
                logger.error("DATA file parsing unsuccessful. Could not find %s file in directory '%s'." % (self.data_file, dir_name))


            list_of_bonds = self.readBonds(bondFunctional)
            list_of_bonds = sorted(list(list_of_bonds.items()), key=operator.itemgetter(0))

            list_of_angles = self.readAngles(angleFunctional)
            list_of_angles = sorted(list(list_of_angles.items()), key=operator.itemgetter(0))

            list_of_dihedrals = self.readDihedrals(dihedralFunctional)
            list_of_dihedrals = sorted(list(list_of_dihedrals.items()), key=operator.itemgetter(0))

            # collecting dispersion interactions ff terms
            updateAtomTypes= self.dataSuperContext.updateAtomTypes
            list_of_ljs, ljs_dict  = self.readPairCoeff(updateAtomTypes, pairFunctional)
            ljs_dict = sorted(list(ljs_dict.items()), key=operator.itemgetter(0))
            list_of_ljs = sorted(list(list_of_ljs.items()), key=operator.itemgetter(0))
            lj_types = len(ljs_dict)

            ####################################################################################################################################################################################################################################
            # FF DEFINITION CAN BE SET WITHIN THE LAMMPS TOPOLOGY (DATA) FILE, here we collect the FF parameters from the .data file
            # (if they are not provided in the input file)

            if not list_of_bonds:
                list_of_bonds = self.dataSuperContext.readBondsFromData(self.dataSuperContext.bond_list, bondFunctional)
                list_of_bonds = sorted(list(list_of_bonds.items()), key=operator.itemgetter(0))


            if not list_of_angles:
                list_of_angles = self.dataSuperContext.readAnglesFromData(self.dataSuperContext.angle_list, angleFunctional)
                list_of_angles = sorted(list(list_of_angles.items()), key=operator.itemgetter(0))


            if not list_of_dihedrals:
                list_of_dihedrals = self.dataSuperContext.readDihedralsFromData(self.dataSuperContext.dihedral_list, dihedralFunctional)
                list_of_dihedrals = sorted(list(list_of_dihedrals.items()), key=operator.itemgetter(0))


            #########################################
            ############### missing
            #########################################
            if not list_of_ljs:
                list_of_ljs, ljs_dict = self.dataSuperContext.readPairCoeffFromData(self.dataSuperContext.lj_list, self.dataSuperContext.updateAtomTypes, pairFunctional)
                ljs_dict = sorted(list(ljs_dict.items()), key=operator.itemgetter(0))
                list_of_ljs = sorted(list(list_of_ljs.items()), key=operator.itemgetter(0))
                lj_types = len(ljs_dict)



            #### BASIC TOPOLOGY INFORMATION IN section_topology ################################################################################################################################################################################
            ####################################################################################################################################################################################################################################
            # number_of_topology_atoms = numberOfTopologyAtoms()
            # p.addValue('number_of_topology_atoms', number_of_topology_atoms)  # backend add number of topology atoms

            number_of_topology_atoms = self.dataSuperContext.at_count
            atom_to_molecule = []
            for i in range(number_of_topology_atoms):
                atom_to_molecule.append([self.dataSuperContext.moleculeInfoResolved[i][1],
                                         self.dataSuperContext.moleculeInfoResolved[i][3]])

            atomic_number = []
            for i in range(number_of_topology_atoms):
                pass

            p.addValue('number_of_topology_molecules', len(self.dataSuperContext.moleculeInfo))  # backend add number of topology molecules
            p.addArrayValues('atom_to_molecule', np.asarray(atom_to_molecule))


            #### ATOM TYPE INFORMATION IN section_atom_type ####################################################################################################################################################################################
            ####################################################################################################################################################################################################################################
            at_types = len(self.dataSuperContext.mass_list)
            for i in range(at_types):

                with o(p, 'section_atom_type'):
                    # p.addValue('atom_type_name', [mass_xyz[i], i+1])  # Here atom_type_name is atomic number plus an integer index identifying the atom type
                    p.addValue('atom_type_name', str(self.dataSuperContext.mass_xyz[i])+' : '+str(i+1))  ### TO BE CHECKED LATER
                    p.addValue('atom_type_mass', self.converter.Mass(self.dataSuperContext.mass_list[i][1]))     # Atomic mass
                    p.addValue('atom_type_charge', self.converter.Charge(self.dataSuperContext.charge_dict[i][1])) # Atomic charge, either partial or ionic
                    pass

                    ####################################################################################################################################################################################################################################

            ####################################################################################################################################################################################################################################
            #### COVALENT BONDS INFORMATION IN section_interaction (number_of_atoms_per_interaction = 2) #######################################################################################################################################
            ####################################################################################################################################################################################################################################

            if self.dataSuperContext.bd_types:  #  COVALENT BONDS

                store = []
                interaction_atoms = []
                for i in self.dataSuperContext.bondTypeList:
                    for j in self.dataSuperContext.bondTypeList:

                        store = [ [x[1], x[2]] for x in self.dataSuperContext.bond_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)


                for i in range(len(self.dataSuperContext.bondTypeList)):

                    with o(p, 'section_interaction'):
                        p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))     # atom indexes of bound pairs for a specific covalent bond
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))             # number of covalent bonds of this type
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0])) # number of atoms involved (2 for covalent bonds)

                        if bondFunctional:
                            p.addValue('interaction_kind', bondFunctional)  # functional form of the interaction

                        int_index_store = self.dataSuperContext.bond_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                        else:
                            for line in int_index_store:
                                temp = sorted([x-1 for x in line])
                                interaction_atom_to_atom_type_ref.append(temp)

                        bondParameters = dict()
                        bondParameters.update({list_of_bonds[i][0] : list_of_bonds[i][1]})

                        p.addArrayValues('x_lammps_interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', bondParameters)  # interaction parameters for the functional


            ####################################################################################################################################################################################################################################
            #### BOND ANGLES INFORMATION IN section_interaction (number_of_atoms_per_interaction = 3) ##########################################################################################################################################
            ####################################################################################################################################################################################################################################
            ag_types = self.dataSuperContext.ag_types
            if ag_types:  # BOND ANGLES

                store = []
                interaction_atoms = []
                for i in self.dataSuperContext.angleTypeList:
                    for j in self.dataSuperContext.angleTypeList:

                        store = [ [x[1], x[2], x[3]] for x in self.dataSuperContext.angle_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(self.dataSuperContext.angleTypeList)):

                    with o(p, 'section_interaction'):
                        p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))     # atom indexes of triplets for a specific bond angle
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))             # number of bond angles of this type
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0])) # number of atoms involved (3 for bond angles)

                        if bondFunctional:
                            p.addValue('interaction_kind', angleFunctional)  # functional form of the interaction

                        int_index_store = self.dataSuperContext.angle_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1]

                        else:
                            for line in int_index_store:
                                temp = [x-1 for x in line]
                                interaction_atom_to_atom_type_ref.append(temp)

                        angleParameters = dict()
                        angleParameters.update({list_of_angles[i][0] : list_of_angles[i][1]})

                        p.addArrayValues('x_lammps_interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', angleParameters)  # interaction parameters for the functional

            ####################################################################################################################################################################################################################################
            #### DIHEDRAL ANGLES INFORMATION IN section_interaction (number_of_atoms_per_interaction = 4) ######################################################################################################################################
            ####################################################################################################################################################################################################################################

            if self.dataSuperContext.dh_types:  # DIHEDRAL ANGLES

                store = []
                interaction_atoms = []
                for i in self.dataSuperContext.dihedralTypeList:
                    for j in self.dataSuperContext.dihedralTypeList:

                        store = [ [x[1], x[2], x[3], x[4]] for x in self.dataSuperContext.dihedral_interaction_atoms if x[0]==i ]
                    interaction_atoms.append(store)

                for i in range(len(self.dataSuperContext.dihedralTypeList)):

                    with o(p, 'section_interaction'):
                        p.addArrayValues('interaction_atoms', np.asarray(interaction_atoms[i]))     # atom indexes of quartets for a specific dihedral angle
                        p.addValue('number_of_interactions', len(interaction_atoms[i]))             # number of dihedral angles of this type
                        p.addValue('number_of_atoms_per_interaction', len(interaction_atoms[0][0])) # number of atoms involved (4 for dihedral angles)

                        if bondFunctional:
                            p.addValue('interaction_kind', dihedralFunctional)  # functional form of the interaction

                        int_index_store = self.dataSuperContext.dihedral_dict[i][1]
                        interaction_atom_to_atom_type_ref = []

                        if all(isinstance(elem, list) for elem in int_index_store) == False:
                            interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1, int_index_store[3]-1]

                        else:
                            for line in int_index_store:
                                temp = [x-1 for x in line]
                                interaction_atom_to_atom_type_ref.append(temp)

                        dihedralParameters = dict()
                        dihedralParameters.update({list_of_dihedrals[i][0] : list_of_dihedrals[i][1]})
                        p.addArrayValues('x_lammps_interaction_atom_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                        p.addValue('interaction_parameters', dihedralParameters)  # interaction parameters for the functional

            ####################################################################################################################################################################################################################################
            #### DISPERSIVE FORCES INFORMATION IN section_interaction (number_of_atoms_per_interaction = 2) ####################################################################################################################################
            ####################################################################################################################################################################################################################################

            if lj_types:  # LJ-like interactions

                with o(p, 'section_interaction'):

                    p.addValue('x_lammps_number_of_defined_pair_interactions', lj_types)  # number of LJ interaction types
                    p.addValue('number_of_atoms_per_interaction', len(ljs_dict[0][1]))  # = 2 for pair interactions

                    if pairFunctionalAndCutOff:
                        p.addValue('interaction_kind', str(pairFunctionalAndCutOff))  # functional form of the interaction (cutoff radius included)

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
                            temp = [x-1 for x in line]
                            interaction_atom_to_atom_type_ref.append(temp)

                    # p.addValue('interaction_atoms', int_index_store)
                    p.addArrayValues('x_lammps_pair_interaction_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                    p.addArrayValues('x_lammps_pair_interaction_parameters', np.asarray(int_param_store))  # interaction parameters for the functional
                    pass

            ####################################################################################################################################################################################################################################
            #### MOLECULE TYPE INFORMATION IN section_molecule_type ############################################################################################################################################################################
            ####################################################################################################################################################################################################################################
            moleculeTypeInfo = self.dataSuperContext.moleculeTypeInfo
            for i in range(len(moleculeTypeInfo)):

                with o(p, 'section_molecule_type'):
                    # gindex = 0
                    p.addValue('molecule_type_name', 'molecule'+'_'+str(moleculeTypeInfo[i][0]))
                    p.addValue('number_of_atoms_in_molecule', len(moleculeTypeInfo[i][1]))

                    p.addArrayValues('atom_in_molecule_to_atom_type_ref', np.asarray([x-1 for x in moleculeTypeInfo[i][1]]))


                    atom_in_molecule_name = []
                    for j in moleculeTypeInfo[i][1]:
                        atom_in_molecule_name.append([ self.dataSuperContext.mass_xyz[j-1], j ] ) # Here atom_in_molecule_name is atomic number plus an integer index

                    # p.addArrayValues('atom_in_molecule_name', np.asarray(atom_in_molecule_name))

                    atom_in_molecule_charge = []
                    for j in moleculeTypeInfo[i][1]:
                        atom_in_molecule_charge.append(self.dataSuperContext.charge_list[j-1][1])

                    p.addValue('atom_in_molecule_charge', atom_in_molecule_charge)

                    ############################################################################################################################################################################################################################
                    #### COVALENT BONDS INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 2) ######################################################################################################################

                    if self.dataSuperContext.bd_types:

                        toMoleculeAtomIndex  = min( moleculeTypeInfo[i][2] )

                        store = []
                        molecule_interaction_atoms = []
                        molecule_interaction_type  = []
                        for h in self.dataSuperContext.bondTypeList:
                            for k in self.dataSuperContext.bondTypeList:

                                store   = [ [x[1] - toMoleculeAtomIndex, x[2] - toMoleculeAtomIndex] for x in self.dataSuperContext.bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]
                                # store   = [ [x[1], x[2]] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]

                            molecule_interaction_atoms.append(store)

                        for l in self.dataSuperContext.bondTypeList:
                            store1  = [ x[0] for x in self.dataSuperContext.bond_interaction_atoms if x[0]==l and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]
                            molecule_interaction_type.append(store1)

                        molecule_interaction_type = [ x for sublist in molecule_interaction_type for x in sublist ] # from list of lists to list

                        for bond in self.dataSuperContext.bondTypeList:
                            if bond in molecule_interaction_type:

                                with o(p, 'section_molecule_interaction'):

                                    p.addArrayValues('molecule_interaction_atoms', np.asarray(molecule_interaction_atoms[bond-1]))
                                    p.addValue('number_of_molecule_interactions', len(molecule_interaction_atoms[bond-1]))
                                    p.addValue('number_of_atoms_per_molecule_interaction', len(molecule_interaction_atoms[0][0]))

                                    if bondFunctional:
                                        p.addValue('molecule_interaction_kind', bondFunctional)


                                    int_index_store = self.dataSuperContext.bond_dict[bond-1][1]

                                    molecule_interaction_atom_to_atom_type_ref = []
                                    if all(isinstance(elem, list) for elem in int_index_store) == False:
                                        molecule_interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1]

                                    else:
                                        for line in int_index_store:
                                            temp = sorted([x-1 for x in line])
                                            molecule_interaction_atom_to_atom_type_ref.append(temp)

                                    moleculeBondParameters = dict()
                                    moleculeBondParameters.update({list_of_bonds[bond-1][0] : list_of_bonds[bond-1][1]})
                                    p.addArrayValues('x_lammps_molecule_interaction_atom_to_atom_type_ref', np.asarray(molecule_interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                                    p.addValue('molecule_interaction_parameters', moleculeBondParameters)

                    ############################################################################################################################################################################################################################
                    #### BOND ANGLE INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 3) #########################################################################################################################

                    if ag_types:

                        toMoleculeAtomIndex  = min( moleculeTypeInfo[i][2] )

                        store = []
                        molecule_interaction_atoms = []
                        molecule_interaction_type  = []
                        for h in self.dataSuperContext.angleTypeList:
                            for k in self.dataSuperContext.angleTypeList:

                                store   = [ [x[1] - toMoleculeAtomIndex, x[2] - toMoleculeAtomIndex, x[3] - toMoleculeAtomIndex] for x in self.dataSuperContext.angle_interaction_atoms
                                            if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] and x[3] in moleculeTypeInfo[i][2] ]
                                # store   = [ [x[1], x[2], x[3]] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]

                            molecule_interaction_atoms.append(store)

                        for l in self.dataSuperContext.angleTypeList:
                            store1  = [ x[0] for x in self.dataSuperContext.angle_interaction_atoms if x[0]==l and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]
                            molecule_interaction_type.append(store1)

                        molecule_interaction_type = [ x for sublist in molecule_interaction_type for x in sublist ] # from list of lists to list

                        for angle in self.dataSuperContext.angleTypeList:
                            if angle in molecule_interaction_type:

                                with o(p, 'section_molecule_interaction'):

                                    p.addArrayValues('molecule_interaction_atoms', np.asarray(molecule_interaction_atoms[angle-1]))
                                    p.addValue('number_of_molecule_interactions', len(molecule_interaction_atoms[angle-1]))
                                    p.addValue('number_of_atoms_per_molecule_interaction', len(molecule_interaction_atoms[0][0]))

                                    if bondFunctional:
                                        p.addValue('molecule_interaction_kind', angleFunctional)


                                    int_index_store = self.dataSuperContext.angle_dict[angle-1][1]
                                    molecule_interaction_atom_to_atom_type_ref = []

                                    # print int_index_store, '######'

                                    if all(isinstance(elem, list) for elem in int_index_store) == False:
                                        molecule_interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1]

                                    else:
                                        for line in int_index_store:
                                            temp = [x-1 for x in line]
                                            molecule_interaction_atom_to_atom_type_ref.append(temp)

                                    moleculeAngleParameters = dict()
                                    moleculeAngleParameters.update({list_of_angles[angle-1][0] : list_of_angles[angle-1][1]})
                                    p.addArrayValues('x_lammps_molecule_interaction_atom_to_atom_type_ref', np.asarray(molecule_interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                                    p.addValue('molecule_interaction_parameters', moleculeAngleParameters)

                    ############################################################################################################################################################################################################################
                    #### DIHEDRAL ANGLE INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 4) ######################################################################################################################

                    if self.dataSuperContext.dh_types:

                        toMoleculeAtomIndex  = min( moleculeTypeInfo[i][2] )

                        store = []
                        molecule_interaction_atoms = []
                        molecule_interaction_type  = []
                        for h in self.dataSuperContext.dihedralTypeList:
                            for k in self.dataSuperContext.dihedralTypeList:

                                store   = [ [x[1] - toMoleculeAtomIndex, x[2] - toMoleculeAtomIndex, x[3] - toMoleculeAtomIndex, x[4] - toMoleculeAtomIndex]
                                            for x in self.dataSuperContext.dihedral_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2]
                                            and x[3] in moleculeTypeInfo[i][2] and x[4] in moleculeTypeInfo[i][2] ]
                                # store   = [ [x[1], x[2], x[3]] for x in bond_interaction_atoms if x[0]==h and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] ]

                            molecule_interaction_atoms.append(store)

                        for l in self.dataSuperContext.dihedralTypeList:
                            store1  = [ x[0] for x in self.dataSuperContext.dihedral_interaction_atoms if x[0]==l
                                        and x[1] in moleculeTypeInfo[i][2] and x[2] in moleculeTypeInfo[i][2] and x[3] in moleculeTypeInfo[i][2] and x[4] in moleculeTypeInfo[i][2] ]
                            molecule_interaction_type.append(store1)

                        molecule_interaction_type = [ x for sublist in molecule_interaction_type for x in sublist ] # from list of lists to list

                        # print molecule_interaction_type,  '#######'
                        # print molecule_interaction_atoms, '#######'

                        for dihedral in self.dataSuperContext.dihedralTypeList:
                            if dihedral in molecule_interaction_type:

                                with o(p, 'section_molecule_interaction'):

                                    p.addArrayValues('molecule_interaction_atoms', np.asarray(molecule_interaction_atoms[dihedral-1]))
                                    p.addValue('number_of_molecule_interactions', len(molecule_interaction_atoms[dihedral-1]))
                                    p.addValue('number_of_atoms_per_molecule_interaction', len(molecule_interaction_atoms[0][0]))

                                    if bondFunctional:
                                        p.addValue('molecule_interaction_kind', dihedralFunctional)


                                    int_index_store = self.dataSuperContext.dihedral_dict[dihedral-1][1]
                                    molecule_interaction_atom_to_atom_type_ref = []

                                    # print int_index_store, '######'

                                    if all(isinstance(elem, list) for elem in int_index_store) == False:
                                        molecule_interaction_atom_to_atom_type_ref = [int_index_store[0]-1, int_index_store[1]-1, int_index_store[2]-1, int_index_store[3]-1]

                                    else:
                                        for line in int_index_store:
                                            temp = [x-1 for x in line]
                                            molecule_interaction_atom_to_atom_type_ref.append(temp)

                                    moleculeDihedralParameters = dict()
                                    moleculeDihedralParameters.update({list_of_dihedrals[dihedral-1][0] : list_of_dihedrals[dihedral-1][1]})
                                    p.addArrayValues('x_lammps_molecule_interaction_atom_to_atom_type_ref', np.asarray(molecule_interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                                    p.addValue('molecule_interaction_parameters', moleculeDihedralParameters)

                    ############################################################################################################################################################################################################################
                    #### DISPERSIVE FORCES INFORMATION IN section_molecule_interaction (number_of_atoms_per_interaction = 2) ###################################################################################################################

                    if lj_types:  # LJ-like interactions

                        with o(p, 'section_molecule_interaction'):

                            if pairFunctionalAndCutOff:
                                p.addValue('molecule_interaction_kind', str(pairFunctionalAndCutOff))  # functional form of the interaction (cutoff radous included)

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
                                    temp = [x-1 for x in line]
                                    interaction_atom_to_atom_type_ref.append(temp)

                            p.addValue('x_lammps_number_of_defined_molecule_pair_interactions', lj_types)  # number of LJ interaction types
                            p.addValue('number_of_atoms_per_molecule_interaction', len(ljs_dict[0][1]))  # = 2 for pair interactions
                            p.addArrayValues('x_lammps_pair_molecule_interaction_to_atom_type_ref', np.asarray(interaction_atom_to_atom_type_ref))  # this points to the relative section_atom_type
                            p.addArrayValues('x_lammps_pair_molecule_interaction_parameters', np.asarray(int_param_store))  # interaction parameters for the functional

                        pass

            ############################################################################################################################################################################################################################
            #### MAPPING THE TOPOLOGY MOLECULES TO THE RELATIVE section_molecule type ##########################################################################################################################################################

            molecule_to_molecule_type_map = []
            for i in range(len(self.dataSuperContext.moleculeInfo)):
                molecule_to_molecule_type_map.append(self.dataSuperContext.moleculeInfo[i][1]-1) # mapping molecules to the relative section_molecule_type

            p.addArrayValues('molecule_to_molecule_type_map', np.asarray(molecule_to_molecule_type_map))

            ####################################################################################################################################################################################################################################



    def onClose_section_run(self, backend, gIndex, section):
        """Trigger called when section_topology is closed.

        Input data is parsed from external file but the result is only stored and written later in lammps_topology
        in section_run.
        """

        p = self.super_backend
        o = open_section


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

        # tmp = section['x_lammps_dummy_text']

        # if section['x_lammps_input_units_store'] is not None:
        #     self.units = section['x_lammps_input_units_store']

        # construct file name
        # data_file = section['x_lammps_data_file_store']

        ####################################################################################################################################################################################################################################
        #### BASIC SAMPLING INFORMATION IN section_topology ################################################################################################################################################################################
        ####################################################################################################################################################################################################################################

        refToSamp = []
        with o(p, 'section_sampling_method') as gid:

            refToSamp.append(gid)

            ensemble, sampling = self.readEnsemble()
            target_t, thermo_tau, langevin_gamma, target_p, baro_tau = self.readTPSettings()
            int_type, tstep, steps = self.readIntegratorSettings()

            p.addValue('x_lammps_integrator_type', int_type)
            p.addValue('x_lammps_integrator_dt', tstep*toTime)
            p.addValue('x_lammps_number_of_steps_requested', steps)

            p.addValue('sampling_method', sampling)
            p.addValue('ensemble_type', ensemble)

            if target_t:
                p.addValue('x_lammps_thermostat_target_temperature', target_t*toTemp)
                p.addValue('x_lammps_thermostat_tau', thermo_tau*toTime)

            if target_p:
                p.addValue('x_lammps_barostat_target_pressure', target_p*toPress)
                p.addValue('x_lammps_barostat_tau', baro_tau*toTime)

            if langevin_gamma:
                p.addValue('x_lammps_thermostat_target_temperature', target_t*toTemp)
                p.addValue('x_lammps_langevin_gamma', langevin_gamma*toTime)

        ########################################################################################################################################################################################################################################
        #### SYSTEM INFORMATION TO section_system ##############################################################################################################################################################################################
        ########################################################################################################################################################################################################################################

        # skipTraj = trajFileOpen()     # if False no trajectory info is available here
        # fNameTraj, stepsPrintFrame, trajDumpStyle = readDumpFileName()
        #
        #
        # #### INTIAL CONFIGURATION AND ATOM LABELS
        #
        # atomPosInit = []
        # atomAtLabel  = []
        # for i in range(len(atomLabelling)):
        #     storeAtNumb = atomLabelling[i][0]    # atomic number from data file
        #     storeAtPos  = [atomLabelling[i][1], atomLabelling[i][2], atomLabelling[i][3]] # atomic position from data file
        #     atomAtLabel.append(storeAtNumb)
        #     atomPosInit.append(storeAtPos)
        #
        #
        # H  = 0
        # C  = 0
        # N  = 0
        # O  = 0
        # P  = 0
        # S  = 0
        # He = 0
        # Ne = 0
        # Ar = 0
        # B  = 0
        # F  = 0
        # Cl = 0
        # Br = 0
        # Na = 0
        # K  = 0
        # Al = 0
        # Si = 0
        # Ti = 0
        # V  = 0
        # Cr = 0
        # Mn = 0
        # Fe = 0
        # Co = 0
        # Ni = 0
        # Cu = 0
        # Zn = 0
        # Ag = 0
        # Au = 0
        #
        # for i, at in enumerate(atomAtLabel):  # converting atomic number to atom_labels
        #     if at == 1:
        #         H += 1
        #         atomAtLabel[i] = 'H' + '  ' + str(H) # hydrogen
        #     if at == 6:
        #         C += 1
        #         atomAtLabel[i] = 'C' + '  ' + str(C) # carbon
        #     if at == 7:
        #         N += 1
        #         atomAtLabel[i] = 'N' + '  ' + str(N) # nitrogen
        #     if at == 8:
        #         O += 1
        #         atomAtLabel[i] = 'O' + '  ' + str(O) # oxygen
        #     if at == 9:
        #         P += 1
        #         atomAtLabel[i] = 'P' + '  ' + str(P) # phosphorus
        #     if at == 16:
        #         S += 1
        #         atomAtLabel[i] = 'S' + '  ' + str(S) # sulfur
        #     if at == 2:
        #         He += 1
        #         atomAtLabel[i] = 'He' + '  ' + str(He) # helium
        #     if at == 10:
        #         Ne += 1
        #         atomAtLabel[i] = 'Ne' + '  ' + str(Ne) # neon
        #     if at == 18:
        #         Ar += 1
        #         atomAtLabel[i] = 'Ar' + '  ' + str(Ar) # argon
        #     if at == 5:
        #         B += 1
        #         atomAtLabel[i] = 'B' + '  ' + str(B) # boron
        #     if at == 9:
        #         F += 1
        #         atomAtLabel[i] = 'F' + '  ' + str(F) # fluorine
        #     if at == 17:
        #         Cl += 1
        #         atomAtLabel[i] = 'Cl' + '  ' + str(Cl) # chlorine
        #     if at == 35:
        #         Br += 1
        #         atomAtLabel[i] = 'Br' + '  ' + str(Br) # bromine
        #     if at == 11:
        #         Na += 1
        #         atomAtLabel[i] = 'Na' + '  ' + str(Na) # sodium
        #     if at == 19:
        #         K += 1
        #         atomAtLabel[i] = 'K' + '  ' + str(K) # potassium
        #     if at == 13:
        #         Al += 1
        #         atomAtLabel[i] = 'Al' + '  ' + str(Al) # aluminium
        #     if at == 14:
        #         Si += 1
        #         atomAtLabel[i] = 'Si' + '  ' + str(Si) # silicon
        #     if at == 22:
        #         Ti += 1
        #         atomAtLabel[i] = 'Ti' + '  ' + str(Ti) # titanium
        #     if at == 23:
        #         V += 1
        #         atomAtLabel[i] = 'V' + '  ' + str(V) # vanadium
        #     if at == 24:
        #         Cr += 1
        #         atomAtLabel[i] = 'Cr' + '  ' + str(Cr) # chromium
        #     if at == 25:
        #         Mn += 1
        #         atomAtLabel[i] = 'Mn' + '  ' + str(Mn) # manganese
        #     if at == 26:
        #         Fe += 1
        #         atomAtLabel[i] = 'Fe' + '  ' + str(Fe) # iron
        #     if at == 27:
        #         Co += 1
        #         atomAtLabel[i] = 'Co' + '  ' + str(Co) # cobalt
        #     if at == 28:
        #         Ni += 1
        #         atomAtLabel[i] = 'Ni' + '  ' + str(Ni) # nickel
        #     if at == 29:
        #         Cu += 1
        #         atomAtLabel[i] = 'Cu' + '  ' + str(Cu) # copper
        #     if at == 30:
        #         Zn += 1
        #         atomAtLabel[i] = 'Zn' + '  ' + str(Zn) # zinc
        #     if at == 47:
        #         Ag += 1
        #         atomAtLabel[i] = 'Ag' + '  ' + str(Ag) # silver
        #     if at == 79:
        #         Au += 1
        #         atomAtLabel[i] = 'Au' + '  ' + str(Au) # gold
        #
        #
        # stepsPrintFrame = 0
        # fNameTraj, stepsPrintFrame, trajDumpStyle = readDumpFileName()
        # frame_length, simulation_length, stepsPrintThermo, integrationSteps = simulationTime()

        # read the trajectory from the lammpstrj file if it is presented
        self.dump_file = self.dump[0].split()[5]
        if self.dump_file is not None:
            self.compile_trj_parser()

            dir_name = os.path.dirname(os.path.abspath(self.fName))
            f_name = os.path.normpath(os.path.join(dir_name, self.dump_file))

            try:
                with open(f_name) as fIn:
                    # construct parser for TRJ file if not present
                    if self.trjSuperContext is None or self.trjParser is None:
                        self.compile_trj_parser()
                    # parse TRJ file
                    self.trjParser.parseFile(fIn)
                    pass
            except IOError:
                logger.error("TRJ file parsing unsuccessful. Could not find %s file in directory '%s' (%s)." % (self.dump_file, dir_name, f_name))


        # make sure that all file is closed properly
        if not self.parser.fIn.fIn.closed:
            self.parser.fIn.fIn.close()

    #===========================================================================
    # adHoc functions that are used to do custom parsing. Primarily these
    # functions are used for data that is formatted as a table or a list.

    # @staticmethod
    # def adHoc_x_lamps_dummy(parser):
    #     line1 = parser.fIn.readline()

    def adHoc_dummy(self, parser):
        line = parser.fIn.readline()

        return None
        # metaName = 'x_lammps_data_file_store'
        # parser.backend.addValue(metaName, line.split()[1])
        # section = self.getSection(parser)

    def adHoc_input_units(self, parser):
        line = parser.fIn.readline()
        self.converter = converter(line.split()[1])

        return None

    def adHoc_boundary(self, parser):
        line = parser.fIn.readline()

        self.boundary = list(map(lambda x: x == "p", line.split()[1:]))

        # metaName = 'configuration_periodic_dimensions'
        # backend.addArray(metaName, self.boundary)

        return None

    def adHoc_read_data(self, parser):
        line = parser.fIn.readline()

        self.data_file = line.split()[1]
        self.read_data.append(line)

        # metaName = 'data_file'
        # backend.addArray(metaName, self.boundary)

        return None

    def adHoc_log(self, parser):
        '''Change the file reader according to the log file
        '''

        # self.compile_log_parser()

        line = parser.fIn.readline()
        filename = line.split()[1]

        dir_name = os.path.dirname(os.path.abspath(self.fName))
        f_name = os.path.normpath(os.path.join(dir_name, filename))

        try:

            fIn = open(f_name, 'r')

        except IOError:
            logger.error("LOG file parsing unsuccessful. Could not find %s file in directory '%s' (%s)." % (f_name, dir_name, f_name))
            return None

        # close the file reader of the actual log file
        parser.fIn.fIn.close()

        # owerwrite and open the file reader of the new log file
        parser.fIn = PushbackLineFile(fIn)




        # dir_name = os.path.dirname(os.path.abspath(self.fName))
        # f_name = os.path.normpath(os.path.join(dir_name, filename))



        return None


    def adHoc_atom_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None

    def adHoc_pair_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None

    def adHoc_bond_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None

    def adHoc_angle_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None

    def adHoc_dihedral_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None
    def adHoc_improper_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None


    def adHoc_thermo_style(self, parser):
        line = parser.fIn.readline()
        self.specs_filter.append(line)
        return None

    def adHoc_pair_modify(self, parser):
        line = parser.fIn.readline()
        self.pm_filter.append(line)
        return None

    def adHoc_special_bonds(self, parser):
        line = parser.fIn.readline()
        self.sb_filter.append(line)
        return None

    def adHoc_bond_coeff(self, parser):
        line = parser.fIn.readline()
        self.bond_coeff.append(line)
        return None

    def adHoc_angle_coeff(self, parser):
        line = parser.fIn.readline()
        self.angle_coeff.append(line)
        return None

    def adHoc_dihedral_coeff(self, parser):
        line = parser.fIn.readline()
        self.dihedral_coeff.append(line)
        return None

    def adHoc_pair_coeff(self, parser):
        line = parser.fIn.readline()
        self.pair_coeff.append(line)
        return None

    def adHoc_run(self, parser):
        line = parser.fIn.readline()
        self.run.append(line)
        return None

    def adHoc_timestep(self, parser):
        line = parser.fIn.readline()
        self.timestep.append(line)
        return None

    def adHoc_thermo(self, parser):
        line = parser.fIn.readline()
        self.thermo.append(line)
        return None

    def adHoc_thermo_style(self, parser):
        line = parser.fIn.readline()
        self.thermo_style.append(line)
        return None

    def adHoc_fix(self, parser):
        line = parser.fIn.readline()
        self.fix.append(line)
        return None

    def adHoc_run_style(self, parser):
        line = parser.fIn.readline()
        self.run_style.append(line)
        return None

    def adHoc_run(self, parser):
        line = parser.fIn.readline()
        self.run.append(line)
        return None

    def adHoc_timestep(self, parser):
        line = parser.fIn.readline()
        self.timestep.append(line)
        return None

    def adHoc_dump(self, parser):
        line = parser.fIn.readline()
        self.dump.append(line)
        return None


    def readStyles(self):  # HERE WE COLLECT CALCULATIONS STYLES (ATOM, BONDS, ANGLES, DIHEDRALS, ELECTROSTATICS, PAIR STYLES)

        specs_filter = self.specs_filter
        pm_filter = self.pm_filter
        sb_filter = self.sb_filter

        list_of_styles = {}
        styles_dict    = {}
        for line in specs_filter:
            line_split = line.split()
            if len(line_split)==2:   # the first 2 terms are strings

                index1 = str(line_split[0])
                index2 = str(line_split[1])

                # creat a list
                styles = [index1, index2]

                # create a dictionary
                styles_dict = { index1 : index2 }

            elif len(line_split)>2 and line_split[0] != "thermo_style":  # this reads lj/cut, lj/class2, lj/charmm, lj/gromacs, etc. with cutoff/relevant parameters

                index1 = str(line_split[0])
                index2 = str(line_split[1])

                cut = []
                for i in range(2, len(line_split)):
                    if "#" in line_split[i]:
                        break

                    try:
                        index = float(line_split[i])
                        cut.append(index)
                    except ValueError:
                        index = line_split[i]
                        cut.append(index)

                # creat a list
                styles = [index1, index2, cut]

                # create a dictionary
                styles_dict = { index1 : [index2, cut] }
            list_of_styles.update(styles_dict)

        modify_dict =[]
        for line in pm_filter:       # reading pair_modify specs
            line_split = line.split()

            modify_dict = { line_split[0] : [ line_split[1:] ] }

        special_dict = []
        for line in sb_filter:       # reading special_bonds specs
            line_split = line.split()

            special_dict = { line_split[0] : [ line_split[1:] ] }

        if modify_dict:
            list_of_styles.update(modify_dict)

        if special_dict:
            list_of_styles.update(special_dict)

        return list_of_styles

    def readBonds(self, bondFunctional):   # HERE WE COLLECT BONDS COEFFICIENTS


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


        bond_filt = self.bond_coeff

        list_of_bonds={}
        for line in bond_filt:
            line_split = line.split()


            if bondFunctional == "harmonic":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*(old_div(toEnergy,(toDistance)**2))
                index3 = float(line_split[3])*toDistance

                bond = [ index2, index3 ]
                bond_dict = {index1 : bond }
                list_of_bonds.update(bond_dict)


            if bondFunctional == "class2":   # COMPASS
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toDistance
                index3 = float(line_split[3])*(old_div(toEnergy,(toDistance)**2))
                index4 = float(line_split[4])*(old_div(toEnergy,(toDistance)**3))
                index5 = float(line_split[5])*(old_div(toEnergy,(toDistance)**4))

                bond = [ index2, index3, index4, index5 ]
                bond_dict = {index1 : bond}
                list_of_bonds.update(bond_dict)


            if bondFunctional == "nonlinear":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toDistance
                index4 = float(line_split[4])*toDistance

                bond = [ index2, index3, index4 ]
                bond_dict = {index1 : bond}
                list_of_bonds.update(bond_dict)


            if bondFunctional == "morse":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*(old_div(1,toDistance))
                index4 = float(line_split[4])*toDistance

                bond = [ index2, index3, index4 ]
                bond_dict = {index1 : bond}
                list_of_bonds.update(bond_dict)


        return list_of_bonds

    ########################################################################################################################



    def readAngles(self, angleFunctional):

        angle_filt = self.angle_coeff

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
        for line in angle_filt:
            line_split = line.split()


            if angleFunctional == "harmonic":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toRadians

                angle = [ index2, index3 ]
                angle_dict = {index1 : angle }
                list_of_angles.update(angle_dict)


            if angleFunctional == "class2":   # COMPASS
                pass


            if angleFunctional == "charmm":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toRadians
                index4 = float(line_split[4])*(old_div(toEnergy,(toDistance)**2))
                index5 = float(line_split[5])*toDistance

                angle = [ index2, index3, index4, index5 ]
                angle_dict = {index1 : angle }
                list_of_angles.update(angle_dict)


            if angleFunctional == "cosine":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy

                angle = [ index2 ]
                angle_dict = {index1 : angle }
                list_of_angles.update(angle_dict)


            if angleFunctional == "cosine/delta":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toRadians

                angle = [ index2, index3 ]
                angle_dict = {index1 : angle }
                list_of_angles.update(angle_dict)


            if angleFunctional == "cosine/periodic":   # DREIDING
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = int(line_split[3])
                index4 = int(line_split[4])

                angle = [ index2, index3, index4 ]
                angle_dict = {index1 : angle }
                list_of_angles.update(angle_dict)


            if angleFunctional == "cosine/squared":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toRadians

                angle = [ index2, index3 ]
                angle_dict = {index1 : angle }
                list_of_angles.update(angle_dict)


        return list_of_angles

    ########################################################################################################################

    def readDihedrals(self, dihedralFunctional):

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


        dihedral_filt = self.dihedral_coeff
        # dihedral_filt = [x for x in storeInput if x.startswith("dihedral_coeff")]

        list_of_dihedrals={}
        for line in dihedral_filt:
            line_split = line.split()


            if dihedralFunctional == "harmonic":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = int(line_split[3])
                index4 = int(line_split[4])

                dihedral = [ index2, index3, index4 ]
                dihedral_dict = {index1 : dihedral }
                list_of_dihedrals.update(dihedral_dict)


            if dihedralFunctional == "class2":   # COMPASS
                pass


            if dihedralFunctional == "multi/harmonic":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toEnergy
                index4 = float(line_split[4])*toEnergy
                index5 = float(line_split[5])*toEnergy
                index6 = float(line_split[6])*toEnergy

                dihedral = [ index2, index3, index4, index5, index6 ]
                dihedral_dict = {index1 : dihedral }
                list_of_dihedrals.update(dihedral_dict)


            if dihedralFunctional == "charmm":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = int(line_split[3])
                index4 = float(line_split[4])*toRadians
                index5 = float(line_split[5])

                dihedral = [ index2, index3, index4, index5 ]
                dihedral_dict = {index1 : dihedral }
                list_of_dihedrals.update(dihedral_dict)


            if dihedralFunctional == "opls":   # OPLS aa
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toEnergy
                index4 = float(line_split[4])*toEnergy
                index5 = float(line_split[5])*toEnergy

                dihedral = [ index2, index3, index4, index5 ]
                dihedral_dict = {index1 : dihedral }
                list_of_dihedrals.update(dihedral_dict)


            if dihedralFunctional == "helix":
                index1 = int(line_split[1])
                index2 = float(line_split[2])*toEnergy
                index3 = float(line_split[3])*toEnergy
                index4 = float(line_split[4])*toEnergy

                dihedral = [ index2, index3, index4 ]
                dihedral_dict = {index1 : dihedral }
                list_of_dihedrals.update(dihedral_dict)


        return list_of_dihedrals


    ################################################################################################################################

    def readPairCoeff(self, updateAtomTypes, pairFunctional):  # HERE WE COLLECT PAIR COEFFICIENTS (LJ)

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

        # this currently gather pair coefficients for pairFunctional reported in the list of strings below

        supportedLJFunct      = ['lj/cut', 'lj/cut/coul/cut', 'lj/cut/coul/long', 'lj/cut/coul/debye',
                                 'lj/cut/coul/dsf', 'lj/cut/coul/msm', 'lj/cut/tip4p/cut', 'lj/cut/tip4p/long']

        supportedCHARMMFunct  = ['lj/charmm/coul/charmm', 'lj/charmm/coul/charmm/implicit', 'lj/charmm/coul/long',
                                 'lj/charmm/coul/msm']

        supportedGROMACSFunct = ['lj/gromacs', 'lj/gromacs/coul/gromacs']

        lj_filt = self.pair_coeff
        # lj_filt = [x for x in storeInput if x.startswith("pair_coeff")]

        list_of_ljs = {}
        ljs_dict     = {}
        at_types    = []
        index       = 0
        for line in lj_filt:
            line_split = line.split()


            if pairFunctional in supportedLJFunct:
                index += 1
                atom1 = int(line_split[1])   # pair atom type 1
                atom2 = int(line_split[2])   # pair atom type 2
                eps   = float(line_split[3])*toEnergy     # epsilon
                sigma = float(line_split[4])*toDistance   # sigma

                coeff = [eps, sigma]

                for i in range(5, len(line_split)):  # if another float is present, it is the pair style cutoff(s) for this pair interaction
                    if "#" in line_split[i]:
                        break

                    try:
                        rad = float(line_split[i])*toDistance
                        coeff.append(rad)
                    except ValueError:
                        pass

                        # creat a list
                lj_coeff = [atom1, atom2, coeff]
                at_types.append(lj_coeff)

                # create dictionaries
                lj_pair = { index : [atom1, atom2] }
                ljs_dict.update(lj_pair)

                lj_param = {index : coeff}
                list_of_ljs.update(lj_param)

            else:
                pass


            if pairFunctional in supportedCHARMMFunct:
                index += 1
                atom1 = int(line_split[1])   # pair atom type 1
                atom2 = int(line_split[2])   # pair atom type 2
                eps   = float(line_split[3])*toEnergy     # epsilon
                sigma = float(line_split[4])*toDistance   # sigma
                eps14   = float(line_split[5])*toEnergy     # epsilon 1-4
                sigma14 = float(line_split[6])*toDistance   # sigma   1-4

                coeff = [eps, sigma, eps14, sigma14]

                # creat a list
                lj_coeff = [atom1, atom2, coeff]
                at_types.append(lj_coeff)

                # create dictionaries
                lj_pair = { index : [atom1, atom2] }
                ljs_dict.update(lj_pair)

                lj_param = {index : coeff}
                list_of_ljs.update(lj_param)

            else:
                pass


            if pairFunctional in supportedGROMACSFunct:
                index += 1
                atom1 = int(line_split[1])   # pair atom type 1
                atom2 = int(line_split[2])   # pair atom type 2
                eps   = float(line_split[3])*toEnergy     # epsilon
                sigma = float(line_split[4])*toDistance   # sigma
                inner = float(line_split[5])*toDistance   # inner sigma
                outer = float(line_split[6])*toDistance   # outer sigma

                coeff = [eps, sigma, inner, outer]

                # creat a list
                lj_coeff = [atom1, atom2, coeff]
                at_types.append(lj_coeff)

                # create dictionaries
                lj_pair = { index : [atom1, atom2] }
                ljs_dict.update(lj_pair)

                lj_param = {index : coeff}
                list_of_ljs.update(lj_param)

            else:
                pass


                # if pairFunctional not in [supportedLJFunct, supportedGROMACSFunct, supportedCHARMMFunct]:
                #
                #     index += 1
                #
                # # creat a list
                #     lj_coeff = ['non supported pair style']
                #     at_types.append(lj_coeff)
                #
                # # create dictionaries
                #     lj_pair = { index : ['non supported pair style'] }
                #     ljs_dict.update(lj_pair)
                #
                #     lj_param = { index : ['non supported pair style']}
                #     list_of_ljs.update(lj_param)


        if updateAtomTypes and lj_filt:  # here I create pair styles including the new atom types (to account for atoms of the same type, but with different partial charges)

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


    ########################################################################################################################

    def readEnsemble(self):  # HERE I READ THE INTEGRATION TYPE AND POTENTIAL CONSTRAINT ALGORITHM

        # ensemble_filter = [x for x in storeInput if fnmatch.fnmatch(x, 'fix*')]
        ensemble_filter = self.fix

        for line in ensemble_filter:
            line_split = line.split()

            if line_split[3] == "langevin":
                sampling = "langevin_dynamics"

        if line_split[3] == "nve":
            sampling = "molecular_dynamics"
            ensemble = "NVE"

        if line_split[3] == "nvt":
            sampling = "molecular_dynamics"
            ensemble = "NVT"

        if line_split[3] == "npt":
            sampling = "molecular_dynamics"
            ensemble = "NPT"

        if line_split[3] == "nph":
            sampling = "molecular_dynamics"
            ensemble = "NPH"

        return (ensemble, sampling)

    ################################################################################################################################

    def readTPSettings(self):  # HERE THERMOSTAT/BAROSTAT TARGETS AND RELAXATION TIMES ARE READ

        target_t       = 0
        target_p       = 0
        thermo_tau     = 0
        baro_tau       = 0
        langevin_gamma = 0

        # ensemble_filter = [x for x in storeInput if fnmatch.fnmatch(x, 'fix*')]

        # ensemble_filter = self.fix
        ensemble_filter = [fix for fix in self.fix if "$" not in fix]

        for line in ensemble_filter:
            line_split = line.split()

            if line_split[3] == "nvt":
                target_t = float(line_split[6])
                thermo_tau = float(line_split[7])

            if line_split[3] == "npt":
                target_t = float(line_split[6])
                thermo_tau = float(line_split[7])
                target_p = float(line_split[10])
                baro_tau = float(line_split[11])

            if line_split[3] == "nph":
                target_p = float(line_split[6])
                baro_tau = float(line_split[7])

            if line_split[3] == "langevin":
                target_t = float(line_split[5])
                langevin_gamma = float(line_split[6])

        return (target_t, thermo_tau, langevin_gamma, target_p, baro_tau)


    ################################################################################################################################

    def readIntegratorSettings(self):  # HERE I READ INTEGRATOR SETTINGS (TYPE, TIME STEP, NUMBER OF STEPS, ...)

        # int_filter = [x for x in storeInput if fnmatch.fnmatch(x, 'run_style*')]
        int_filter = self.run_style

        if int_filter:
            for line in int_filter:
                line_split = line.split()
                int_type = line_split[1]
        else:
            int_type = "verlet"  # if no run_style command, the integrator is standard Verlet


        # run_filt = [x for x in storeInput if x.startswith("run")]  # OK FOR A SINGLE RUN INPUT SCRIPT
        run_filt = self.run

        steps = 0
        for line in run_filt:
            line_split = line.split()
            steps = int(line_split[1])


        # ts_filter = [x for x in storeInput if fnmatch.fnmatch(x, 'timestep*')]
        # ts_filter = self.timestep
        ts_filter = [timestep for timestep in self.timestep if "$" not in timestep]


        for line in ts_filter:
            line_split = line.split()

            tstep = float(line_split[1])

        return (int_type, tstep, steps)

    ################################################################################################################################
