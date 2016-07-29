# import setup_paths
# import os, sys, json
import logging

from nomadcore.baseclasses import ParserInterface, ParserContext
from LammpsLogParser import LammpsMainParser

# from nomadcore.simple_parser import SimpleMatcher, mainFunction
# from LammpsParserCommon import get_metaInfo, get_parseInfo

logger = logging.getLogger(name="nomad.LammpsParser")


class LammpsContext(ParserContext):
    def __init__(self):
        super(LammpsContext, self).__init__()


class LammpsParser(ParserInterface):
    """This class handles the initial setup before any parsing can happen. It
    determines which version of the software was used to generate the output
    and then sets up a correct main parser.

    After the implementation has been setup, you can parse the files with
    parse().
    """
    def __init__(self, main_file, metainfo_to_keep=None, backend=None, default_units=None, metainfo_units=None):
        super(LammpsParser, self).__init__(main_file, metainfo_to_keep, backend, default_units, metainfo_units)

    def setup_version(self):
        """Setups the version by looking at the output file and the version
        specified in it.
        """
        # Search for the version specification and initialize a correct
        # main parser for this version.

        self.main_parser = LammpsMainParser(self.parser_context.main_file, self.parser_context)

    def get_metainfo_filename(self):
        return "lammps.nomadmetainfo.json"

    def get_parser_info(self):
        return {'name': 'lammps-parser', 'version': '1.0'}


# class LammpsParserContext(object):
#     def __init__(self):
#         """ Initialise variables used within the current superContext """
#         self.dummy = []
#
#     def initialize_values(self):
#         """ Initializes the values of variables in superContexts that are used to parse different files """
#         self.dummy = None
#
#     def startedParsing(self, fInName, parser):
#         """Function is called when the parsing starts.
#
#         Get compiled parser, filename and metadata.
#
#         Args:
#             fInName: The file name on which the current parser is running.
#             parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
#         """
#         self.parser = parser
#         self.fName = fInName
#         # save metadata
#         self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
#         # allows to reset values if the same superContext is used to parse different files
#         self.initialize_values()
#
# # description of the input
#
# mainFileDescription = SimpleMatcher(
#     name='root',
#     weak=True,
#     startReStr="",
#     subMatchers=[
#         SimpleMatcher(
#             name='newRun',
#             startReStr=r"\s*# SampleParser #\s*",
#             repeats=True,
#             required=True,
#             forwardMatch=True,
#             sections=['section_run'],
#             subMatchers=[
#                 SimpleMatcher(
#                     name='header',
#                     startReStr=r"\s*# SampleParser #\s*")
#             ])
#     ])
#
#
# def get_cachingLevelForMetaName(metaInfoEnv):
#     """Sets the caching level for the metadata.
#
#     Args:
#         metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
#
#     Returns:
#         Dictionary with metaname as key and caching level as value.
#     """
#     # manually adjust caching of metadata
#     pass
#     # cachingLevelForMetaName = { 'x_castep_ir_intensity_store' : CachingLevel.Cache,
#     #                             'x_castep_ir_store' : CachingLevel.Cache,
#     #                             'x_castep_vibrationl_frequencies_store' : CachingLevel.Cache,
#     #                             'x_castep_n_iterations_phonons' : CachingLevel.Cache,
#     #                             'x_castep_mulliken_charge_store' : CachingLevel.Cache,
#     #                             'x_castep_orbital_contributions' : CachingLevel.Cache,
#     #                             # 'x_castep_section_cell_optim': CachingLevel.Cache,
#     #                             'x_castep_section_atom_positions_optim' : CachingLevel.Cache,
#     #                             'x_castep_section_eigenvalues':CachingLevel.Cache,
#     #                             'x_castep_section_k_points':CachingLevel.Cache,
#     #                             'x_castep_section_k_band':CachingLevel.Cache,
#     #                             # 'band_energies' : CachingLevel.Cache,
#     #                             # 'band_k_points' : CachingLevel.Cache,
#     #                             'x_castep_basis_set_planewave_cutoff' : CachingLevel.Cache,
#     #                             # 'eigenvalues_values': CachingLevel.Cache,
#     #                             # 'eigenvalues_kpoints':CachingLevel.Cache,
#     #                             'x_castep_total_energy_corrected_for_finite_basis_store': CachingLevel.Cache,
#     #                             'x_castep_frame_time':CachingLevel.Cache,
#     #                             'x_castep_section_SCF_iteration_frame':CachingLevel.Cache,
#     #                             'x_castep_raman_activity_store': CachingLevel.Cache,
#     #                             'x_castep_SCF_frame_energy_gain':CachingLevel.Cache,
#     #                             'x_castep_frame_energy_free':CachingLevel.Cache,
#     #                             'x_castep_frame_energy_total_T0':CachingLevel.Cache,
#     #                             'x_castep_frame_time_scf_iteration_wall_end':CachingLevel.Cache,
#     #                             'x_castep_SCF_frame_energy':CachingLevel.Cache}
#     #
#     # # Set caching for temparary storage variables
#     # for name in metaInfoEnv.infoKinds:
#     #     if (   name.startswith('x_castep_store_')
#     #            or name.startswith('x_castep_cell_')):
#     #         cachingLevelForMetaName[name] = CachingLevel.Cache
#     # return cachingLevelForMetaName
#
#
#
#
# def main():
#     """Main function.
#
#     Set up everything for the parsing of the CASTEP main file and run the parsing.
#     """
#     # get main file description
#     CastepMainFileSimpleMatcher = build_CastepMainFileSimpleMatcher()
#
#     # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json
#     metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
#     metaInfoEnv = get_metaInfo(metaInfoPath)
#
#     # set parser info
#     parserInfo = {'name':'castep-parser', 'version': '1.0'}
#
#     # get caching level for metadata
#     cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv)
#
#     # start parsing
#     mainFunction(mainFileDescription = CastepMainFileSimpleMatcher,
#                  metaInfoEnv = get_metaInfo(),
#                  parserInfo = get_parseInfo(),
#                  cachingLevelForMetaName = cachingLevelForMetaName,
#                  superContext = CastepParserContext(),
#                  defaultSectionCachingLevel = True)
#
#     mainFunction(mainFileDescription, get_metaInfo(), get_parseInfo())
#
#
# if __name__ == "__main__":
#     main()