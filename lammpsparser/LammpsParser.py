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

# import setup_paths
# import os, sys, json
import logging

from nomadcore.baseclasses import ParserInterface, ParserContext
from .LammpsLogParser import LammpsMainParser

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
    def __init__(self, metainfo_to_keep=None, backend=None, default_units=None, metainfo_units=None):
        super(LammpsParser, self).__init__(metainfo_to_keep, backend, default_units, metainfo_units)

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
