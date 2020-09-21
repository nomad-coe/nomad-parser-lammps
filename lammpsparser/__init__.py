from .metainfo import m_env

from nomad.parsing.parser import FairdiParser
from lammpsparser.lammps_parser import LammpsOutput


class LammpsParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/lammps', code_name='LAMMPS', code_homepage='https://lammps.sandia.gov/',
            domain='dft', mainfile_contents_re=r'^LAMMPS')

    def parse(self, filepath, archive, logger=None):
        self._metainfo_env = m_env

        parser = LammpsOutput(filepath, archive, logger)

        parser.parse()
