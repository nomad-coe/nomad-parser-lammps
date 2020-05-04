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

import LAMMPS.AtomStyle as Atom
from LAMMPS.AtomStyle import atom_style
from LAMMPS import Units
from LAMMPS.BaseClasses import Lammps, AtomStyle


class LAMMPS(Lammps):

    def __init__(self):
        self.atom_style = Atom.Atomic()
        self.units      = Units.lj()

    @property
    def atom_style(self):
        return self.__atom_style

    @atom_style.setter
    def atom_style(self, atom_style):
        if isinstance(atom_style, AtomStyle):
            self.__atom_style = atom_style
        else:
            print("wrong")


    def __str__(self):

        ret_str = '\n'.join([
            self._row_format.format('Program name:', 'LAMMPS'),
            str(self.atom_style),
            str(self.units)
        ])

        return ret_str


if __name__ == '__main__':
    l = LAMMPS()
    l.atom_style = atom_style('bond')

    print(l)