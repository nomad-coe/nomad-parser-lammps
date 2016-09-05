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