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

from abc import ABCMeta, abstractproperty


# Class hierarchy in LAMMPS source code:
# LAMMPS
# - Memory
# - Error
# - Universe
# - Input
#    - Variable
#    - Command
#       - Finish
#       - Special
# - Atom
#    - AtomVec
# - Update
#    - Integrate
#    - Min
# - Neighbor
#    - NeighList
#    - NeighRequest
# - Comm
#    - Irregular
# - Domain
#    - Region
#    - Lattice
# - Force
#    - Pair
#    - Bond
#    - Angle
#    - Dihedral
#    - Improper
#    - KSpace
#       - FFT3D
#       - Remap
# - Modify
#    - Fix
#    - Compute
# - Group
# - Output
#    - Thermo
#    - Dump
#    - WriteRestart
# - Timer


# Build a class hierarchy for easier operator handling

class Lammps(metaclass=ABCMeta):
    """Abstract superclass of all Lammps related classes"""
    _row_format = '{:>15} {}'
    pass


class Data(Lammps):
    def __init__(self):
        from LAMMPS.Domain import Domain
        from LAMMPS.Atom import AtomAtomic
        self.domain = Domain()   # simulation box
        self.atom = AtomAtomic()

        self.system = None
        self.communication = None
        self.computes = None
        self.dumps = None
        self.fixes = None
        self.groups = None
        self.variables = None

    # def __str__(self):
        # return self._row_format.format('Units:', self.style)

    def info(self):
        out = (
            "\n"
            "LAMMPS object:\n"
            "Domain:\n==============\n{domain}\n").format(
            domain=self.domain)
        return out


class Lattice(Lammps):
    """Define a lattice for use by other commands. In LAMMPS, a lattice is simply a set of points in space, determined by a unit cell with basis atoms, that is replicated infinitely in all dimensions. The arguments of the lattice command can be used to define a wide variety of crystallographic lattices."""
    pass


class Region(Lammps):
    pass


class Group(Lammps):
    pass


class Neighbor(Lammps):
    """This command sets parameters that affect the building of pairwise neighbor lists. All atom pairs within a neighbor cutoff distance equal to the their force cutoff plus the skin distance are stored in the list. """
    pass


class Units(Lammps):
    @abstractproperty
    def style(self): pass

    @abstractproperty
    def mass(self): pass

    @abstractproperty
    def distance(self): pass

    @abstractproperty
    def time(self): pass

    @abstractproperty
    def energy(self): pass

    @abstractproperty
    def velocity(self): pass

    @abstractproperty
    def force(self): pass

    @abstractproperty
    def torque(self): pass

    @abstractproperty
    def temperature(self): pass

    @abstractproperty
    def pressure(self): pass

    @abstractproperty
    def dynamic_viscosity(self): pass

    @abstractproperty
    def charge(self): pass

    @abstractproperty
    def dipole(self): pass

    @abstractproperty
    def electric_field(self): pass

    @abstractproperty
    def density(self, dim): pass

    def __str__(self):
        return self._row_format.format('Units:', self.style)

    def info(self):
        out = (
            "\n"
            "Unit sytle:            {style} \n"
            "   mass:               {mass:.6e} \n"
            "   distance:           {distance:.6e} \n"
            "   time:               {time:.6e} \n"
            "   energy:             {energy:.6e} \n"
            "   velocity:           {velocity:.6e} \n"
            "   force:              {force:.6e} \n"
            "   torque:             {torque:.6e} \n"
            "   temperature:        {temperature:.6e} \n"
            "   pressure:           {pressure:.6e} \n"
            "   dynamic_viscosity:  {dynamic_viscosity:.6e} \n"
            "   charge:             {charge:.6e} \n"
            "   dipole:             {dipole:.6e} \n"
            "   electric_field:     {electric_field:.6e} \n"
            "   density:            {density:.6e} \n").format(
            style=self.style,
            mass=self.mass,
            distance=self.distance,
            time=self.time,
            energy=self.energy,
            velocity=self.velocity,
            force=self.force,
            torque=self.torque,
            temperature=self.temperature,
            pressure=self.pressure,
            dynamic_viscosity=self.dynamic_viscosity,
            charge=self.charge,
            dipole=self.dipole,
            electric_field=self.electric_field,
            density=self.density)
        return out


###################################################
# Force fields:
###################################################


class AtomStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('Atom style:', self.name)


class PairStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('Pair style:', self.name)


class BondStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('Bond style:', self.name)


class AngleStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('Angle style:', self.name)


class DihedralStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('Dihedral style:', self.name)


class ImproperStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('Improper style:', self.name)


class KSpaceStyle(Lammps):
    @abstractproperty
    def name(self): pass

    def __str__(self):
        return self._row_format.format('KSpace style:', self.name)


class SpetialBonds(Lammps):
    pass


###################################################
# Abstract class hierarchy:
###################################################

class Input(Lammps):
    """The Input class reads an input script, stores variables,
    and invokes stand- alone commands that are child classes
    of the Command class"""
    pass


class Variable(Input):
    """The Command class is a parent class for certain input script
    commands that perform a one-time operation before/after/between
    simulations or which invoke a simulation."""
    pass


class Command(Input):
    """The Command class is a parent class for certain input script
    commands that perform a one-time operation before/after/between
    simulations or which invoke a simulation."""
    pass


class Finish(Command):
    """The Finish class is instantiated to print statistics to the screen
    after a simulation is performed, by commands like run and minimize."""
    pass


class Special(Command):
    """The Special class walks the bond topology of a molecular system to
    find 1st, 2nd, 3rd neighbors of each atom. It is invoked by several commands,
    like read data, read restart, and replicate."""
    pass


class Atom(Lammps):
    """The Atom class stores all per-atom arrays.
    More precisely, they are allo- cated and stored by the AtomVec class,
    and the Atom class simply stores a pointer to them.
    The AtomVec class is a parent class for atom styles,
    defined by the atom style command."""
    pass


class Neighbor(Lammps):
    """The Neighbor class builds and stores neighbor lists. The NeighList class stores
    a single list (for all atoms). The NeighRequest class is called by pair, fix, or
    compute styles when they need a particular kind of neighbor list."""
    pass


class Domain(Lammps):
    """The Domain class stores the simulation box geometry, as well as geometric Regions
    and any user definition of a Lattice. The latter are defined by region and
    lattice commands in an input script."""
    pass


class Modify(Lammps):
    """The Modify class stores lists of Fix and Compute classes, both of which are parent styles."""
    pass


class Group(Lammps):
    """The Group class manipulates groups that atoms are assigned to via the group command.
    It also computes various attributes of groups of atoms."""
    pass


class Output(Lammps):
    """The Output class is used to generate 3 kinds of output from a LAMMPS simulation:
    thermodynamic information printed to the screen and log file, dump file snapshots, and restart files.
    These correspond to the Thermo, Dump, and WriteRestart classes respectively.
    The Dump class is a parent style."""
    pass


if __name__ == '__main__':
    lmp = Data()
    print(lmp)
    print(lmp.info())
