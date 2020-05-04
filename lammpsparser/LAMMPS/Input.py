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

from abc import ABCMeta
from LAMMPS import BaseClasses



class AtomVec(metaclass=ABCMeta):
    pass




class Atom(object):
    def __init__(self, atom_style, ):
        self.atom_style = atom_style




class AtomStyle(metaclass=ABCMeta):
    pass


class AtomStyleAngle(AtomStyle):
    pass


class Commands(BaseClasses.Command):

    def __init__(self):
        self.CMD_DICT={

            # input script commands
            'clear': self.clear,
            'echo': self.echo,
            'ifthenelse': self.ifthenelse,
            'include': self.include,
            'jump': self.jump,
            'label': self.label,
            'log': self.log,
            'next_command': self.next_command,
            'partition': self.partition,
            'print': self.print_cmd,
            'python': self.python,
            'quit': self.quit,
            'shell': self.shell,
            'variable_command': self.variable_command,

            # LAMMPS commands
            'angle_coeff': self.angle_coeff,
            'angle_style': self.angle_style,
            'atom_modify': self.atom_modify,
            'atom_style': self.atom_style,
            'bond_coeff': self.bond_coeff,
            'bond_style': self.bond_style,
            'bond_write': self.bond_write,
            'boundary': self.boundary,
            'box': self.box,
            'comm_modify': self.comm_modify,
            'comm_style': self.comm_style,
            'compute': self.compute,
            'compute_modify': self.compute_modify,
            'dielectric': self.dielectric,
            'dihedral_coeff': self.dihedral_coeff,
            'dihedral_style': self.dihedral_style,
            'dimension': self.dimension,
            'dump': self.dump,
            'dump_modify': self.dump_modify,
            'fix': self.fix,
            'fix_modify': self.fix_modify,
            'group_command': self.group_command,
            'improper_coeff': self.improper_coeff,
            'improper_style': self.improper_style,
            'kspace_modify': self.kspace_modify,
            'kspace_style': self.kspace_style,
            'lattice': self.lattice,
            'mass': self.mass,
            'min_modify': self.min_modify,
            'min_style': self.min_style,
            'molecule': self.molecule,
            'neigh_modify': self.neigh_modify,
            'neighbor_command': self.neighbor_command,
            'newton': self.newton,
            'package': self.package,
            'pair_coeff': self.pair_coeff,
            'pair_modify': self.pair_modify,
            'pair_style': self.pair_style,
            'pair_write': self.pair_write,
            'processors': self.processors,
            'region': self.region,
            'reset_timestep': self.reset_timestep,
            'restart': self.restart,
            'run_style': self.run_style,
            'special_bonds': self.special_bonds,
            'suffix': self.suffix,
            'thermo': self.thermo,
            'thermo_modify': self.thermo_modify,
            'thermo_style': self.thermo_style,
            'timestep': self.timestep,
            'timer_command': self.timer_command,
            'uncompute': self.uncompute,
            'undump': self.undump,
            'unfix': self.unfix,
            'units': self.units,
        }


    def parse(self, line):
        pass


    # input script commands
    def clear(self, *args):
        print(args)
        pass
    def echo(self, *args):
        print(args)
        pass

    def ifthenelse(self, *args):
        print(args)
        pass

    def include(self, *args):
        print(args)
        pass

    def jump(self, *args):
        print(args)
        pass

    def label(self, *args):
        print(args)
        pass

    def log(self, *args):
        print(args)
        pass

    def next_command(self, *args):
        print(args)
        pass

    def partition(self, *args):
        print(args)
        pass

    def print_cmd(self, *args):
        print(args)
        pass

    def python(self, *args):
        print(args)
        pass

    def quit(self, *args):
        print(args)
        pass

    def shell(self, *args):
        print(args)
        pass

    def variable_command(self, *args):
        print(args)
        pass


    # LAMMPS commands




    def angle_coeff(self, *args):
        print(args)
        pass

    def angle_style(self, *args):
        print(args)
        pass

    def atom_modify(self, *args):
        print(args)
        pass

    def atom_style(self, *args):
        print(args)
        pass

    def bond_coeff(self, *args):
        print(args)
        pass

    def bond_style(self, *args):
        print(args)
        pass

    def bond_write(self, *args):
        print(args)
        pass

    def boundary(self, *args):
        print(args)
        pass

    def box(self, *args):
        print(args)
        pass

    def comm_modify(self, *args):
        print(args)
        pass

    def comm_style(self, *args):
        print(args)
        pass

    def compute(self, *args):
        print(args)
        pass

    def compute_modify(self, *args):
        print(args)
        pass

    def dielectric(self, *args):
        print(args)
        pass

    def dihedral_coeff(self, *args):
        print(args)
        pass

    def dihedral_style(self, *args):
        print(args)
        pass

    def dimension(self, *args):
        print(args)
        pass

    def dump(self, *args):
        print(args)
        pass

    def dump_modify(self, *args):
        print(args)
        pass

    def fix(self, *args):
        print(args)
        pass

    def fix_modify(self, *args):
        print(args)
        pass

    def group_command(self, *args):
        print(args)
        pass

    def improper_coeff(self, *args):
        print(args)
        pass

    def improper_style(self, *args):
        print(args)
        pass

    def kspace_modify(self, *args):
        print(args)
        pass

    def kspace_style(self, *args):
        print(args)
        pass

    def lattice(self, *args):
        print(args)
        pass

    def mass(self, *args):
        print(args)
        pass

    def min_modify(self, *args):
        print(args)
        pass

    def min_style(self, *args):
        print(args)
        pass

    def molecule(self, *args):
        print(args)
        pass

    def neigh_modify(self, *args):
        print(args)
        pass

    def neighbor_command(self, *args):
        print(args)
        pass

    def newton(self, *args):
        print(args)
        pass

    def package(self, *args):
        print(args)
        pass

    def pair_coeff(self, *args):
        print(args)
        pass

    def pair_modify(self, *args):
        print(args)
        pass

    def pair_style(self, *args):
        print(args)
        pass

    def pair_write(self, *args):
        print(args)
        pass

    def processors(self, *args):
        print(args)
        pass

    def region(self, *args):
        print(args)
        pass

    def reset_timestep(self, *args):
        print(args)
        pass

    def restart(self, *args):
        print(args)
        pass

    def run_style(self, *args):
        print(args)
        pass

    def special_bonds(self, *args):
        print(args)
        pass

    def suffix(self, *args):
        print(args)
        pass

    def thermo(self, *args):
        print(args)
        pass

    def thermo_modify(self, *args):
        print(args)
        pass

    def thermo_style(self, *args):
        print(args)
        pass

    def timestep(self, *args):
        print(args)
        pass

    def timer_command(self, *args):
        print(args)
        pass

    def uncompute(self, *args):
        print(args)
        pass

    def undump(self, *args):
        print(args)
        pass

    def unfix(self, *args):
        print(args)
        pass

    def units(self, *args):
        print(args)
        pass

    def parse(self, line):

        line_list = line.split()

        if len(line_list) > 0:
            cmd = line_list[0]
        else:
            return

        if len(line_list) > 1:
            args = line_list[1:]
        else:
            args = []

        cmd_fn = self.CMD_DICT.get(cmd, None)

        if cmd_fn is None:
            return

        cmd_fn(args)

        # if cmd == "clear":
        #     self.clear(self, *args)
        # elif cmd == "echo":
        #     self.echo(self, *args)
        # elif cmd == "if":
        #     self.ifthenelse(self, *args)
        # elif cmd == "include":
        #     self.include(self, *args)
        # elif cmd == "jump":
        #     self.jump(self, *args)
        # elif cmd == "label":
        #     self.label(self, *args)
        # elif cmd == "log":
        #     self.log(self, *args)
        # elif cmd == "next":
        #     self.next_command(self, *args)
        # elif cmd == "partition":
        #     self.partition(self, *args)
        # elif cmd == "print":
        #     self.print_cmd(self, *args)
        # elif cmd == "python":
        #     self.python(self, *args)
        # elif cmd == "quit":
        #     self.quit(self, *args)
        # elif cmd == "shell":
        #     self.shell(self, *args)
        # elif cmd == "variable":
        #     self.variable_command(self, *args)
        #
        # elif cmd == "angle_coeff":
        #     self.angle_coeff(self, *args)
        # elif cmd == "angle_style":
        #     self.angle_style(self, *args)
        # elif cmd == "atom_modify":
        #     self.atom_modify(self, *args)
        # elif cmd == "atom_style":
        #     self.atom_style(self, *args)
        # elif cmd == "bond_coeff":
        #     self.bond_coeff(self, *args)
        # elif cmd == "bond_style":
        #     self.bond_style(self, *args)
        # elif cmd == "bond_write":
        #     self.bond_write(self, *args)
        # elif cmd == "boundary":
        #     self.boundary(self, *args)
        # elif cmd == "box": self.box(self, *args)
        # elif cmd == "comm_modify":
        #     self.comm_modify(self, *args)
        # elif cmd == "comm_style":
        #     self.comm_style(self, *args)
        # elif cmd == "compute":
        #     self.compute(self, *args)
        # elif cmd == "compute_modify":
        #     self.compute_modify(self, *args)
        # elif cmd == "dielectric":
        #     self.dielectric(self, *args)
        # elif cmd == "dihedral_coeff":
        #     self.dihedral_coeff(self, *args)
        # elif cmd == "dihedral_style":
        #     self.dihedral_style(self, *args)
        # elif cmd == "dimension":
        #     self.dimension(self, *args)
        # elif cmd == "dump":
        #     self.dump(self, *args)
        # elif cmd == "dump_modify":
        #     self.dump_modify(self, *args)
        # elif cmd == "fix":
        #     self.fix(self, *args)
        # elif cmd == "fix_modify":
        #     self.fix_modify(self, *args)
        # elif cmd == "group":
        #     self.group_command(self, *args)
        # elif cmd == "improper_coeff":
        #     self.improper_coeff(self, *args)
        # elif cmd == "improper_style":
        #     self.improper_style(self, *args)
        # elif cmd == "kspace_modify":
        #     self.kspace_modify(self, *args)
        # elif cmd == "kspace_style":
        #     self.kspace_style(self, *args)
        # elif cmd == "lattice":
        #     self.lattice(self, *args)
        # elif cmd == "mass":
        #     self.mass(self, *args)
        # elif cmd == "min_modify":
        #     self.min_modify(self, *args)
        # elif cmd == "min_style":
        #     self.min_style(self, *args)
        # elif cmd == "molecule":
        #     self.molecule(self, *args)
        # elif cmd == "neigh_modify":
        #     self.neigh_modify(self, *args)
        # elif cmd == "neighbor":
        #     self.neighbor_command(self, *args)
        # elif cmd == "newton":
        #     self.newton(self, *args)
        # elif cmd == "package":
        #     self.package(self, *args)
        # elif cmd == "pair_coeff":
        #     self.pair_coeff(self, *args)
        # elif cmd == "pair_modify":
        #     self.pair_modify(self, *args)
        # elif cmd == "pair_style":
        #     self.pair_style(self, *args)
        # elif cmd == "pair_write":
        #     self.pair_write(self, *args)
        # elif cmd == "processors":
        #     self.processors(self, *args)
        # elif cmd == "region":
        #     self.region(self, *args)
        # elif cmd == "reset_timestep":
        #     self.reset_timestep(self, *args)
        # elif cmd == "restart":
        #     self.restart(self, *args)
        # elif cmd == "run_style":
        #     self.run_style(self, *args)
        # elif cmd == "special_bonds":
        #     self.special_bonds(self, *args)
        # elif cmd == "suffix":
        #     self.suffix(self, *args)
        # elif cmd == "thermo":
        #     self.thermo(self, *args)
        # elif cmd == "thermo_modify":
        #     self.thermo_modify(self, *args)
        # elif cmd == "thermo_style":
        #     self.thermo_style(self, *args)
        # elif cmd == "timestep":
        #     self.timestep(self, *args)
        # elif cmd == "timer":
        #     self.timer_command(self, *args)
        # elif cmd == "uncompute":
        #     self.uncompute(self, *args)
        # elif cmd == "undump":
        #     self.undump(self, *args)
        # elif cmd == "unfix":
        #     self.unfix(self, *args)
        # elif cmd == "units":
        #     self.units(self, *args)

        print(cmd)


if __name__ == '__main__':

    input = Commands()
    file_name = '../test/methane/input.methane_nvt_thermo_style_one'


    with open(file_name) as file:
        for line in file.readlines():

            input.parse(line)

            print(line)
