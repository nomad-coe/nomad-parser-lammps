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
from collections import OrderedDict
import numpy as np

from LAMMPS.BaseClasses import AtomStyle

from enum import Enum, unique
import unittest


def read_atoms(atom_style, ):
    pass


class Angle(AtomStyle):
    name = 'angle'
    keys = ['atom_id', 'molecule_id', 'atom_type', 'x', 'y', 'z']

class Atomic(AtomStyle):
    name = 'atomic'
    keys = ['atom_id', 'atom_type', 'x', 'y', 'z']

    def __init__(self, n_atoms = 0):
        self.n_atoms = n_atoms
        self.data = dict()



    def add_atom(self, atom_type, x, y, z):
        pass


    def __str__(self):
        row_format ='{:>15} {}'

        ret_str = '\n'.join([
            "Atoms",
            row_format.format('atom_style:', self.name),
            row_format.format('# of atoms:', self.n_atoms)
        ])
        return ret_str

class Body(AtomStyle):
    name = 'body'
    keys = ['atom_id', 'atom_type', 'bodyflag', 'mass', 'x', 'y', 'z']

class Bond(AtomStyle):
    name = 'bond'
    keys = ['atom_id', 'molecule_id', 'atom_type', 'x', 'y', 'z']

class Charge(AtomStyle):
    name = 'charge'
    keys = ['atom_id', 'atom_type', 'q', 'x', 'y', 'z']

class Dipole(AtomStyle):
    name = 'dipole'
    keys = ['atom_id', 'atom_type', 'q', 'x', 'y', 'z', 'mux', 'muy', 'muz']

class Dpd(AtomStyle):
    name = 'dpd'
    keys = ['atom_id', 'atom_type', 'theta', 'x', 'y', 'z']

class Electron(AtomStyle):
    name = 'electron'
    keys = ['atom_id', 'atom_type', 'q', 'spin', 'eradius', 'x', 'y', 'z']

class Ellipsoid(AtomStyle):
    name = 'ellipsoid'
    keys = ['atom_id', 'atom_type', 'ellipsoidflag', 'density', 'x', 'y', 'z']

class Full(AtomStyle):
    name = 'full'
    keys = ['atom_id', 'molecule_id', 'atom_type', 'q', 'x', 'y', 'z']

class Line(AtomStyle):
    name = 'line'
    keys = ['atom_id', 'molecule_id', 'atom_type', 'lineflag', 'density', 'x', 'y', 'z']

class Meso(AtomStyle):
    name = 'meso'
    keys = ['atom_id', 'atom_type', 'rho', 'e', 'cv', 'x', 'y', 'z']

class Molecular(AtomStyle):
    name = 'molecular'
    keys = ['atom_id', 'molecule_id', 'atom_type', 'x', 'y', 'z']

class Peri(AtomStyle):
    name = 'peri'
    keys = ['atom_id', 'atom_type', 'volume', 'density', 'x', 'y', 'z']

class Smd(AtomStyle):
    name = 'smd'
    keys = ['atom_id', 'atom_type', 'molecule', 'volume', 'mass', 'kernel-radius', 'contact-radius', 'x', 'y', 'z']

class Sphere(AtomStyle):
    name = 'sphere'
    keys = ['atom_id', 'atom_type', 'diameter', 'density', 'x', 'y', 'z']

class Template(AtomStyle):
    name = 'template'
    keys = ['atom_id', 'molecule_id', 'template-index', 'template-atom', 'atom_type', 'x', 'y', 'z']

class Tri(AtomStyle):
    name = 'tri'
    keys = ['atom_id', 'molecule_id', 'atom_type', 'triangleflag', 'density', 'x', 'y', 'z']

class Wavepacket(AtomStyle):
    name = 'wavepacket'
    keys = ['atom_id', 'atom_type', 'charge', 'spin', 'eradius', 'etag', 'cs_re', 'cs_im', 'x', 'y', 'z']

class Hybrid(AtomStyle):
    name = 'hybrid'
    keys = ['atom_id', 'atom_type', 'x', 'y', 'z']







def atom_style(style, *args):
    """Wrapper for atom_style command"""

    if style == 'body':
        return body(*args)

    elif style == 'template':
        return template(*args)

    elif style == 'hybrid':

        list = hybrid()
        args_iter = args.__iter__()

        while True:
            try:
                style = args_iter.__next__()

                if style == 'body':
                    bstyle = args_iter.__next__()
                    nmin = args_iter.__next__()
                    nmax = args_iter.__next__()
                    list.append(body(bstyle, nmin, nmax))

                elif style == 'template':
                    templat_id = args_iter.__next__()
                    list.append(template(templat_id))

                else:
                    list.append(eval(style)())

            except StopIteration:
                break
        return list
    else:
        return eval(style)()


    # # return eval(style)(*args)
    # cls = globals()[style]
    # return cls(*args)

    # # create an itarator
    # args_iter = args.__iter__()
    # while True:
    #     try:
    #         style = args_iter.__next__()
    #
    #
    #     except StopIteration:
    #         break
    #


class _atom_style(metaclass=ABCMeta):
    """Abstract class for atom_style"""
    pass

class _body_style(metaclass=ABCMeta): pass
class angle(_atom_style): pass
class atomic(_atom_style): pass

class body(_atom_style):
    def __init__(self, bstyle, Nmin, Nmax):
        self.bstyle = bstyle
        self.nmin = Nmin
        self.nmax = Nmax

class bond(_atom_style): pass
class charge(_atom_style): pass
class dipole(_atom_style): pass
class dpd(_atom_style): pass
class electron(_atom_style): pass
class ellipsoid(_atom_style): pass
class full(_atom_style): pass
class line(_atom_style): pass
class meso(_atom_style): pass
class molecular(_atom_style): pass
class peri(_atom_style): pass
class smd(_atom_style): pass
class sphere(_atom_style): pass
class tri(_atom_style): pass
# class wavepacket(_atom_style): pass
# not defined yet

class template(_atom_style):
    def __init__(self, id):
        self.id = id


class hybrid(_atom_style, list):
    pass
    # def __init__(self, *args):
    #
    #     # self.list = []
    #
    #     # create an itarator
    #     args_iter = args.__iter__()
    #     while True:
    #         try:
    #             style = args_iter.__next__()
    #
    #             if style == 'body':
    #                 bstyle = args_iter.__next__()
    #                 nmin = args_iter.__next__()
    #                 nmax = args_iter.__next__()
    #                 self.list.append(body(bstyle, nmin, nmax))
    #             elif style == 'template':
    #                 templat_id = args_iter.__next__()
    #                 self.list.append(template(templat_id))
    #             else:
    #                 self.list.append(eval(style))
    #
    #         except StopIteration:
    #             break




# class _Atoms(metaclass=ABCMeta):
#     """
#
#     """
#     pass
#
# @unique
# class _atom_styles(Enum):
#     angle = 1
#     atomic = 2
#     body = 3
#     bond = 4
#     charge = 5
#     dipole = 6
#     dpd = 7
#     electron = 8
#     ellipsoid = 9
#     full = 10
#     line = 11
#     meso = 12
#     molecular = 13
#     peri = 14
#     smd = 15
#     sphere = 16
#     tri = 17
#     template = 18
#     wavepacket = 19
#     hybrid = 20
#
#
# class atom_style():
#     # STYLES = ['angle', 'atomic', 'body', 'bond', 'charge', 'dipole', 'dpd', 'electron', 'ellipsoid', 'full', 'line', 'meso', 'molecular', 'peri', 'smd', 'sphere', 'template', 'tri', 'wavepacket']
#
#     def __init__(self, style, *args):
#         if style not in _atom_styles.__members__:
#             raise NameError('Wrong style type!')
#
#         self.style = _atom_styles[style]
#
#         if style is 'body':
#             self.bstyle = args[0]
#             self.bstyle_args = args[1:]
#
#         self.args = args
#         # for arg in args:
#         #     print(arg)
#
#
#
#
# class s(_Atoms):
#     def __init__(self):
#         self.style = ''




if __name__ == '__main__':

    a = np.array([1,2,3])
    print(a)


    a = Atomic(10)
    print(a)


    pass

    # atom_style('atomic')
    # atom_style('bond')
    # atom_style('full')
    # a = atom_style('body', 'nparticle', 2, 10)
    # print(a.bstyle)
    # print(a.nmin)
    # print(a.nmax)
    #
    # atom_style('hybrid', 'charge', 'bond')
    # a = atom_style('hybrid', 'charge', 'body', 'nparticle', 2, 5)
    # print(a.list)
    # atom_style('template', 'myMols')
