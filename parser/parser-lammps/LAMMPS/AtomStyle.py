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

from LAMMPS.BaseClasses import AtomStyle
import warnings

    #
    # if style == 'body':
    #     return body(*args)
    #
    # elif style == 'template':
    #     return template(*args)
    #
    # elif style == 'hybrid':
    #
    #     list = hybrid()
    #     args_iter = args.__iter__()
    #
    #     while True:
    #         try:
    #             style = args_iter.__next__()
    #
    #             if style == 'body':
    #                 bstyle = args_iter.__next__()
    #                 nmin = args_iter.__next__()
    #                 nmax = args_iter.__next__()
    #                 list.append(body(bstyle, nmin, nmax))
    #
    #             elif style == 'template':
    #                 templat_id = args_iter.__next__()
    #                 list.append(template(templat_id))
    #
    #             else:
    #                 list.append(eval(style)())
    #
    #         except StopIteration:
    #             break
    #     return list
    # else:
    #     return eval(style)()



class Base(AtomStyle):
    pass


class Angle(AtomStyle):
    name = 'angle'
    pass

class Atomic(AtomStyle):
    name = 'atomic'

class Angle(AtomStyle):
    name = 'angle'
    pass
class Atomic(AtomStyle):
    name = 'atomic'
    pass

class Body(AtomStyle):
    name = 'body'

    def __init__(self, bstyle, Nmin, Nmax):
        self.bstyle = bstyle
        self.nmin = Nmin
        self.nmax = Nmax

class Bond(AtomStyle):
    name = 'bond'
    pass
class Charge(AtomStyle):
    name = 'charge'
    pass
class Dipole(AtomStyle):
    name = 'dipole'
    pass
class Dpd(AtomStyle):
    name = 'dpd'
    pass
class Electron(AtomStyle):
    name = 'electron'
    pass
class Ellipsoid(AtomStyle):
    name = 'ellipsoid'
    pass
class Full(AtomStyle):
    name = 'full'
    pass
class Line(AtomStyle):
    name = 'line'
    pass
class Meso(AtomStyle):
    name = 'meso'
    pass
class Molecular(AtomStyle):
    name = 'molecular'
    pass
class Peri(AtomStyle):
    name = 'peri'
    pass
class Smd(AtomStyle):
    name = 'smd'
    pass
class Sphere(AtomStyle):
    name = 'sphere'
    pass
class Tri(AtomStyle):
    name = 'tri'
    pass

class Wavepacket(AtomStyle):
    pass

class Template(AtomStyle):
    name = 'template'

    def __init__(self, id):
        self.id = id


class Hybrid(AtomStyle):
    name = 'hybrid'

    def __init__(self, *args):
        pass
    #     super().__init__()
    #     self.list = []
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

_STYLE_DICT = {
    'angle': Angle,
    'atomic': Atomic,
    'body': Body,
    'bond': Bond,
    'charge': Charge,
    'dipole': Dipole,
    'dpd': Dpd,
    'electron': Electron,
    'ellipsoid': Ellipsoid,
    'full': Full,
    'line': Line,
    'meso': Meso,
    'molecular': Molecular,
    'peri': Peri,
    'smd': Smd,
    'sphere': Sphere,
    'tri': Tri,
    'wavepacket': Wavepacket,
    'template': Template,
    'hybrid': Hybrid,
}


def atom_style(style, *args):

    cls = _STYLE_DICT.get(style, None)

    if not cls:
        print("Unsupported atom_style: {}".format(style))
        return None

    return cls(*args)




if __name__ == '__main__':
    a = Bond()
    print(a)

    a = atom_style('atomic')
    print(a)

    a = atom_style('bond')
    print(a)

    atom_style('full')
    atom_style('body', 'nparticle', 2, 10)
    atom_style('hybrid', 'charge', 'bond')
    atom_style('hybrid', 'charge', 'body', 'nparticle', 2, 5)
    atom_style('template', 'myMols')