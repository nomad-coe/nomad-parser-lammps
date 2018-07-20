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

from LAMMPS.BaseClasses import Lattice


# Default
# lattice none 1.0

# class none(Lattice):
#     name = 'none'
#     pass


class Sc(Lattice):
    name = 'sc'
    pass


class Bcc(Lattice):
    name = 'bcc'
    pass


class Fcc(Lattice):
    name = 'fcc'
    pass


class Hcp(Lattice):
    name = 'hcp'
    pass


class Diamond(Lattice):
    name = 'diamond'
    pass


class Sq(Lattice):
    name = 'sq'
    pass


class sq2(Lattice):
    name = 'sq2'
    pass


class hex(Lattice):
    name = 'hex'
    pass


class custom(Lattice):
    name = 'custom'
    pass

