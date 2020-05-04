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

from LAMMPS.BaseClasses import Group

# Default
#
# All atoms belong to the “all” group.


class Delete(Group):
    name = 'delete'
    pass


class Region(Group):
    name = 'region'
    pass


class Type(Group):
    name = 'type'
    pass


class Id(Group):
    name = 'id'
    pass


class Molecule(Group):
    name = 'molecule'
    pass


class Variable(Group):
    name = 'variable'
    pass


class Include(Group):
    name = 'include'
    pass


class Subtract(Group):
    name = 'subtract'
    pass


class Union(Group):
    name = 'union'
    pass


class Intersect(Group):
    name = 'intersect'
    pass


class Dynamic(Group):
    name = 'dynamic'
    pass


class Static(Group):
    name = 'static'
    pass

