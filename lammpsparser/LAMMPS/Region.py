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

from LAMMPS.BaseClasses import Region


class Block(Region):
    name = 'block'

    def __init__(self, *args):
        pass


class Cone(Region):
    name = 'cone'

    def __init__(self, *args):
        pass


class Cylinder(Region):
    name = 'cylinder'

    def __init__(self, *args):
        pass


class Plane(Region):
    name = 'plane'

    def __init__(self, *args):
        pass


class Prism(Region):
    name = 'prism'

    def __init__(self, *args):
        pass


class Sphere(Region):
    name = 'sphere'

    def __init__(self, *args):
        pass


class Union(Region):
    name = 'union'

    def __init__(self, *args):
        pass


class Intersect(Region):
    name = 'intersect'

    def __init__(self, *args):
        pass


if __name__ == '__main__':
    pass
