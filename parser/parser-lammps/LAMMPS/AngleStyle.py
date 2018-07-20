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

from LAMMPS.BaseClasses import AngleStyle


# class none(AngleStyle):
#     name = 'none'
#     pass


class Hybrid(AngleStyle):
    name = 'hybrid'
    pass


class Charmm(AngleStyle):
    name = 'charmm'
    pass


class Class2(AngleStyle):
    name = 'class2'
    pass


class Cosine(AngleStyle):
    name = 'cosine'
    pass


class CosineSquared(AngleStyle):
    name = 'cosine/squared'
    pass


class Harmonic(AngleStyle):
    name = 'harmonic'
    pass

