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

from LAMMPS.BaseClasses import Neighbor

# Default
#
# 0.3 bin for units = lj, skin = 0.3 sigma
# 2.0 bin for units = real or metal, skin = 2.0 Angstroms
# 0.001 bin for units = si, skin = 0.001 meters = 1.0 mm
# 0.1 bin for units = cgs, skin = 0.1 cm = 1.0 mm


class Bin(Neighbor):
    name = 'bin'
    pass


class Nsq(Neighbor):
    name = 'nsq'
    pass


class Multi(Neighbor):
    name = 'multi'
    pass

