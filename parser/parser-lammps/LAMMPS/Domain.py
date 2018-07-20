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

from LAMMPS.BaseClasses import Domain as AbsDomain


class Domain(AbsDomain):

    def __init__(self, dimension=3):
        self.dimension = dimension
        self.periodicity = None
        self.boundary = None

        self.minxlo = 0.0
        self.minxhi = 0.0
        self.minylo = 0.0
        self.minyhi = 0.0
        self.minzlo = 0.0
        self.minzhi = 0.0

if __name__ == '__main__':
    pass
