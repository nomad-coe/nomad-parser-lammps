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

from LAMMPS.BaseClasses import KSpaceStyle


class none(KSpaceStyle):
    name = 'none'
    pass


class Ewald(KSpaceStyle):
    name = 'ewald'
    pass


class EwaldDisp(KSpaceStyle):
    name = 'ewald/disp'
    pass


class EwaldOmp(KSpaceStyle):
    name = 'ewald/omp'
    pass


class Pppm(KSpaceStyle):
    name = 'pppm'
    pass


class PppmCg(KSpaceStyle):
    name = 'pppm/cg'
    pass


class PppmDisp(KSpaceStyle):
    name = 'pppm/disp'
    pass


class PppmTip4p(KSpaceStyle):
    name = 'pppm/tip4p'
    pass


class PppmStagger(KSpaceStyle):
    name = 'pppm/stagger'
    pass


class PppmDispTip4p(KSpaceStyle):
    name = 'pppm/disp/tip4p'
    pass


class PppmGpu(KSpaceStyle):
    name = 'pppm/gpu'
    pass


class PppmOmp(KSpaceStyle):
    name = 'pppm/omp'
    pass


class PppmCgOmp(KSpaceStyle):
    name = 'pppm/cg/omp'
    pass


class PppmTip4pOmp(KSpaceStyle):
    name = 'pppm/tip4p/omp'
    pass


class Msm(KSpaceStyle):
    name = 'msm'
    pass


class MsmCg(KSpaceStyle):
    name = 'msm/cg'
    pass


class MsmOmp(KSpaceStyle):
    name = 'msm/omp'
    pass


class MsmCgOmp(KSpaceStyle):
    name = 'msm/cg/omp'
    pass

