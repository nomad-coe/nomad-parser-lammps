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

from LAMMPS.BaseClasses import PairStyle


# class None(PairStyle):
#     name = 'none'
#     pass


class Hybrid(PairStyle):
    name = 'hybrid'
    pass


class HybridOverlay(PairStyle):
    name = 'hybrid/overlay'
    pass


class Zero(PairStyle):
    name = 'zero'
    pass


class Adp(PairStyle):
    name = 'adp'
    pass


class Airebo(PairStyle):
    name = 'airebo'
    pass


class AireboMorse(PairStyle):
    name = 'airebo/morse'
    pass


class Beck(PairStyle):
    name = 'beck'
    pass


class Body(PairStyle):
    name = 'body'
    pass


class Bop(PairStyle):
    name = 'bop'
    pass


class Born(PairStyle):
    name = 'born'
    pass


class BornCoulLong(PairStyle):
    name = 'born/coul/long'
    pass


class BornCoulLongCs(PairStyle):
    name = 'born/coul/long/cs'
    pass


class BornCoulMsm(PairStyle):
    name = 'born/coul/msm'
    pass


class BornCoulWolf(PairStyle):
    name = 'born/coul/wolf'
    pass


class Brownian(PairStyle):
    name = 'brownian'
    pass


class BrownianPoly(PairStyle):
    name = 'brownian/poly'
    pass


class Buck(PairStyle):
    name = 'buck'
    pass


class BuckCoulCut(PairStyle):
    name = 'buck/coul/cut'
    pass


class BuckCoulLong(PairStyle):
    name = 'buck/coul/long'
    pass


class BuckCoulLongCs(PairStyle):
    name = 'buck/coul/long/cs'
    pass


class BuckCoulMsm(PairStyle):
    name = 'buck/coul/msm'
    pass


class BuckLongCoulLong(PairStyle):
    name = 'buck/long/coul/long'
    pass


class Colloid(PairStyle):
    name = 'colloid'
    pass


class Comb(PairStyle):
    name = 'comb'
    pass


class Comb3(PairStyle):
    name = 'comb3'
    pass


class CoulCut(PairStyle):
    name = 'coul/cut'
    pass


class CoulDebye(PairStyle):
    name = 'coul/debye'
    pass


class CoulDsf(PairStyle):
    name = 'coul/dsf'
    pass


class CoulLong(PairStyle):
    name = 'coul/long'
    pass


class CoulLongCs(PairStyle):
    name = 'coul/long/cs'
    pass


class CoulMsm(PairStyle):
    name = 'coul/msm'
    pass


class CoulStreitz(PairStyle):
    name = 'coul/streitz'
    pass


class CoulWolf(PairStyle):
    name = 'coul/wolf'
    pass


class Dpd(PairStyle):
    name = 'dpd'
    pass


class DpdTstat(PairStyle):
    name = 'dpd/tstat'
    pass


class Dsmc(PairStyle):
    name = 'dsmc'
    pass


class Eam(PairStyle):
    name = 'eam'
    pass


class EamAlloy(PairStyle):
    name = 'eam/alloy'
    pass


class EamFs(PairStyle):
    name = 'eam/fs'
    pass


class Eim(PairStyle):
    name = 'eim'
    pass


class Gauss(PairStyle):
    name = 'gauss'
    pass


class Gayberne(PairStyle):
    name = 'gayberne'
    pass


class GranHertzHistory(PairStyle):
    name = 'gran/hertz/history'
    pass


class GranHooke(PairStyle):
    name = 'gran/hooke'
    pass


class GranHookeHistory(PairStyle):
    name = 'gran/hooke/history'
    pass


class HbondDreidingLj(PairStyle):
    name = 'hbond/dreiding/lj'
    pass


class HbondDreidingMorse(PairStyle):
    name = 'hbond/dreiding/morse'
    pass


class Kim(PairStyle):
    name = 'kim'
    pass


class Lcbop(PairStyle):
    name = 'lcbop'
    pass


class LineLj(PairStyle):
    name = 'line/lj'
    pass


class LjCharmmCoulCharmm(PairStyle):
    name = 'lj/charmm/coul/charmm'
    pass


class LjCharmmCoulCharmmImplicit(PairStyle):
    name = 'lj/charmm/coul/charmm/implicit'
    pass


class LjCharmmCoulLong(PairStyle):
    name = 'lj/charmm/coul/long'
    pass


class LjCharmmCoulMsm(PairStyle):
    name = 'lj/charmm/coul/msm'
    pass


class LjClass2(PairStyle):
    name = 'lj/class2'
    pass


class LjClass2CoulCut(PairStyle):
    name = 'lj/class2/coul/cut'
    pass


class LjClass2CoulLong(PairStyle):
    name = 'lj/class2/coul/long'
    pass


class LjCubic(PairStyle):
    name = 'lj/cubic'
    pass


class LjCut(PairStyle):
    name = 'lj/cut'
    pass


class LjCutCoulCut(PairStyle):
    name = 'lj/cut/coul/cut'
    pass


class LjCutCoulDebye(PairStyle):
    name = 'lj/cut/coul/debye'
    pass


class LjCutCoulDsf(PairStyle):
    name = 'lj/cut/coul/dsf'
    pass


class LjCutCoulLong(PairStyle):
    name = 'lj/cut/coul/long'
    pass


class LjCutCoulLongCs(PairStyle):
    name = 'lj/cut/coul/long/cs'
    pass


class LjCutCoulMsm(PairStyle):
    name = 'lj/cut/coul/msm'
    pass


class LjCutDipoleCut(PairStyle):
    name = 'lj/cut/dipole/cut'
    pass


class LjCutDipoleLong(PairStyle):
    name = 'lj/cut/dipole/long'
    pass


class LjCutTip4pCut(PairStyle):
    name = 'lj/cut/tip4p/cut'
    pass


class LjCutTip4pLong(PairStyle):
    name = 'lj/cut/tip4p/long'
    pass


class LjExpand(PairStyle):
    name = 'lj/expand'
    pass


class LjGromacs(PairStyle):
    name = 'lj/gromacs'
    pass


class LjGromacsCoulGromacs(PairStyle):
    name = 'lj/gromacs/coul/gromacs'
    pass


class LjLongCoulLong(PairStyle):
    name = 'lj/long/coul/long'
    pass


class LjLongDipoleLong(PairStyle):
    name = 'lj/long/dipole/long'
    pass


class LjLongTip4pLong(PairStyle):
    name = 'lj/long/tip4p/long'
    pass


class LjSmooth(PairStyle):
    name = 'lj/smooth'
    pass


class LjSmoothLinear(PairStyle):
    name = 'lj/smooth/linear'
    pass


class Lj96Cut(PairStyle):
    name = 'lj96/cut'
    pass


class Lubricate(PairStyle):
    name = 'lubricate'
    pass


class LubricatePoly(PairStyle):
    name = 'lubricate/poly'
    pass


class LubricateU(PairStyle):
    name = 'lubricate/u'
    pass


class LubricateUPoly(PairStyle):
    name = 'lubricate/u/poly'
    pass


class Meam(PairStyle):
    name = 'meam'
    pass


class MieCut(PairStyle):
    name = 'mie/cut'
    pass


class Morse(PairStyle):
    name = 'morse'
    pass


class Nb3bHarmonic(PairStyle):
    name = 'nb3b/harmonic'
    pass


class NmCut(PairStyle):
    name = 'nm/cut'
    pass


class NmCutCoulCut(PairStyle):
    name = 'nm/cut/coul/cut'
    pass


class NmCutCoulLong(PairStyle):
    name = 'nm/cut/coul/long'
    pass


class PeriEps(PairStyle):
    name = 'peri/eps'
    pass


class PeriLps(PairStyle):
    name = 'peri/lps'
    pass


class PeriPmb(PairStyle):
    name = 'peri/pmb'
    pass


class PeriVes(PairStyle):
    name = 'peri/ves'
    pass


class Polymorphic(PairStyle):
    name = 'polymorphic'
    pass


class Reax(PairStyle):
    name = 'reax'
    pass


class Rebo(PairStyle):
    name = 'rebo'
    pass


class Resquared(PairStyle):
    name = 'resquared'
    pass


class Snap(PairStyle):
    name = 'snap'
    pass


class Soft(PairStyle):
    name = 'soft'
    pass


class Sw(PairStyle):
    name = 'sw'
    pass


class Table(PairStyle):
    name = 'table'
    pass


class Tersoff(PairStyle):
    name = 'tersoff'
    pass


class TersoffMod(PairStyle):
    name = 'tersoff/mod'
    pass


class TersoffZbl(PairStyle):
    name = 'tersoff/zbl'
    pass


class Tip4pCut(PairStyle):
    name = 'tip4p/cut'
    pass


class Tip4pLong(PairStyle):
    name = 'tip4p/long'
    pass


class TriLj(PairStyle):
    name = 'tri/lj'
    pass


class Vashishta(PairStyle):
    name = 'vashishta'
    pass


class Yukawa(PairStyle):
    name = 'yukawa'
    pass


class YukawaColloid(PairStyle):
    name = 'yukawa/colloid'
    pass


class Zbl(PairStyle):
    name = 'zbl'
    pass

