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

