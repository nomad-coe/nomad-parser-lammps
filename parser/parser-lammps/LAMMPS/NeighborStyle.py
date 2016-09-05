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

