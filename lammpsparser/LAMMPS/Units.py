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

import logging
from LAMMPS.BaseClasses import Units

logger = logging.getLogger(__name__)


# """
# The units command also sets the timestep size and neighbor skin distance to default values for each style:
# For style lj these are dt = 0.005 tau and skin = 0.3 sigma.
# For style real these are dt = 1.0 fmsec and skin = 2.0 Angstroms.
# For style metal these are dt = 0.001 psec and skin = 2.0 Angstroms.
# For style si these are dt = 1.0e-8 sec and skin = 0.001 meters.
# For style cgs these are dt = 1.0e-8 sec and skin = 0.1 cm.
# For style electron these are dt = 0.001 fmsec and skin = 2.0 Bohr.
# For style micro these are dt = 2.0 microsec and skin = 0.1 micrometers.
# For style nano these are dt = 0.00045 nanosec and skin = 0.1 nanometers.
# """

def units(style, parent=None):
    """Wrapper for unit style command"""
    if not parent:
        print("hehe")


class _units(object):
    # According Lammps soucre code the physical constants from:
    # http://physics.nist.gov/cuu/Constants/Table/allascii.txt
    # using thermochemical calorie = 4.184 J

    # constants
    mole        = 6.022140857e23    # Avogadro constant [mol^-1]
    Kb          = 1.38064852e-23    # Boltzmann constant [J K^-1]
    perm0       = 8.854187817e-12   # vacuum permittivity / electric constant [F m^-1]

    # mass [kilograms]
    kilogram    = 1.
    gram        = 1e-3
    picogram    = 1e-15
    attogram    = 1e-21
    atomicmass  = 1.660539040e-27   # atomic mass units [kg]

    # distance [meters]
    meter       = 1.
    cm          = 1e-2
    micrometer  = 1e-6
    nanometer   = 1e-9
    Angstrom    = 1.00001495e-10
    Bohr        = 0.52917721067e-10 # Bohr radius [m]

    # time [seconds]
    second      = 1.
    microsecond = 1e-6
    nanosecond  = 1e-9
    picosecond  = 1e-12
    femtosecond = 1e-15
    atomictime  = 2.418884326509e-17  # 1.03275e-15       # atomic time units

    # energy [Joules]
    Joule       = 1.
    Kcal        = 4.184             # Thermochemical calorie
    ergs        = 1e-7
    eV          = 1.6021766208e-19  # electron volt [J]
    Hartree     = 4.359744650e-18   # Hartree energy [J]


    # velocity [meters/second]


    # force [Newtons]
    Newton      = 1.
    dyne        = 1e-5

    # torque [Newton-meters]

    # temperature [Kelvin]
    Kelvin      = 1.

    # pressure [Pascals]
    Pascal      = 1.
    atmosphere  = 101325.
    bar         = 1e5

    # dynamic_viscosity [Pascal*second]
    Poise       = 0.1

    # charge [Coulombs (1.6021765e-19 is a proton)]
    Coulomb     = 1.
    picocoulomb = 1e-6
    charge_e    = 1.6021766208e-19
    statcoulomb = 1/2997924580      # statcoulomb (statC) or franklin (Fr) or electrostatic unit of charge (esu)

    # dipole [Coulombs*meters]
    # electric_field [volts/meter]
    volt        = 1.
    statvolt    = 299.792458
    Debye       = 1/299792458 * 1e-21

    # density [kilograms/meter^dim]





class lj(Units):
    """For style lj, all quantities are unitless. Without loss of generality, LAMMPS sets the fundamental quantities mass, sigma, epsilon, and the Boltzmann constant = 1. The masses, distances, energies you specify are multiples of these fundamental values. The formulas relating the reduced or unitless quantity (with an asterisk) to the same quantity with units is also given. Thus you can use the mass & sigma & epsilon values for a specific material and convert the results from a unitless LJ simulation into physical quantities.
        - mass = mass or m
        - distance = sigma, where x* = x / sigma
        - time = tau, where t* = t (epsilon / m / sigma^2)^1/2
        - energy = epsilon, where E* = E / epsilon
        - velocity = sigma/tau, where v* = v tau / sigma
        - force = epsilon/sigma, where f* = f sigma / epsilon
        - torque = epsilon, where t* = t / epsilon
        - temperature = reduced LJ temperature, where T* = T Kb / epsilon
        - pressure = reduced LJ pressure, where P* = P sigma^3 / epsilon
        - dynamic viscosity = reduced LJ viscosity, where eta* = eta sigma^3 / epsilon / tau
        - charge = reduced LJ charge, where q* = q / (4 pi perm0 sigma epsilon)^1/2
        - dipole = reduced LJ dipole, moment where *mu = mu / (4 pi perm0 sigma^3 epsilon)^1/2
        - electric field = force/charge, where E* = E (4 pi perm0 sigma epsilon)^1/2 sigma / epsilon
        - density = mass/volume, where rho* = rho sigma^dim"""

    style = "lj"

    def __init__(self, mass=1., sigma=1., epsilon=1.):
        self._mass = mass
        self._sigma = sigma
        self._epsilon = epsilon

    @property
    def _tau(self):
        from math import sqrt
        return self._sigma * sqrt(self.mass/self._epsilon)

    @property
    def mass(self):
        return self._mass

    @property
    def distance(self):
        return self._sigma

    @property
    def time(self):
        return self._tau

    @property
    def energy(self):
        return self._epsilon

    @property
    def velocity(self):
        return self._sigma/self._tau

    @property
    def force(self):
        return self._epsilon/self._sigma

    @property
    def torque(self):
        return self._epsilon

    @property
    def temperature(self):
        return self._epsilon/_units.Kb

    @property
    def pressure(self):
        return self._epsilon/self._sigma**3

    @property
    def dynamic_viscosity(self):
        return self._tau * self._epsilon / self._sigma**3

    @property
    def charge(self):
        from math import pi, sqrt
        return sqrt(4 * pi * _units.perm0 * self._sigma * self._epsilon)

    @property
    def dipole(self):
        from math import pi, sqrt
        return sqrt(4 * pi * _units.perm0 * self._sigma**3 * self._epsilon)

    @property
    def electric_field(self):
        from math import pi, sqrt
        return self._epsilon / self._sigma / sqrt(4 * pi * _units.perm0 * self._sigma * self._epsilon)

    @property
    def density(self, dim=3):
        return self._mass / self._sigma**dim

    def parse(selfs):
        logger.warn("NOT SUPPORTED: For LJ style all quantities are unitless")

class real(Units):
    """For style real, these are the units:
        - mass = grams/mole
        - distance = Angstroms
        - time = femtoseconds
        - energy = Kcal/mole
        - velocity = Angstroms/femtosecond
        - force = Kcal/mole-Angstrom
        - torque = Kcal/mole
        - temperature = Kelvin
        - pressure = atmospheres
        - dynamic viscosity = Poise
        - charge = multiple of electron charge (1.0 is a proton)
        - dipole = charge*Angstroms
        - electric field = volts/Angstrom
        - density = gram/cm^dim"""

    style = "real"

    def __init__(self):
        pass

    @property
    def mass(self):
        return _units.gram/_units.mole

    @property
    def distance(self):
        return _units.Angstrom

    @property
    def time(self):
        return _units.femtosecond

    @property
    def energy(self):
        return _units.Kcal / _units.mole

    @property
    def velocity(self):
        return _units.Angstrom / _units.femtosecond

    @property
    def force(self):
        return _units.Kcal / _units.mole / _units.Angstrom

    @property
    def torque(self):
        return _units.Kcal / _units.mole

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.atmosphere

    @property
    def dynamic_viscosity(self):
        return _units.Poise

    @property
    def charge(self):
        return _units.charge_e

    @property
    def dipole(self):
        return _units.charge_e * _units.Angstrom

    @property
    def electric_field(self):
        return _units.volt / _units.Angstrom

    @property
    def density(self, dim=3):
        return _units.gram / _units.cm**dim

class metal(Units):
    """For style metal, these are the units:
        - mass = grams/mole
        - distance = Angstroms
        - time = picoseconds
        - energy = eV
        - velocity = Angstroms/picosecond
        - force = eV/Angstrom
        - torque = eV
        - temperature = Kelvin
        - pressure = bars
        - dynamic viscosity = Poise
        - charge = multiple of electron charge (1.0 is a proton)
        - dipole = charge*Angstroms
        - electric field = volts/Angstrom
        - density = gram/cm^dim"""

    style = "metal"

    def __init__(self):
        pass

    @property
    def mass(self):
        return _units.gram/_units.mole

    @property
    def distance(self):
        return _units.Angstrom

    @property
    def time(self):
        return _units.picosecond

    @property
    def energy(self):
        return _units.eV

    @property
    def velocity(self):
        return _units.Angstrom / _units.picosecond

    @property
    def force(self):
        return _units.eV / _units.Angstrom

    @property
    def torque(self):
        return _units.eV

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.bar

    @property
    def dynamic_viscosity(self):
        return _units.Poise

    @property
    def charge(self):
        return _units.charge_e

    @property
    def dipole(self):
        return _units.charge_e * _units.Angstrom

    @property
    def electric_field(self):
        return _units.volt / _units.Angstrom

    @property
    def density(self, dim=3):
        return _units.gram / _units.cm**dim


class si(Units):
    """For style si, these are the units:
        - mass = kilograms
        - distance = meters
        - time = seconds
        - energy = Joules
        - velocity = meters/second
        - force = Newtons
        - torque = Newton-meters
        - temperature = Kelvin
        - pressure = Pascals
        - dynamic viscosity = Pascal*second
        - charge = Coulombs (1.6021765e-19 is a proton)
        - dipole = Coulombs*meters
        - electric field = volts/meter
        - density = kilograms/meter^dim"""

    style = "si"

    def __init__(self):
        pass


    @property
    def mass(self):
        return _units.kilogram

    @property
    def distance(self):
        return _units.meter

    @property
    def time(self):
        return _units.second

    @property
    def energy(self):
        return _units.Joule

    @property
    def velocity(self):
        return _units.meter / _units.second

    @property
    def force(self):
        return _units.Newton

    @property
    def torque(self):
        return _units.Newton * _units.meter

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.Pascal

    @property
    def dynamic_viscosity(self):
        return _units.Pascal * _units.second

    @property
    def charge(self):
        return _units.Coulomb

    @property
    def dipole(self):
        return _units.Coulomb * _units.meter

    @property
    def electric_field(self):
        return _units.volt / _units.meter

    @property
    def density(self, dim=3):
        return _units.kilogram / _units.meter**dim


class cgs(Units):
    """For style cgs, these are the units:
        - mass = grams
        - distance = centimeters
        - time = seconds
        - energy = ergs
        - velocity = centimeters/second
        - force = dynes
        - torque = dyne-centimeters
        - temperature = Kelvin
        - pressure = dyne/cm^2 or barye = 1.0e-6 bars
        - dynamic viscosity = Poise
        - charge = statcoulombs or esu (4.8032044e-10 is a proton)
        - dipole = statcoul-cm = 10^18 debye
        - electric field = statvolt/cm or dyne/esu
        - density = grams/cm^dim"""

    style = "cgs"

    def __init__(self):
        pass

    @property
    def mass(self):
        return _units.gram

    @property
    def distance(self):
        return _units.cm

    @property
    def time(self):
        return _units.second

    @property
    def energy(self):
        return _units.ergs

    @property
    def velocity(self):
        return _units.cm/_units.second

    @property
    def force(self):
        return _units.dyne

    @property
    def torque(self):
        return _units.dyne * _units.cm

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.dyne/_units.cm**2

    @property
    def dynamic_viscosity(self):
        return _units.Poise

    @property
    def charge(self):
        return _units.statcoulomb

    @property
    def dipole(self):
        return _units.statcoulomb * _units.cm

    @property
    def electric_field(self):
        return _units.statvolt/_units.cm

    @property
    def density(self, dim=3):
        return _units.gram/_units.cm**dim

class electron(Units):
    """For style electron, these are the units:
        - mass = atomic mass units
        - distance = Bohr
        - time = femtoseconds
        - energy = Hartrees
        - velocity = Bohr/atomic time units [1.03275e-15 seconds]
        - force = Hartrees/Bohr
        - temperature = Kelvin
        - pressure = Pascals
        - charge = multiple of electron charge (1.0 is a proton)
        - dipole moment = Debye
        - electric field = volts/cm"""

    style = "electron"

    def __init__(self):
        pass

    @property
    def mass(self):
        return _units.atomicmass

    @property
    def distance(self):
        return _units.Bohr

    @property
    def time(self):
        return _units.femtosecond

    @property
    def energy(self):
        return _units.Hartree

    @property
    def velocity(self):
        return _units.Bohr / _units.atomictime

    @property
    def force(self):
        return _units.Hartree / _units.Bohr

    @property
    def torque(self):
        return _units.Hartree

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.Pascal

    @property
    def dynamic_viscosity(self):
        return _units.Poise

    @property
    def charge(self):
        return _units.charge_e

    @property
    def dipole(self):
        return _units.Debye

    @property
    def electric_field(self):
        return _units.volt / _units.cm

    @property
    def density(self, dim=3):
        return _units.gram / _units.cm**dim

class micro(Units):
    """For style micro, these are the units:
        - mass = picograms
        - distance = micrometers
        - time = microseconds
        - energy = picogram-micrometer^2/microsecond^2
        - velocity = micrometers/microsecond
        - force = picogram-micrometer/microsecond^2
        - torque = picogram-micrometer^2/microsecond^2
        - temperature = Kelvin
        - pressure = picogram/(micrometer-microsecond^2)
        - dynamic viscosity = picogram/(micrometer-microsecond)
        - charge = picocoulombs (1.6021765e-7 is a proton)
        - dipole = picocoulomb-micrometer
        - electric field = volt/micrometer
        - density = picograms/micrometer^dim"""

    style = "micro"

    def __init__(self):
        pass

    @property
    def mass(self):
        return _units.picogram

    @property
    def distance(self):
        return _units.micrometer

    @property
    def time(self):
        return _units.microsecond

    @property
    def energy(self):
        return _units.picogram * _units.micrometer**2 / _units.microsecond**2

    @property
    def velocity(self):
        return _units.micrometer / _units.microsecond

    @property
    def force(self):
        return _units.picogram * _units.micrometer / _units.microsecond**2

    @property
    def torque(self):
        return _units.picogram * _units.micrometer**2 / _units.microsecond**2

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.picogram / (_units.micrometer * _units.microsecond**2)

    @property
    def dynamic_viscosity(self):
        return _units.picogram / (_units.micrometer * _units.microsecond)

    @property
    def charge(self):
        return _units.picocoulomb

    @property
    def dipole(self):
        return _units.picocoulomb * _units.micrometer

    @property
    def electric_field(self):
        return _units.volt / _units.micrometer

    @property
    def density(self, dim=3):
        return _units.picogram / _units.micrometer**dim

class nano(Units):
    """For style nano, these are the units:
        - mass = attograms
        - distance = nanometers
        - time = nanoseconds
        - energy = attogram-nanometer^2/nanosecond^2
        - velocity = nanometers/nanosecond
        - force = attogram-nanometer/nanosecond^2
        - torque = attogram-nanometer^2/nanosecond^2
        - temperature = Kelvin
        - pressure = attogram/(nanometer-nanosecond^2)
        - dynamic viscosity = attogram/(nanometer-nanosecond)
        - charge = multiple of electron charge (1.0 is a proton)
        - dipole = charge-nanometer
        - electric field = volt/nanometer
        - density = attograms/nanometer^dim"""

    style = "nano"

    def __init__(self):
        pass

    @property
    def mass(self):
        return _units.attogram

    @property
    def distance(self):
        return _units.nanometer

    @property
    def time(self):
        return _units.nanometer

    @property
    def energy(self):
        return _units.attogram * _units.nanometer**2 / _units.nanosecond**2

    @property
    def velocity(self):
        return _units.nanometer / _units.nanosecond

    @property
    def force(self):
        return _units.attogram * _units.nanometer / _units.nanosecond**2

    @property
    def torque(self):
        return _units.attogram * _units.nanometer**2 / _units.nanosecond**2

    @property
    def temperature(self):
        return _units.Kelvin

    @property
    def pressure(self):
        return _units.attogram / (_units.nanometer * _units.nanosecond**2)

    @property
    def dynamic_viscosity(self):
        return _units.attogram / (_units.nanometer * _units.nanosecond)

    @property
    def charge(self):
        return _units.charge_e

    @property
    def dipole(self):
        return _units.charge_e * _units.nanometer

    @property
    def electric_field(self):
        return _units.volt / _units.nanometer

    @property
    def density(self, dim=3):
        return _units.attogram / _units.nanometer**dim


if __name__ == '__main__':

    # u = lj()
    # u.parse()
    # print(u)
    # print(u.info())

    styles = ['lj', 'real', 'metal', 'si', 'cgs', 'electron', 'micro', 'nano']

    for style in styles:
        u = eval(style)()
        print(u)
        print(u.info())

