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

import os

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl

############################################################
# This file contains functions that are needed
# by more than one parser.
############################################################

# loading metadata from nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))


def get_metaInfo(filePath=metaInfoPath):
    """Loads metadata.

    Args:
        filePath: Location of metadata.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    metaInfoEnv, warnings = loadJsonFile(
        filePath=filePath,
        dependencyLoader=None,
        extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
        uri=None)

    return metaInfoEnv


#  List of unit conversion coefficient for all LAMMPS units styles
#  Multiply relevant parsed values by "toProperty" to convert to NOMAD units
class converter(object):

    def __init__(self, unitsType = 'lj'):
        self.name = unitsType

        if unitsType == 'lj':
            # mass = mass or m
            # distance = sigma, where x* = x / sigma
            # time = tau, where t* = t (epsilon / m / sigma^2)^1/2
            # energy = epsilon, where E* = E / epsilon
            # velocity = sigma/tau, where v* = v tau / sigma
            # force = epsilon/sigma, where f* = f sigma / epsilon
            # torque = epsilon, where t* = t / epsilon
            # temperature = reduced LJ temperature, where T* = T Kb / epsilon
            # pressure = reduced LJ pressure, where P* = P sigma^3 / epsilon
            # dynamic viscosity = reduced LJ viscosity, where eta* = eta sigma^3 / epsilon / tau
            # charge = reduced LJ charge, where q* = q / (4 pi perm0 sigma epsilon)^1/2
            # dipole = reduced LJ dipole, moment where *mu = mu / (4 pi perm0 sigma^3 epsilon)^1/2
            # electric field = force/charge, where E* = E (4 pi perm0 sigma epsilon)^1/2 sigma / epsilon
            # density = mass/volume, where rho* = rho sigma^dim
            ########################################################################################################
            # Conversion not supported!
            ########################################################################################################
            self.ratioMass     = 1
            self.ratioDistance = 1
            self.ratioTime     = 1
            self.ratioEnergy   = 1
            self.ratioVelocity = 1
            self.ratioForce    = 1
            self.ratioTorque   = 1
            self.ratioTemp     = 1
            self.ratioPress    = 1
            self.ratioDVisc    = 1
            self.ratioCharge   = 1
            self.ratioDipole   = 1
            self.ratioEField   = 1
            self.ratioDensity  = 1

        if unitsType == 'real':

            # Define unit conversion parameters for "units = real"
            # For style real, these are the default LAMMPS units:
            #
            #     mass = grams/mole
            #     distance = Angstroms
            #     time = femtoseconds
            #     energy = Kcal/mole
            #     velocity = Angstroms/femtosecond
            #     force = Kcal/mole-Angstrom
            #     torque = Kcal/mole
            #     temperature = Kelvin
            #     pressure = atmospheres
            #     dynamic viscosity = Poise
            #     charge = multiple of electron charge (1.0 is a proton)
            #     dipole = charge*Angstroms
            #     electric field = volts/Angstrom
            #     density = gram/cm^dim

            self.ratioMass     = 1.66054e-27
            self.ratioDistance = 1e-10
            self.ratioTime     = 1e-15
            self.ratioEnergy   = 6.94786e-21
            self.ratioVelocity = 1e5
            self.ratioForce    = 6.94786e-11
            self.ratioTorque   = 6.94786e-21
            self.ratioTemp     = 1
            self.ratioPress    = 1.01325e5
            self.ratioDVisc    = 0
            self.ratioCharge   = 1.602176565e-19
            self.ratioDipole   = 0
            self.ratioEField   = 0
            self.ratioDensity  = 0



        if unitsType == 'metal':

            # Define unit conversion parameters for "units = metal"
            # For style metal, these are the default LAMMPS units:
            #
            #   mass = grams/mole
            #   distance = Angstroms
            #   time = picoseconds
            #   energy = eV
            #   velocity = Angstroms/picosecond
            #   force = eV/Angstrom
            #   torque = eV
            #   temperature = Kelvin
            #   pressure = bars
            #   dynamic viscosity = Poise
            #   charge = multiple of electron charge (1.0 is a proton)
            #   dipole = charge*Angstroms
            #   electric field = volts/Angstrom
            #   density = gram/cm^dim

            self.ratioMass     = 1.66054e-27
            self.ratioDistance = 1e-10
            self.ratioTime     = 1e-12
            self.ratioEnergy   = 1.60218e-19
            self.ratioVelocity = 1e2
            self.ratioForce    = 1.60218e-9
            self.ratioTorque   = 1.60218e-19
            self.ratioTemp     = 1
            self.ratioPress    = 1e5
            self.ratioDVisc    = 0
            self.ratioCharge   = 1.602176565e-19
            self.ratioDipole   = 0
            self.ratioEField   = 0
            self.ratioDensity  = 0



        if unitsType == 'si':

            # Define unit conversion parameters for "units = si"
            # For style si, these are the default LAMMPS units:
            #
            #   mass = kilograms
            #   distance = meters
            #   time = seconds
            #   energy = Joules
            #   velocity = meters/second
            #   force = Newtons
            #   torque = Newton-meters
            #   temperature = Kelvin
            #   pressure = Pascals
            #   dynamic viscosity = Pascal*second
            #   charge = Coulombs (1.6021765e-19 is a proton)
            #   dipole = Coulombs*meters
            #   electric field = volts/meter
            #   density = kilograms/meter^dim

            self.ratioMass     = 1
            self.ratioDistance = 1
            self.ratioTime     = 1
            self.ratioEnergy   = 1
            self.ratioVelocity = 1
            self.ratioForce    = 1
            self.ratioTorque   = 1
            self.ratioTemp     = 1
            self.ratioPress    = 1
            self.ratioDVisc    = 1
            self.ratioCharge   = 1
            self.ratioDipole   = 1
            self.ratioEField   = 1
            self.ratioDensity  = 1


        if unitsType == 'cgs':

            # Define unit conversion parameters for "units = cgs"
            # For style cgs, these are the default LAMMPS units:
            #
            #   mass = grams
            #   distance = centimeters
            #   time = seconds
            #   energy = ergs
            #   velocity = centimeters/second
            #   force = dynes
            #   torque = dyne-centimeters
            #   temperature = Kelvin
            #   pressure = dyne/cm^2 or barye = 1.0e-6 bars
            #   dynamic viscosity = Poise
            #   charge = statcoulombs or esu (4.8032044e-10 is a proton)
            #   dipole = statcoul-cm = 10^18 debye
            #   electric field = statvolt/cm or dyne/esu
            #   density = grams/cm^dim

            self.ratioMass     = 1e-3
            self.ratioDistance = 1e-2
            self.ratioTime     = 1
            self.ratioEnergy   = 1e-7
            self.ratioVelocity = 1e-2
            self.ratioForce    = 1e-5
            self.ratioTorque   = 1e-7
            self.ratioTemp     = 1
            self.ratioPress    = 1e-1
            self.ratioDVisc    = 1
            self.ratioCharge   = 3.335640951982e-10
            self.ratioDipole   = 0
            self.ratioEField   = 0
            self.ratioDensity  = 0

        if unitsType == 'electron':

            # Define unit conversion parameters for "units = electron"
            # For style electron, these are the default LAMMPS units:
            #
            # mass = atomic mass units
            # distance = Bohr
            # time = femtoseconds
            # energy = Hartrees
            # velocity = Bohr/atomic time units [1.03275e-15 seconds]
            # force = Hartrees/Bohr
            # temperature = Kelvin
            # pressure = Pascals
            # charge = multiple of electron charge (1.0 is a proton)
            # dipole moment = Debye
            # electric field = volts/cm

            self.ratioMass     = 1.66054e-27
            self.ratioDistance = 5.29177249e-11
            self.ratioTime     = 1e-15
            self.ratioEnergy   = 4.35974e-18
            self.ratioVelocity = 1e-2
            self.ratioForce    = 8.2387264e-8
            self.ratioTorque   = 4.35974e-18
            self.ratioTemp     = 1
            self.ratioPress    = 1
            self.ratioDVisc    = 0
            self.ratioCharge   = 1.602176565e-19
            self.ratioDipole   = 0
            self.ratioEField   = 0
            self.ratioDensity  = 0

        if unitsType == 'micro':

            # Define unit conversion parameters for "units = micro"
            # For style micro, these are the default LAMMPS units:
            #
            #   mass = picograms
            #   distance = micrometers
            #   time = microseconds
            #   energy = picogram-micrometer^2/microsecond^2
            #   velocity = micrometers/microsecond
            #   force = picogram-micrometer/microsecond^2
            #   torque = picogram-micrometer^2/microsecond^2
            #   temperature = Kelvin
            #   pressure = picogram/(micrometer-microsecond^2)
            #   dynamic viscosity = picogram/(micrometer-microsecond)
            #   charge = picocoulombs (1.6021765e-7 is a proton)
            #   dipole = picocoulomb-micrometer
            #   electric field = volt/micrometer
            #   density = picograms/micrometer^dim

            self.ratioMass     = 1e-15
            self.ratioDistance = 1e-6
            self.ratioTime     = 1e-6
            self.ratioEnergy   = 1e-15
            self.ratconveioVelocity = 1
            self.ratioForce    = 1e-9
            self.ratioTorque   = 1e-15
            self.ratioTemp     = 1
            self.ratioPress    = 1e3
            self.ratioDVisc    = 0
            self.ratioCharge   = 1e-12
            self.ratioDipole   = 0
            self.ratioEField   = 0
            self.ratioDensity  = 0

        if unitsType == 'nano':

            # Define unit conversion parameters for "units = nano"
            # For style nano, these are the default LAMMPS units:
            #
            #   mass = attograms
            #   distance = nanometers
            #   time = nanoseconds
            #   energy = attogram-nanometer^2/nanosecond^2
            #   velocity = nanometers/nanosecond
            #   force = attogram-nanometer/nanosecond^2
            #   torque = attogram-nanometer^2/nanosecond^2
            #   temperature = Kelvin
            #   pressure = attogram/(nanometer-nanosecond^2)
            #   dynamic viscosity = attogram/(nanometer-nanosecond)
            #   charge = multiple of electron charge (1.0 is a proton)
            #   dipole = charge-nanometer
            #   electric field = volt/nanometer
            #   density = attograms/nanometer^dim

            self.ratioMass     = 1e-21
            self.ratioDistance = 1e-9
            self.ratioTime     = 1e-9
            self.ratioEnergy   = 1e-21
            self.ratioVelocity = 1
            self.ratioForce    = 1e-12
            self.ratioTorque   = 1e-21
            self.ratioTemp     = 1
            self.ratioPress    = 1e6
            self.ratioDVisc    = 0
            self.ratioCharge   = 1.602176565e-19
            self.ratioDipole   = 0
            self.ratioEField   = 0
            self.ratioDensity  = 0

    def Mass(self, val):
        return self.ratioMass * val

    def Distance(self, val):
        return self.ratioDistance * val

    def Time(self, val):
        return self.ratioTime * val

    def Energy(self, val):
        return self.ratioEnergy * val

    def Velocity(self, val):
        return self.ratioVelocity * val

    def Force(self, val):
        return self.ratioForce * val

    def Torque(self, val):
        return self.ratioTorque * val

    def Temp(self, val):
        return self.ratioTemp * val

    def Press(self, val):
        return self.ratioPress * val

    def DVisc(self, val):
        return self.ratioDVisc * val

    def Charge(self, val):
        return self.ratioCharge * val

    def Dipole(self, val):
        return self.ratioDipole * val

    def EField(self, val):
        return self.ratioEField * val

    def Density(self, val):
        return self.ratioDensity * val
