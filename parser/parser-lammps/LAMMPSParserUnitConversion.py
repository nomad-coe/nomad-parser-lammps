#  List of unit conversion coefficient for all LAMMPS units styles
#
#  Multiply relevant parsed values by "toProperty" to convert to NOMAD units
#
def unitConversion(unitsType):

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

        toMass     = 1.66054e-27
        toDistance = 1e-10
        toTime     = 1e-15
        toEnergy   = 6.94786e-21
        toVelocity = 1e5
        toForce    = 6.94786e-11
        toTorque   = 6.94786e-21
        toTemp     = 1
        toPress    = 1.01325e5
        toDynVisc  = 0
        toCharge   = 1.602176565e-19
        toDipole   = 0
        toElField  = 0
        toDensity  = 0



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

        toMass     = 1.66054e-27
        toDistance = 1e-10
        toTime     = 1e-12
        toEnergy   = 1.60218e-19
        toVelocity = 1e2
        toForce    = 1.60218e-9
        toTorque   = 1.60218e-19
        toTemp     = 1
        toPress    = 1e5
        toDynVisc  = 0
        toCharge   = 1.602176565e-19
        toDipole   = 0
        toElField  = 0
        toDensity  = 0



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

        toMass     = 1
        toDistance = 1
        toTime     = 1
        toEnergy   = 1
        toVelocity = 1
        toForce    = 1
        toTorque   = 1
        toTemp     = 1
        toPress    = 1
        toDynVisc  = 1
        toCharge   = 1
        toDipole   = 1
        toElField  = 1
        toDensity  = 1



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

        toMass     = 1e-3
        toDistance = 1e-2
        toTime     = 1
        toEnergy   = 1e-7
        toVelocity = 1e-2
        toForce    = 1e-5
        toTorque   = 1e-7
        toTemp     = 1
        toPress    = 1e-1
        toDynVisc  = 1
        toCharge   = 3.335640951982e-10
        toDipole   = 0
        toElField  = 0
        toDensity  = 0



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

        toMass     = 1.66054e-27
        toDistance = 5.29177249e-11
        toTime     = 1e-15
        toEnergy   = 4.35974e-18
        toVelocity = 1e-2
        toForce    = 8.2387264e-8
        toTorque   = 4.35974e-18
        toTemp     = 1
        toPress    = 1
        toDynVisc  = 0
        toCharge   = 1.602176565e-19
        toDipole   = 0
        toElField  = 0
        toDensity  = 0



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

        toMass     = 1e-15
        toDistance = 1e-6
        toTime     = 1e-6
        toEnergy   = 1e-15
        toVelocity = 1
        toForce    = 1e-9
        toTorque   = 1e-15
        toTemp     = 1
        toPress    = 1e3
        toDynVisc  = 0
        toCharge   = 1e-12
        toDipole   = 0
        toElField  = 0
        toDensity  = 0


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

        toMass     = 1e-21
        toDistance = 1e-9
        toTime     = 1e-9
        toEnergy   = 1e-21
        toVelocity = 1
        toForce    = 1e-12
        toTorque   = 1e-21
        toTemp     = 1
        toPress    = 1e6
        toDynVisc  = 0
        toCharge   = 1.602176565e-19
        toDipole   = 0
        toElField  = 0
        toDensity  = 0




    return(toMass,toDistance,toTime,toEnergy,toVelocity,toForce,toTorque,toTemp,toPress,toDynVisc,toCharge,toDipole,toElField,toDensity)