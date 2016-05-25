def unitConversion(unitsType):

    if unitsType == 'real':

        # Define unit conversion parameters for "units = real"
        # For style real, these are default LAMMPS units:
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

        toMass     = 1.6605e-27
        toDistance = 1e-10
        toTime     = 1e-15
        toEnergy   = 6.94786e-21
        toVelocity = 1e5
        toForce    = 6.94786e-11
        toTorque   = 6.94786e-21
        toTemp     = 1
        toPress    = 101325
        toDynVisc  = 0
        toCharge   = 1.602176565e-19
        toDipole   = 0
        toElField  = 0
        toDensity  = 0

    return(toMass,toDistance,toTime,toEnergy,toVelocity,toForce,toTorque,toTemp,toPress,toDynVisc,toCharge,toDipole,toElField,toDensity)