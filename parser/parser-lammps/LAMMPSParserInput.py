import fnmatch
import os, sys
from LAMMPSParserUnitConversion import unitConversion


#examplesPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../test/examples/methane"))
########################################################################################################################
# FIRST I OPEN THE LAMMPS INPUT FILE TO READ UNITS STYLE AND THE LIST OF LOGGED THERMO VARIABLES
########################################################################################################################

n = sys.argv[1]
examplesPath = os.path.dirname(os.path.abspath(n))

lines = open(n).readlines()  # read LAMMPS input file

########################################################################################################################
# HERE I FIND USER DEFINED VARIABLES AND SUBSTITUTE THEIR NUMERIC VALUE THROUGHOUT
# example: variable tempe equal 300.0
########################################################################################################################
storeInput = []
for line in lines:
    line = line.strip('\n' + '').split(' ')
    line = filter(None, line)

    # If line is just empty
    if line != []:
        pass
        storeInput.append(line)


var_name  = []
var_value = []

for line in storeInput:
    if "variable" and "equal" in line:
        name = '$'+line[1]
        try:
            numb = float(line[3])
        except ValueError:
            continue

    else:
        continue

    var_name.append(name)
    var_value.append(numb)

# print var_name, var_value

for i in range(0, len(var_name)):
    storeInput = [[w.replace(var_name[i],str(var_value[i])) for w in line] for line in storeInput]

storeInput = map(' '.join, storeInput)
storeInput = [line for line in storeInput if not line.startswith('#')]  # CLEAR COMMENTED LINES




################################################################################################################################

################################
##### LAMMPS INPUT PARSING #####
################################

################################################################################################################################

def readUnits():  # HERE I READ THE UNITS STYLE

    units_filter = filter(lambda x: fnmatch.fnmatch(x, '*units*'), storeInput)

    units = []
    for line in units_filter:
            line_split = line.split()
            unitsKey  = line_split[0]
            unitsType = line_split[1]

    units = [unitsKey, unitsType]

    unitsDict = { unitsKey : unitsType }

    return (unitsDict, unitsType)

unitsDict, unitsType = readUnits()


#########################################################################################################################################################
# LOAD UNIT CONVERSION

toMass,toDistance,toTime,toEnergy,toVelocity,toForce,toTorque,toTemp,toPress,toDynVisc,toCharge,toDipole,toElField,toDensity = unitConversion(unitsType)
toRadians = 0.0174533

#########################################################################################################################################################


def readLogFileName():   ### here I pick the name of the log file storing the logged thermodynamic information

    log_filter = filter(lambda x: x.startswith("log"), storeInput)

    if not log_filter:
        logFileName = 'log.lammps'  # thermo output is put in the log.lammps file
                                    # with input echo if a file name is not given from the input

    for line in log_filter:
        line_split = line.split()

        logFileName = line_split[1]

    return logFileName


################################################################################################################################

def readDumpFileName():   ### NOTA BENE: we might have more that one dump file (for now we consider just one,
                          ### which is the trajectiory file)

    dump_filter = filter(lambda x: x.startswith("dump "), storeInput)

    dumpFileName    = None
    stepsPrintFrame = 0
    trajDumpStyle   = None
    for line in dump_filter:
        line_split = line.split()

        stepsPrintFrame = int(line_split[4])
        dumpFileName    = line_split[5]
        trajDumpStyle   = line_split[3]

    return (dumpFileName, stepsPrintFrame, trajDumpStyle)


################################################################################################################################

def readDataFileName():   ### here I pick the name of the LAMMPS topology data file

    data_filter = filter(lambda x: x.startswith("read_data "), storeInput)

    dataFileName = None
    for line in data_filter:
        line_split = line.split()

        dataFileName = line_split[1]

    return (dataFileName)


################################################################################################################################

def readStyles():  # HERE WE COLLECT CALCULATIONS STYLES (ATOM, BONDS, ANGLES, DIHEDRALS, ELECTROSTATICS, PAIR STYLES)

    specs_filter = filter(lambda x: fnmatch.fnmatch(x, '*_style*'), lines)
    pm_filter    = filter(lambda x: x.startswith("pair_modify"), lines)
    sb_filter    = filter(lambda x: x.startswith("special_bonds"), lines)

    list_of_styles = {}
    styles_dict    = {}
    for line in specs_filter:
        line_split = line.split()
        if len(line_split)==2:   # the first 2 terms are strings

            index1 = str(line_split[0])
            index2 = str(line_split[1])

    # creat a list
            styles = [index1, index2]

    # create a dictionary
            styles_dict = { index1 : index2 }

        elif len(line_split)>2 and line_split[0] != "thermo_style":  # this reads lj/cut, lj/class2, lj/charmm, lj/gromacs, etc. with cutoff

            index1 = str(line_split[0])
            index2 = str(line_split[1])

            cut = []
            for i in range(2, len(line_split)):
                if "#" in line_split[i]:
                        break

                try:
                    index = float(line_split[i])
                    cut.append(index)
                except ValueError:
                    index = line_split[i]
                    cut.append(index)

            # creat a list
            styles = [index1, index2, cut]

        # create a dictionary
            styles_dict = { index1 : [index2, cut] }
        list_of_styles.update(styles_dict)

    for line in pm_filter:       # reading pair_modify specs
        line_split = line.split()

        modify_dict = { line_split[0] : [ line_split[1:] ] }

    for line in sb_filter:       # reading special_bonds specs
        line_split = line.split()

        special_dict = { line_split[0] : [ line_split[1:] ] }

    list_of_styles.update(modify_dict)
    list_of_styles.update(special_dict)

    return list_of_styles


################################################################################################################################

def readEnsemble():  # HERE I READ THE INTEGRATION TYPE AND POTENTIAL CONSTRAINT ALGORITHM

	ensemble_filter = filter(lambda x: fnmatch.fnmatch(x, 'fix*'), storeInput)

	for line in ensemble_filter:
		line_split = line.split()

		if line_split[3] == "langevin":
			sampling = "langevin_dynamics"

        if line_split[3] == "nve":
            sampling = "molecular_dynamics"
            ensemble = "NVE"

        if line_split[3] == "nvt":
            sampling = "molecular_dynamics"
            ensemble = "NVT"

        if line_split[3] == "npt":
            sampling = "molecular_dynamics"
            ensemble = "NPT"

        if line_split[3] == "nph":
            sampling = "molecular_dynamics"
            ensemble = "NPH"

	return (ensemble, sampling)


################################################################################################################################

def readTPSettings():  # HERE THERMOSTAT/BAROSTAT TARGETS AND RELAXATION TIMES ARE READ

    target_t       = 0
    target_p       = 0
    thermo_tau     = 0
    baro_tau       = 0
    langevin_gamma = 0

    ensemble_filter = filter(lambda x: fnmatch.fnmatch(x, 'fix*'), storeInput)

    for line in ensemble_filter:
        line_split = line.split()

        if line_split[3] == "nvt":
            target_t = float(line_split[6])
            thermo_tau = float(line_split[7])

        if line_split[3] == "npt":
            target_t = float(line_split[6])
            thermo_tau = float(line_split[7])
            target_p = float(line_split[10])
            baro_tau = float(line_split[11])

        if line_split[3] == "nph":
            target_p = float(line_split[6])
            baro_tau = float(line_split[7])

        if line_split[3] == "langevin":
            target_t = float(line_split[5])
            langevin_gamma = float(line_split[6])

    return (target_t, thermo_tau, langevin_gamma, target_p, baro_tau)


################################################################################################################################

def readIntegratorSettings():  # HERE I READ INTEGRATOR SETTINGS (TYPE, TIME STEP, NUMBER OF STEPS, ...)

    int_filter = filter(lambda x: fnmatch.fnmatch(x, 'run_style*'), storeInput)

    if int_filter:
        for line in int_filter:
            line_split = line.split()
            int_type = line_split[1]
    else:
        int_type = "verlet"  # if no run_style command, the integrator is standard Verlet


    run_filt = filter(lambda x: x.startswith("run"), storeInput)  # OK FOR A SINGLE RUN INPUT SCRIPT

    for line in run_filt:
        line_split = line.split()
        steps = int(line_split[1])


    ts_filter = filter(lambda x: fnmatch.fnmatch(x, 'timestep*'), storeInput)

    for line in ts_filter:
        line_split = line.split()

        tstep = float(line_split[1])

    return (int_type, tstep, steps)


int_type, tstep, steps = readIntegratorSettings()


################################################################################################################################

def readPairCoeff(updateAtomTypes, pairFunctional):  # HERE WE COLLECT PAIR COEFFICIENTS (LJ)

    # this currently gather pair coefficients for pairFunctional reported in the list of strings below

    supportedLJFunct      = ['lj/cut', 'lj/cut/coul/cut', 'lj/cut/coul/long', 'lj/cut/coul/debye',
                             'lj/cut/coul/dsf', 'lj/cut/coul/msm', 'lj/cut/tip4p/cut', 'lj/cut/tip4p/long']

    supportedCHARMMFunct  = ['lj/charmm/coul/charmm', 'lj/charmm/coul/charmm/implicit', 'lj/charmm/coul/long',
                             'lj/charmm/coul/msm']

    supportedGROMACSFunct = ['lj/gromacs', 'lj/gromacs/coul/gromacs']

    lj_filt = filter(lambda x: x.startswith("pair_coeff"), storeInput)

    list_of_ljs = {}
    ljs_dict     = {}
    at_types    = []
    index       = 0
    for line in lj_filt:
        line_split = line.split()


        if pairFunctional in supportedLJFunct:
            index += 1
            atom1 = int(line_split[1])   # pair atom type 1
            atom2 = int(line_split[2])   # pair atom type 2
            eps   = float(line_split[3])*toEnergy     # epsilon
            sigma = float(line_split[4])*toDistance   # sigma

            coeff = [eps, sigma]

            for i in range(5, len(line_split)):  # if another float is present, it is the pair style cutoff(s) for this pair interaction
                if "#" in line_split[i]:
                            break

                try:
                    rad = float(line_split[i])*toDistance
                    coeff.append(rad)
                except ValueError:
                        pass

        # creat a list
            lj_coeff = [atom1, atom2, coeff]
            at_types.append(lj_coeff)

        # create dictionaries
            lj_pair = { index : [atom1, atom2] }
            ljs_dict.update(lj_pair)

            lj_param = {index : coeff}
            list_of_ljs.update(lj_param)

        else:
            pass


        if pairFunctional in supportedCHARMMFunct:
            index += 1
            atom1 = int(line_split[1])   # pair atom type 1
            atom2 = int(line_split[2])   # pair atom type 2
            eps   = float(line_split[3])*toEnergy     # epsilon
            sigma = float(line_split[4])*toDistance   # sigma
            eps14   = float(line_split[5])*toEnergy     # epsilon 1-4
            sigma14 = float(line_split[6])*toDistance   # sigma   1-4

            coeff = [eps, sigma, eps14, sigma14]

        # creat a list
            lj_coeff = [atom1, atom2, coeff]
            at_types.append(lj_coeff)

        # create dictionaries
            lj_pair = { index : [atom1, atom2] }
            ljs_dict.update(lj_pair)

            lj_param = {index : coeff}
            list_of_ljs.update(lj_param)

        else:
            pass


        if pairFunctional in supportedGROMACSFunct:
            index += 1
            atom1 = int(line_split[1])   # pair atom type 1
            atom2 = int(line_split[2])   # pair atom type 2
            eps   = float(line_split[3])*toEnergy     # epsilon
            sigma = float(line_split[4])*toDistance   # sigma
            inner = float(line_split[5])*toDistance   # inner sigma
            outer = float(line_split[6])*toDistance   # outer sigma

            coeff = [eps, sigma, eps14, sigma14]

        # creat a list
            lj_coeff = [atom1, atom2, coeff]
            at_types.append(lj_coeff)

        # create dictionaries
            lj_pair = { index : [atom1, atom2] }
            ljs_dict.update(lj_pair)

            lj_param = {index : coeff}
            list_of_ljs.update(lj_param)

        else:
            pass


        # if pairFunctional not in [supportedLJFunct, supportedGROMACSFunct, supportedCHARMMFunct]:
        #
        #     index += 1
        #
        # # creat a list
        #     lj_coeff = ['non supported pair style']
        #     at_types.append(lj_coeff)
        #
        # # create dictionaries
        #     lj_pair = { index : ['non supported pair style'] }
        #     ljs_dict.update(lj_pair)
        #
        #     lj_param = { index : ['non supported pair style']}
        #     list_of_ljs.update(lj_param)


    if updateAtomTypes:  # here I create pair styles including the new atom types (to account for atoms of the same type, but with different partial charges)

        for line in updateAtomTypes:
            if line[0] != line[1]:

                list_of_ljs.setdefault(line[1], [])
                try:
                    list_of_ljs[line[1]].append(list_of_ljs[line[0]][0])
                except KeyError:
                    pass

                try:
                    list_of_ljs[line[1]].append(list_of_ljs[line[0]][1])
                except KeyError:
                    pass

                ljs_dict.setdefault(line[1], [])
                ljs_dict[line[1]].append(line[1])
                ljs_dict[line[1]].append(line[1])


    return (list_of_ljs, ljs_dict)


########################################################################################################################

def readBonds(bondFunctional):   # HERE WE COLLECT BONDS COEFFICIENTS

    bond_filt = filter(lambda x: x.startswith("bond_coeff"), storeInput)

    list_of_bonds={}
    for line in bond_filt:
        line_split = line.split()


        if bondFunctional == "harmonic":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*(toEnergy/(toDistance)**2)
            index3 = float(line_split[3])*toDistance

            bond = [ index2, index3 ]
            bond_dict = {index1 : bond }
            list_of_bonds.update(bond_dict)


        if bondFunctional == "class2":   # COMPASS
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toDistance
            index3 = float(line_split[3])*(toEnergy/(toDistance)**2)
            index4 = float(line_split[4])*(toEnergy/(toDistance)**3)
            index5 = float(line_split[5])*(toEnergy/(toDistance)**4)

            bond = [ index2, index3, index4, index5 ]
            bond_dict = {index1 : bond}
            list_of_bonds.update(bond_dict)


        if bondFunctional == "nonlinear":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toDistance
            index4 = float(line_split[4])*toDistance

            bond = [ index2, index3, index4 ]
            bond_dict = {index1 : bond}
            list_of_bonds.update(bond_dict)


        if bondFunctional == "morse":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*(1/toDistance)
            index4 = float(line_split[4])*toDistance

            bond = [ index2, index3, index4 ]
            bond_dict = {index1 : bond}
            list_of_bonds.update(bond_dict)


    return list_of_bonds


########################################################################################################################

def readAngles(angleFunctional):

    angle_filt = filter(lambda x: x.startswith("angle_coeff"), storeInput)

    list_of_angles={}
    for line in angle_filt:
        line_split = line.split()


        if angleFunctional == "harmonic":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toRadians

            angle = [ index2, index3 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "class2":   # COMPASS
            pass


        if angleFunctional == "charmm":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toRadians
            index4 = float(line_split[4])*(toEnergy/(toDistance)**2)
            index5 = float(line_split[5])*toDistance

            angle = [ index2, index3, index4, index5 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy

            angle = [ index2 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine/delta":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toRadians

            angle = [ index2, index3 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine/periodic":   # DREIDING
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = int(line_split[3])
            index4 = int(line_split[4])

            angle = [ index2, index3, index4 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


        if angleFunctional == "cosine/squared":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toRadians

            angle = [ index2, index3 ]
            angle_dict = {index1 : angle }
            list_of_angles.update(angle_dict)


    return list_of_angles


########################################################################################################################

def readDihedrals(dihedralFunctional):

    dihedral_filt = filter(lambda x: x.startswith("dihedral_coeff"), storeInput)

    list_of_dihedrals={}
    for line in dihedral_filt:
        line_split = line.split()


        if dihedralFunctional == "harmonic":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = int(line_split[3])
            index4 = int(line_split[4])

            dihedral = [ index2, index3, index4 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "class2":   # COMPASS
            pass


        if dihedralFunctional == "multi/harmonic":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toEnergy
            index4 = float(line_split[4])*toEnergy
            index5 = float(line_split[5])*toEnergy
            index6 = float(line_split[6])*toEnergy

            dihedral = [ index2, index3, index4, index5, index6 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "charmm":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = int(line_split[3])
            index4 = float(line_split[4])*toRadians
            index5 = float(line_split[5])

            dihedral = [ index2, index3, index4, index5 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "opls":   # OPLS aa
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toEnergy
            index4 = float(line_split[4])*toEnergy
            index5 = float(line_split[5])*toEnergy

            dihedral = [ index2, index3, index4, index5 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


        if dihedralFunctional == "helix":
            index1 = int(line_split[1])
            index2 = float(line_split[2])*toEnergy
            index3 = float(line_split[3])*toEnergy
            index4 = float(line_split[4])*toEnergy

            dihedral = [ index2, index3, index4 ]
            dihedral_dict = {index1 : dihedral }
            list_of_dihedrals.update(dihedral_dict)


    return list_of_dihedrals


################################################################################################################################

def readLoggedThermoOutput():  # HERE I READ THE LIST OF THERMO VARIABLES THAT ARE LOGGED (excludind user-defined ones)

    logvars_filter = filter(lambda x: fnmatch.fnmatch(x, 'thermo_style*'), storeInput)

    var          = []
    thermo_style = str
    for line in logvars_filter:
        line_split = line.split()

        if line_split[1] == "one" or logvars_filter:
            var = ['step', 'temp', 'epair', 'emol', 'etotal', 'press']
            thermo_style = 'one'

        if line_split[1] == "custom":
            var = line_split[2:]
            thermo_style = 'custom'

        if line_split[1] == "multi":
            thermo_style = 'multi'

    var = [i for i in var if not i.startswith('f_') and not i.startswith('v_') and not i.startswith('c_')]

    return (var, thermo_style)


################################################################################################################################

def simulationTime():

    run_filt   = filter(lambda x: x.startswith("run"), storeInput)
    frame_filt = filter(lambda  x: x.startswith("thermo "), storeInput)

    for line in run_filt:
        line_split = line.split()
        integrationSteps = int(line_split[1]) # total number of integration steps (single run)

        if unitsType == "real":
            time_length = integrationSteps * tstep
            #time_length = str(time_length) + " ns"
            #time_length = { "Simulation time" : time_length }

    for line in frame_filt:
        line_split = line.split()
        stepsPrintThermo = int(line_split[1])  # thermo output is printed every stepsPrintThermo steps
        frame_length = stepsPrintThermo * tstep


    return (frame_length, time_length, stepsPrintThermo, integrationSteps)

################################################################################################################################


