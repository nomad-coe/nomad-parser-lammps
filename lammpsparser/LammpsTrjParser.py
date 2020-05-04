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
import os
import sys
import operator
import numpy as np


from .LammpsCommon import get_metaInfo
from .LammpsCommon import converter as lconverter
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.simple_parser import mainFunction
import nomadcore.ActivateLogging
from contextlib import contextmanager


from re import escape as esc




class sformat(object):
    float = r"[-+]?\d*\.\d+"
    int = r"[-+]?\d+"

nformat = {
    'float': r"[-+]?\d*\.\d+",
    'int': r"[-+]?\d+"
}


# import numpy as np
# import setup_paths
# import logging, os, re, sys
# import nomadcore.ActivateLogging

############################################################
# This is the parser for the input file of LAMMPS.
############################################################

logger = logging.getLogger(name="nomad.LammpsDataParser")


@contextmanager # SECTIONS ARE CLOSED AUTOMATICALLY
def open_section(p, name):
    gid = p.openSection(name)
    yield gid
    p.closeSection(name, gid)



class LammpsTrjParserContext(object):

    """Context for parsing LAMMPS input file.

    Attributes:
        # dos_energies: Stores parsed energies.
        # dos_values: Stores parsed DOS values.

    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of writing and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """

    def __init__(self, converter, dump_style, writeMetaData=True):
        """Args:
            writeMetaData: Determines if metadata is written or stored in class attributes.
        """
        self.parser = None
        self.converter = converter
        self.dump_style = dump_style
        self.writeMetaData = writeMetaData

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts and the compiled parser is obtained.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        # allows to reset values if the same superContext is used to parse different files
        # self.dos_energies = None
        # self.dos_values = None

        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()


    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        # self.masses = {}
        self.pbcBool = []

        self.simulationCell = []
        self.atomPositionScaled = []
        self.atomPositionScaledBool = []
        self.atomPosition = []
        self.atomPositionBool = []
        self.atomPositionWrapped = []
        self.atomPositionWrappedBool = []
        self.imageFlagIndex = []
        self.imageFlagIndexBool = []
        self.atomVelocity = []

        pass





    def onClose_section_run(self, backend, gIndex, section):

        pass

    def onClose_section_system(self, backend, gIndex, section):

        p = backend.superBackend
        o = open_section

        toMass     = self.converter.ratioMass
        toDistance = self.converter.ratioDistance
        toTime     = self.converter.ratioTime
        toEnergy   = self.converter.ratioEnergy
        toVelocity = self.converter.ratioVelocity
        toForce    = self.converter.ratioForce
        toTorque   = self.converter.ratioTorque
        toTemp     = self.converter.ratioTemp
        toPress    = self.converter.ratioPress
        toDVisc    = self.converter.ratioDVisc
        toCharge   = self.converter.ratioCharge
        toDipole   = self.converter.ratioDipole
        toEField   = self.converter.ratioEField
        toDensity  = self.converter.ratioDensity

        with o(p, 'section_system'):


            box_bound = section["x_lammps_trj_box_bound_store"][0]
            boundary = list(map(lambda x: x == "pp", box_bound.split()))

            self.pbcBool.append(boundary)


            # Calculating simulation cell vectors frame by frame
            simulationCelldata = []
            for line in section['x_lammps_trj_box_bounds_store']:
                simulationCelldata.append([float(x) for x in line.split()])


            xx = [ simulationCelldata[0][0] - simulationCelldata[0][1], 0, 0 ]
            yy = [ 0, simulationCelldata[1][0] - simulationCelldata[1][1], 0 ]
            zz = [ 0, 0, simulationCelldata[2][0] - simulationCelldata[2][1] ]
            simulationCellstore = [xx, yy, zz]

            self.simulationCell.append(simulationCellstore)
            temp_simulation_cell = [ [ dim*toDistance for dim in box ] for box in simulationCellstore ]
            p.addArrayValues('simulation_cell', np.array(temp_simulation_cell))

            variables = section['x_lammps_trj_variables_store'][0]
            variables = variables.split()

            frameAtomInfo = []
            for line in section["x_lammps_trj_atoms_store"]:
                frameAtomInfo.append(line.split())

            if 'id' in variables:
                idInd = variables.index('id')
                frameAtomInfo.sort(key=lambda x: int(x[idInd]))


            isAtomPosition = False
            if 'x' in variables and 'y' in variables and 'z' in variables:     # if true, unwrapped coord are dumped
                isAtomPosition = True

            self.atomPositionBool.append(isAtomPosition)



            isAtomPositionScaled = False
            if 'xs' in variables and 'ys' in variables and 'zs' in variables:  # if true, scaled coord are dumped
                isAtomPositionScaled = True

            self.atomPositionScaledBool.append(isAtomPositionScaled)

            isImageFlagIndex = False
            if 'ix' in variables and 'iy' in variables and 'iz' in variables:  # if true, image flag indexes are dumped
                isImageFlagIndex = True

            self.imageFlagIndexBool.append(isImageFlagIndex)
            # self.imageFlagIndex = imageFlagIndex


            # Atom velocities
            isAtomVelocity = False
            if 'vx' in variables and 'vy' in variables and 'vz' in variables:  # if true, scaled coord are dumped
                isAtomVelocity = True

            # Atom forces
            isAtomForces = False
            if 'fx' in variables and 'fy' in variables and 'fz' in variables:  # if true, scaled coord are dumped
                isAtomForces = True

            # Atom position (unwrapped)
            if isAtomPosition == True:

                xInd = variables.index('x')
                yInd = variables.index('y')
                zInd = variables.index('z')

                store_single =[]
                frame = section["x_lammps_trj_atoms_store"]
                for line in frame:
                    line = line.split()
                    if isImageFlagIndex:
                        ixInd = variables.index('ix')
                        iyInd = variables.index('iy')
                        izInd = variables.index('iz')

                        store = [float(line[xInd]) + int(line[ixInd]) * np.linalg.norm(simulationCellstore[0]),
                                 float(line[yInd]) + int(line[iyInd]) * np.linalg.norm(simulationCellstore[1]),
                                 float(line[zInd]) + int(line[izInd]) * np.linalg.norm(simulationCellstore[2])]
                    else:
                        store = [float(line[xInd]), float(line[yInd]), float(line[zInd])]

                    store_single.append(store)

                self.atomPosition.append(store_single)

                temp_atom_positions = [ [ crd*toDistance for crd in atom ] for atom in store_single ]
                p.addArrayValues('atom_positions', np.asarray(temp_atom_positions))
                # p.addArrayValues('atom_labels', np.asarray(atomAtLabel))

            # else:
            #     atomPosition = []
            #     atomPositionBool = False


            # Atom position (scaled) [0, 1]
            if isAtomPositionScaled == True:

                xsInd = variables.index('xs')
                ysInd = variables.index('ys')
                zsInd = variables.index('zs')

                atomPositionScaled = []

                store_single_scaled =[]

                frame = section["x_lammps_trj_atoms_store"]
                for line in frame:
                    line = line.split()
                    store = [float(line[xsInd]), float(line[ysInd]), float(line[zsInd])]
                    store_single_scaled.append(store)

                self.atomPositionScaled.append(store_single_scaled)



            # Atom position (converted from scaled positions)

                store = []
                store_single = []
                for atom in store_single_scaled:

                    x = atom[0] * np.linalg.norm(simulationCellstore[0]) + np.min(simulationCellstore[0])
                    y = atom[1] * np.linalg.norm(simulationCellstore[1]) + np.min(simulationCellstore[1])
                    z = atom[2] * np.linalg.norm(simulationCellstore[2]) + np.min(simulationCellstore[2])
                    store = [x, y, z]
                    store_single.append(store)

                self.atomPosition.append(store_single)

                temp_atom_positions = [ [ crd*toDistance for crd in atom ] for atom in store_single ]
                p.addArrayValues('atom_positions', np.asarray(temp_atom_positions))

            # else:
            #     atomPositionScaled = []
            #     atomPositionScaledBool = False
            #     atomPosition = []
            #     atomPositionBool = False


            # Atom velocities
            if isAtomVelocity == True:

                vxInd = variables.index('vx')
                vyInd = variables.index('vy')
                vzInd = variables.index('vz')

                # atomVelocity = []

                store_single_velocity =[]

                frame = section["x_lammps_trj_atoms_store"]
                for line in frame:
                    line = line.split()
                    store = [float(line[vxInd]), float(line[vyInd]), float(line[vzInd])]
                    store_single_velocity.append(store)

                self.atomVelocity.append(store_single_velocity)
                temp_atom_velocities = [ [ vi*toVelocity for vi in atom ] for atom in store_single_velocity ]
                p.addArrayValues('atom_velocities', np.asarray(temp_atom_velocities))

            # Atom Forces
            if isAtomForces == True:

                fxInd = variables.index('fx')
                fyInd = variables.index('fy')
                fzInd = variables.index('fz')

                # atomVelocity = []

                store_single_forces =[]

                frame = section["x_lammps_trj_atoms_store"]
                for line in frame:
                    line = line.split()
                    store = [float(line[fxInd]), float(line[fyInd]), float(line[fzInd])]
                    store_single_forces.append(store)

                self.atomVelocity.append(store_single_forces)




                # atomForce = []
                # for frame in frameAtomInfo:
                #     store_1 = []
                #     for line in frame:
                #
                #         if fxInd and fyInd and fzInd:
                #             store = [ float(line[fxInd]), float(line[fyInd]), float(line[fzInd]) ]
                #             store_1.append(store)
                #
                #     atomForce.append(store_1)

            # self.atomPositionWrapped = atomPositionWrapped
            # self.atomPositionWrappedBool  = atomPositionWrappedBool



            # return (simulationCell, atomPositionScaled, atomPositionScaledBool, atomPosition, atomPositionBool,
            #             atomPositionWrapped, atomPositionWrappedBool , imageFlagIndex, imageFlagIndexBool )




                #
            # xx = [ header[5][0] - header[5][1], 0, 0 ]
            # yy = [ 0, header[6][0] - header[6][1], 0 ]
            # zz = [ 0, 0, header[7][0] - header[7][1] ]
            # store = [xx, yy, zz]

        # self.simulationCell.append()

        pass




def build_LammpsTrjFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the TRJ file of Lammps.



    Returns:
       SimpleMatcher that parses TRJ file of Lammps.
    """
    # Define the output parsing tree for this version
    r_float = "[-+]?\d*\.\d+"  # Regex for a floating point value
    r_int = "[-+]?\d+"         # Regex for an integer
    r_file = "[0-9a-zA-Z_.]+"
    r_word = "[\w]+"
    r_words = "[\w\s]+"

    timestep = SM(
        name='lammps-trj-timestep',
        startReStr=r"ITEM:\sTIMESTEP",
        # sections=['section_topology'],
        subMatchers=[
            SM(
                r"(?P<x_lammps_trj_timestep_store>{0})".format(r_int)
            )
        ])

    numberofatoms = SM (
        name = "lampps-trj-number-of-atoms",
        startReStr=r"ITEM:\sNUMBER\sOF\sATOMS",
        subMatchers=[
            SM(
                r"(?P<x_lammps_trj_number_of_atoms_store>{0})".format(r_int)
            )
        ],
    )

    boxbound = SM (
        name = "lampps-trj-box-bound",
        startReStr=r"ITEM:\sBOX\sBOUNDS\s(?P<x_lammps_trj_box_bound_store>{0}\s{0}\s{0})".format(r_word),
        subMatchers=[
            SM(r"(?P<x_lammps_trj_box_bounds_store>[\s\d.eE-]+)".format(r_int, r_float), repeats=True)
        ],
    )


    atoms = SM (
        name = "lampps-trj-box-bound",
        startReStr=r"ITEM:\sATOMS\s(?P<x_lammps_trj_variables_store>{0})+".format(r_words),
        subMatchers=[
            SM(r"(?P<x_lammps_trj_atoms_store>[\s\deE.-]+)".format(r_int, r_float), repeats=True)
        ],
    )

    return SM(
        name = 'root',
        startReStr = r"",
        sections = ['section_run'],
        forwardMatch = True,
        weak = True,
        subMatchers = [
            SM(
                name = 'root-timestep',
                startReStr = r"ITEM:\sTIMESTEP",
                sections = ['section_system'],
                forwardMatch=True,
                repeats=True,
                subMatchers=[
                    timestep,
                    numberofatoms,
                    boxbound,
                    atoms
                ]
            )
        ]
    )




def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
        CachingLvl: Sets the CachingLevel for the sections dos, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
        'section_run': CachingLvl,
        'section_topology': CachingLvl,
    }

    # Set all temporarly metadata to Cache as they need post-processsing.
    for name in metaInfoEnv.infoKinds:
        if name.endswith('_store'):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main(converter = lconverter(), dump_style="atom", CachingLvl = CachingLevel.Forward):
    """Main function.

    Set up everything for the parsing of the FHI-aims DOS file and run the parsing.

    Args:
        CachingLvl: Sets the CachingLevel for the sections dos, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.
    """
    # get dos file description
    LammpsTrjSimpleMatcher = build_LammpsTrjFileSimpleMatcher()

    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/lammps.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)

    # set parser info
    parserInfo = {'name': 'lammps-data-parser', 'version': '1.0'}

    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)

    # start parsing
    mainFunction(
        mainFileDescription=LammpsTrjSimpleMatcher,
        metaInfoEnv=metaInfoEnv,
        parserInfo=parserInfo,
        cachingLevelForMetaName=cachingLevelForMetaName,
        superContext=LammpsTrjParserContext(converter, dump_style)
    )


if __name__ == "__main__":
    main(dump_style="custom", CachingLvl=CachingLevel.Forward)
