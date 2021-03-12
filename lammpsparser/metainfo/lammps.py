#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import public
from nomad.datamodel.metainfo import common

m_package = Package(
    name='lammps_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='lammps.nomadmetainfo.json'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_lammps_atom_positions_image_index = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms', 3],
        unit='dimensionless',
        description='''
        PBC image flag index.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_atom_positions_image_index'))

    x_lammps_atom_positions_scaled = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='dimensionless',
        description='''
        Position of the atoms in a scaled format [0, 1].
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_atom_positions_scaled'))

    x_lammps_atom_positions_wrapped = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='meter',
        description='''
        Position of the atoms wrapped back to the periodic box.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_atom_positions_wrapped'))

    x_lammps_trj_timestep_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_trj_timestep_store'))

    x_lammps_trj_number_of_atoms_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_trj_number_of_atoms_store'))

    x_lammps_trj_box_bound_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_trj_box_bound_store'))

    x_lammps_trj_box_bounds_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_trj_box_bounds_store'))

    x_lammps_trj_variables_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_trj_variables_store'))

    x_lammps_trj_atoms_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_trj_atoms_store'))


class section_sampling_method(public.section_sampling_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_sampling_method'))

    x_lammps_barostat_target_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        MD barostat target pressure.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_barostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_barostat_target_pressure'))

    x_lammps_barostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD barostat relaxation time.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_barostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_barostat_tau'))

    x_lammps_barostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD barostat type, valid values are defined in the barostat_type wiki page.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_barostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_barostat_type'))

    x_lammps_integrator_dt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD integration time step.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_sampling, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_lammps_integrator_dt'))

    x_lammps_integrator_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD integrator type, valid values are defined in the integrator_type wiki page.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_sampling, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_lammps_integrator_type'))

    x_lammps_langevin_gamma = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Langevin thermostat damping factor.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_langevin_gamma'))

    x_lammps_number_of_steps_requested = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of requested MD integration time steps.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_sampling, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_lammps_number_of_steps_requested'))

    x_lammps_thermostat_level = Quantity(
        type=str,
        shape=[],
        description='''
        MD thermostat level (see wiki: single, multiple, regional).
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_thermostat_level'))

    x_lammps_thermostat_target_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kelvin',
        description='''
        MD thermostat target temperature.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_thermostat_target_temperature'))

    x_lammps_thermostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD thermostat relaxation time.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_thermostat_tau'))

    x_lammps_thermostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD thermostat type, valid values are defined in the thermostat_type wiki page.
        ''',
        categories=[public.settings_molecular_dynamics, public.settings_thermostat, public.settings_sampling],
        a_legacy=LegacyDefinition(name='x_lammps_thermostat_type'))


class section_interaction(common.section_interaction):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_interaction'))

    x_lammps_interaction_atom_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type of each interaction atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_interaction_atom_to_atom_type_ref'))

    x_lammps_number_of_defined_pair_interactions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of defined pair interactions (L-J pairs).
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_number_of_defined_pair_interactions'))

    x_lammps_pair_interaction_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['x_lammps_number_of_defined_pair_interactions', 'number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type for pair interactions.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_pair_interaction_atom_type_ref'))

    x_lammps_pair_interaction_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=['x_lammps_number_of_defined_pair_interactions', 2],
        description='''
        Pair interactions parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_pair_interaction_parameters'))


class section_molecule_interaction(common.section_molecule_interaction):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_molecule_interaction'))

    x_lammps_molecule_interaction_atom_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type of each molecule interaction atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_molecule_interaction_atom_to_atom_type_ref'))

    x_lammps_number_of_defined_molecule_pair_interactions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of defined pair interactions within a molecule (L-J pairs).
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_number_of_defined_molecule_pair_interactions'))

    x_lammps_pair_molecule_interaction_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_defined_molecule_pair_interactions', 2],
        description='''
        Molecule pair interactions parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_pair_molecule_interaction_parameters'))

    x_lammps_pair_molecule_interaction_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['x_lammps_number_of_defined_pair_interactions', 'number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type for pair interactions within a molecule.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_pair_molecule_interaction_to_atom_type_ref'))


class x_lammps_section_input_output_files(MSection):
    '''
    Section to store input and output file names
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_lammps_section_input_output_files'))

    x_lammps_inout_file_data = Quantity(
        type=str,
        shape=[],
        description='''
        Lammps input data file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_file_data'))

    x_lammps_inout_file_trajectory = Quantity(
        type=str,
        shape=[],
        description='''
        Lammps input trajectory file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_file_trajectory'))


class x_lammps_section_control_parameters(MSection):
    '''
    Section to store the input and output control parameters
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_lammps_section_control_parameters'))

    x_lammps_inout_control_anglecoeff = Quantity(
        type=str,
        shape=[],
        description='''
        Specify the angle force field coefficients for one or more angle types.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_anglecoeff'))

    x_lammps_inout_control_anglestyle = Quantity(
        type=str,
        shape=[],
        description='''
        Set the formula(s) LAMMPS uses to compute angle interactions between triplets of
        atoms, which remain in force for the duration of the simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_anglestyle'))

    x_lammps_inout_control_atommodify = Quantity(
        type=str,
        shape=[],
        description='''
        Modify certain attributes of atoms defined and stored within LAMMPS, in addition
        to what is specified by the atom_style command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_atommodify'))

    x_lammps_inout_control_atomstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Define what style of atoms to use in a simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_atomstyle'))

    x_lammps_inout_control_balance = Quantity(
        type=str,
        shape=[],
        description='''
        This command adjusts the size and shape of processor sub-domains within the
        simulation box, to attempt to balance the number of atoms or particles and thus
        indirectly the computational cost (load) more evenly across processors.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_balance'))

    x_lammps_inout_control_bondcoeff = Quantity(
        type=str,
        shape=[],
        description='''
        Specify the bond force field coefficients for one or more bond type.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_bondcoeff'))

    x_lammps_inout_control_bondstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Set the formula(s) LAMMPS uses to compute bond interactions between pairs of atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_bondstyle'))

    x_lammps_inout_control_bondwrite = Quantity(
        type=str,
        shape=[],
        description='''
        Write energy and force values to a file as a function of distance for the
        currently defined bond potential.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_bondwrite'))

    x_lammps_inout_control_boundary = Quantity(
        type=str,
        shape=[],
        description='''
        Set the style of boundaries for the global simulation box in each dimension.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_boundary'))

    x_lammps_inout_control_box = Quantity(
        type=str,
        shape=[],
        description='''
        Set attributes of the simulation box.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_box'))

    x_lammps_inout_control_changebox = Quantity(
        type=str,
        shape=[],
        description='''
        Change the volume and/or shape and/or boundary conditions for the simulation box.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_change_box'))

    x_lammps_inout_control_clear = Quantity(
        type=str,
        shape=[],
        description='''
        This command deletes all atoms, restores all settings to their default values,
        and frees all memory allocated by LAMMPS.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_clear'))

    x_lammps_inout_control_commmodify = Quantity(
        type=str,
        shape=[],
        description='''
        This command sets parameters that affect the inter-processor communication of atom
        information that occurs each timestep as coordinates and other properties are
        exchanged between neighboring processors and stored as properties of ghost atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_commmodify'))

    x_lammps_inout_control_commstyle = Quantity(
        type=str,
        shape=[],
        description='''
        This command sets the style of inter-processor communication of atom information
        that occurs each timestep as coordinates and other properties are exchanged
        between neighboring processors and stored as properties of ghost atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_commstyle'))

    x_lammps_inout_control_compute = Quantity(
        type=str,
        shape=[],
        description='''
        Define a computation that will be performed on a group of atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_compute'))

    x_lammps_inout_control_computemodify = Quantity(
        type=str,
        shape=[],
        description='''
        Modify one or more parameters of a previously defined compute.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_computemodify'))

    x_lammps_inout_control_createatoms = Quantity(
        type=str,
        shape=[],
        description='''
        This command creates atoms (or molecules) on a lattice, or a single atom
        (or molecule), or a random collection of atoms (or molecules), as an alternative
        to reading in their coordinates explicitly via a read_data or read_restart command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_createatoms'))

    x_lammps_inout_control_createbonds = Quantity(
        type=str,
        shape=[],
        description='''
        Create bonds between pairs of atoms that meet a specified distance criteria.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_createbonds'))

    x_lammps_inout_control_createbox = Quantity(
        type=str,
        shape=[],
        description='''
        This command creates a simulation box based on the specified region.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_createbox'))

    x_lammps_inout_control_deleteatoms = Quantity(
        type=str,
        shape=[],
        description='''
        Delete the specified atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_deleteatoms'))

    x_lammps_inout_control_deletebonds = Quantity(
        type=str,
        shape=[],
        description='''
        Turn off (or on) molecular topology interactions, i.e. bonds, angles, dihedrals,
        impropers.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_deletebonds'))

    x_lammps_inout_control_dielectric = Quantity(
        type=str,
        shape=[],
        description='''
        Set the dielectric constant for Coulombic interactions (pairwise and long-range)
        to this value.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_dielectric'))

    x_lammps_inout_control_dihedralcoeff = Quantity(
        type=str,
        shape=[],
        description='''
        Specify the dihedral force field coefficients for one or more dihedral types.
        ''',
        a_legacy=LegacyDefinition(name=''))

    x_lammps_inout_control_dihedralstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Set the formula(s) LAMMPS uses to compute dihedral interactions between quadruplets
        of atoms, which remain in force for the duration of the simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_dihedralstyle'))

    x_lammps_inout_control_dimension = Quantity(
        type=str,
        shape=[],
        description='''
        Set the dimensionality of the simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_dimension'))

    x_lammps_inout_control_displaceatoms = Quantity(
        type=str,
        shape=[],
        description='''
        Displace a group of atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_displaceatoms'))

    x_lammps_inout_control_dump = Quantity(
        type=str,
        shape=[],
        description='''
        Dump a snapshot of atom quantities to one or more files every N timesteps in one
        of several styles.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_dump'))

    x_lammps_inout_control_dynamicalmatrix = Quantity(
        type=str,
        shape=[],
        description='''
        Calculate the dynamical matrix by finite difference of the selected group.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_dynamicalmatrix'))

    x_lammps_inout_control_echo = Quantity(
        type=str,
        shape=[],
        description='''
        This command determines whether LAMMPS echoes each input script command to the
        screen and/or log file as it is read and processed.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_echo'))

    x_lammps_inout_control_fix = Quantity(
        type=str,
        shape=[],
        description='''
        Set a fix that will be applied to a group of atoms. In LAMMPS, a “fix” is any
        operation that is applied to the system during timestepping or minimization
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_fix'))

    x_lammps_inout_control_fixmodify = Quantity(
        type=str,
        shape=[],
        description='''
        Modify one or more parameters of a previously defined fix.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_fixmodify'))

    x_lammps_inout_control_group = Quantity(
        type=str,
        shape=[],
        description='''
        Identify a collection of atoms as belonging to a group.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_group'))

    x_lammps_inout_control_group2ndx = Quantity(
        type=str,
        shape=[],
        description='''
        Write or read a Gromacs style index file in text format that associates atom IDs
        with the corresponding group definitions. The group2ndx command will write group
        definitions to an index file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_group2ndx'))

    x_lammps_inout_control_ndx2group = Quantity(
        type=str,
        shape=[],
        description='''
        Write or read a Gromacs style index file in text format that associates atom IDs
        with the corresponding group definitions. The ndx2group command will create of
        update group definitions from those stored in an index file.
        ''',
        a_legacy=LegacyDefinition(name=''))

    x_lammps_inout_control_hyper = Quantity(
        type=str,
        shape=[],
        description='''
        Run a bond-boost hyperdynamics (HD) simulation where time is accelerated by
        application of a bias potential to one or more pairs of nearby atoms in the system.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_hyper'))

    x_lammps_inout_control_if = Quantity(
        type=str,
        shape=[],
        description='''
        This command provides an if-then-else capability within an input script.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_if'))

    x_lammps_inout_control_impropercoeff = Quantity(
        type=str,
        shape=[],
        description='''
        Specify the improper force field coefficients for one or more improper types.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_impropercoeff'))

    x_lammps_inout_control_improperstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Set the formula(s) LAMMPS uses to compute improper interactions between
        quadruplets of atoms, which remain in force for the duration of the simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_improperstyle'))

    x_lammps_inout_control_include = Quantity(
        type=str,
        shape=[],
        description='''
        This command opens a new input script file and begins reading LAMMPS commands
        from that file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_include'))

    x_lammps_inout_control_info = Quantity(
        type=str,
        shape=[],
        description='''
        Print out information about the current internal state of the running LAMMPS process.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_info'))

    x_lammps_inout_control_jump = Quantity(
        type=str,
        shape=[],
        description='''
        This command closes the current input script file, opens the file with the
        specified name, and begins reading LAMMPS commands from that file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_jump'))

    x_lammps_inout_control_kiminit = Quantity(
        type=str,
        shape=[],
        description='''
        The set of kim_commands provide a high-level wrapper around the Open Knowledgebase
        of Interatomic Models (OpenKIM) repository of interatomic models (IMs)
        (potentials and force fields), so that they can be used by LAMMPS scripts.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kiminit'))

    x_lammps_inout_control_kiminteractions = Quantity(
        type=str,
        shape=[],
        description='''
        The set of kim_commands provide a high-level wrapper around the Open Knowledgebase
        of Interatomic Models (OpenKIM) repository of interatomic models (IMs)
        (potentials and force fields), so that they can be used by LAMMPS scripts.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kiminteractions'))

    x_lammps_inout_control_kimquery = Quantity(
        type=str,
        shape=[],
        description='''
        The set of kim_commands provide a high-level wrapper around the Open Knowledgebase
        of Interatomic Models (OpenKIM) repository of interatomic models (IMs)
        (potentials and force fields), so that they can be used by LAMMPS scripts.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kimquery'))

    x_lammps_inout_control_kimparam = Quantity(
        type=str,
        shape=[],
        description='''
        The set of kim_commands provide a high-level wrapper around the Open Knowledgebase
        of Interatomic Models (OpenKIM) repository of interatomic models (IMs)
        (potentials and force fields), so that they can be used by LAMMPS scripts.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kimparam'))

    x_lammps_inout_control_kimproperty = Quantity(
        type=str,
        shape=[],
        description='''
        The set of kim_commands provide a high-level wrapper around the Open Knowledgebase
        of Interatomic Models (OpenKIM) repository of interatomic models (IMs)
        (potentials and force fields), so that they can be used by LAMMPS scripts.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kimproperty'))

    x_lammps_inout_control_kspacemodify = Quantity(
        type=str,
        shape=[],
        description='''
        Set parameters used by the kspace solvers defined by the kspace_style command.
        Not all parameters are relevant to all kspace styles.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kspacemodify'))

    x_lammps_inout_control_kspacestyle = Quantity(
        type=str,
        shape=[],
        description='''
        Define a long-range solver for LAMMPS to use each timestep to compute long-range
        Coulombic interactions or long-range interactions.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_kspacestyle'))

    x_lammps_inout_control_label = Quantity(
        type=str,
        shape=[],
        description='''
        Label this line of the input script with the chosen ID.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_label'))

    x_lammps_inout_control_lattice = Quantity(
        type=str,
        shape=[],
        description='''
        Define a lattice for use by other commands.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_lattice'))

    x_lammps_inout_control_log = Quantity(
        type=str,
        shape=[],
        description='''
        This command closes the current LAMMPS log file, opens a new file with the
        specified name, and begins logging information to it.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_log'))

    x_lammps_inout_control_mass = Quantity(
        type=str,
        shape=[],
        description='''
        Set the mass for all atoms of one or more atom types.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_mass'))

    x_lammps_inout_control_message = Quantity(
        type=str,
        shape=[],
        description='''
        Establish a messaging protocol between LAMMPS and another code for the purpose of
        client/server coupling.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_message'))

    x_lammps_inout_control_minmodify = Quantity(
        type=str,
        shape=[],
        description='''
        This command sets parameters that affect the energy minimization algorithms
        selected by the min_style command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_minmodify'))

    x_lammps_inout_control_minstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Apply a minimization algorithm to use when a minimize command is performed.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_minstyle'))

    x_lammps_inout_control_minimize = Quantity(
        type=str,
        shape=[],
        description='''
        Perform an energy minimization of the system, by iteratively adjusting atom
        coordinates.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_minimize'))

    x_lammps_inout_control_molecule = Quantity(
        type=str,
        shape=[],
        description='''
        Define a molecule template that can be used as part of other LAMMPS commands,
        typically to define a collection of particles as a bonded molecule or a rigid body.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_molecule'))

    x_lammps_inout_control_neb = Quantity(
        type=str,
        shape=[],
        description='''
        Perform a nudged elastic band (NEB) calculation using multiple replicas of a system.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_neb'))

    x_lammps_inout_control_neighmodify = Quantity(
        type=str,
        shape=[],
        description='''
        This command sets parameters that affect the building and use of pairwise neighbor
        lists.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_neighmodify'))

    x_lammps_inout_control_neighbor = Quantity(
        type=str,
        shape=[],
        description='''
        This command sets parameters that affect the building of pairwise neighbor lists.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_neighbor'))

    x_lammps_inout_control_newton = Quantity(
        type=str,
        shape=[],
        description='''
        This command turns Newton’s third law on or off for pairwise and bonded interactions.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_newton'))

    x_lammps_inout_control_next = Quantity(
        type=str,
        shape=[],
        description='''
        This command is used with variables defined by the variable command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_next'))

    x_lammps_inout_control_package = Quantity(
        type=str,
        shape=[],
        description='''
        This command invokes package-specific settings for the various accelerator
        packages available in LAMMPS.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_package'))

    x_lammps_inout_control_paircoeff = Quantity(
        type=str,
        shape=[],
        description='''
        Specify the pairwise force field coefficients for one or more pairs of
        atom types.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_paircoeff'))

    x_lammps_inout_control_pairmodify = Quantity(
        type=str,
        shape=[],
        description='''
        Modify the parameters of the currently defined pair style.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_pairmodify'))

    x_lammps_inout_control_pairstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Set the formula(s) LAMMPS uses to compute pairwise interactions.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_pairstyle'))

    x_lammps_inout_control_pairwrite = Quantity(
        type=str,
        shape=[],
        description='''
        Write energy and force values to a file as a function of distance for
        the currently defined pair potential.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_pairwrite'))

    x_lammps_inout_control_partition = Quantity(
        type=str,
        shape=[],
        description='''
        This command invokes the specified command on a subset of the
        partitions of processors you have defined via the -partition command-line switch.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_partition'))

    x_lammps_inout_control_prd = Quantity(
        type=str,
        shape=[],
        description='''
        Run a parallel replica dynamics (PRD) simulation using multiple
        replicas of a system.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_prd'))

    x_lammps_inout_control_print = Quantity(
        type=str,
        shape=[],
        description='''
        Print a text string to the screen and logfile.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_print'))

    x_lammps_inout_control_processors = Quantity(
        type=str,
        shape=[],
        description='''
        Specify how processors are mapped as a regular 3d grid to the global
        simulation box.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_processors'))

    x_lammps_inout_control_quit = Quantity(
        type=str,
        shape=[],
        description='''
        This command causes LAMMPS to exit, after shutting down all output
        cleanly.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_quit'))

    x_lammps_inout_control_readdata = Quantity(
        type=str,
        shape=[],
        description='''
        Read in a data file containing information LAMMPS needs to run a
        simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_readdata'))

    x_lammps_inout_control_readdump = Quantity(
        type=str,
        shape=[],
        description='''
        Read atom information from a dump file to overwrite the current atom
        coordinates, and optionally the atom velocities and image flags and
        the simulation box dimensions.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_readdump'))

    x_lammps_inout_control_readrestart = Quantity(
        type=str,
        shape=[],
        description='''
        Read in a previously saved system configuration from a restart file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_readrestart'))

    x_lammps_inout_control_region = Quantity(
        type=str,
        shape=[],
        description='''
        This command defines a geometric region of space.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_region'))

    x_lammps_inout_control_replicate = Quantity(
        type=str,
        shape=[],
        description='''
        Replicate the current simulation one or more times in each dimension.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_replicate'))

    x_lammps_inout_control_rerun = Quantity(
        type=str,
        shape=[],
        description='''
        Perform a pseudo simulation run where atom information is read one
        snapshot at a time from a dump file(s), and energies and forces are
        computed on the shapshot to produce thermodynamic or other output.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_rerun'))

    x_lammps_inout_control_resetatomids = Quantity(
        type=str,
        shape=[],
        description='''
        Reset atom IDs for the system, including all the global IDs stored
        for bond, angle, dihedral, improper topology data.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_resetatomids'))

    x_lammps_inout_control_resetmolids = Quantity(
        type=str,
        shape=[],
        description='''
        Reset molecule IDs for a group of atoms based on current bond
        connectivity.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_resetmolids'))

    x_lammps_inout_control_resettimestep = Quantity(
        type=str,
        shape=[],
        description='''
        Set the timestep counter to the specified value.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_resettimestep'))

    x_lammps_inout_control_restart = Quantity(
        type=str,
        shape=[],
        description='''
        Write out a binary restart file with the current state of the
        simulation every so many timesteps, in either or both of two modes, as
        a run proceeds.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_restart'))

    x_lammps_inout_control_run = Quantity(
        type=str,
        shape=[],
        description='''
        Run or continue dynamics for a specified number of timesteps.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_run'))

    x_lammps_inout_control_runstyle = Quantity(
        type=str,
        shape=[],
        description='''
        Choose the style of time integrator used for molecular dynamics
        simulations performed by LAMMPS.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_runstyle'))

    x_lammps_inout_control_server = Quantity(
        type=str,
        shape=[],
        description='''
        This command starts LAMMPS running in “server” mode, where it receives
        messages from a separate “client” code and responds by sending a reply
        message back to the client.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_server'))

    x_lammps_inout_control_set = Quantity(
        type=str,
        shape=[],
        description='''
        Set one or more properties of one or more atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_set'))

    x_lammps_inout_control_shell = Quantity(
        type=str,
        shape=[],
        description='''
        Execute a shell command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_shell'))

    x_lammps_inout_control_specialbonds = Quantity(
        type=str,
        shape=[],
        description='''
        Set weighting coefficients for pairwise energy and force contributions
        between pairs of atoms that are also permanently bonded to each other,
        either directly or via one or two intermediate bonds.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_specialbonds'))

    x_lammps_inout_control_suffix = Quantity(
        type=str,
        shape=[],
        description='''
        This command allows you to use variants of various styles if they
        exist.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_suffix'))

    x_lammps_inout_control_tad = Quantity(
        type=str,
        shape=[],
        description='''
        Run a temperature accelerated dynamics (TAD) simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_tad'))

    x_lammps_inout_control_tempergrem = Quantity(
        type=str,
        shape=[],
        description='''
        Run a parallel tempering or replica exchange simulation in LAMMPS
        partition mode using multiple generalized replicas (ensembles) of a
        system defined by fix grem, which stands for the
        generalized replica exchange method (gREM) originally developed by
        (Kim).
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_tempergrem'))

    x_lammps_inout_control_tempernpt = Quantity(
        type=str,
        shape=[],
        description='''
        Run a parallel tempering or replica exchange simulation using multiple
        replicas (ensembles) of a system in the isothermal-isobaric (NPT)
        ensemble.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_tempernpt'))

    x_lammps_inout_control_thermo = Quantity(
        type=str,
        shape=[],
        description='''
        Compute and print thermodynamic info (e.g. temperature, energy,
        pressure) on timesteps that are a multiple of N and at the beginning
        and end of a simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_thermo'))

    x_lammps_inout_control_thermomodify = Quantity(
        type=str,
        shape=[],
        description='''
        Set options for how thermodynamic information is computed and printed
        by LAMMPS.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_thermomodify'))

    x_lammps_inout_control_thermostyle = Quantity(
        type=str,
        shape=[],
        description='''
        Set the style and content for printing thermodynamic data to the
        screen and log file.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_thermostyle'))

    x_lammps_inout_control_thirdorder = Quantity(
        type=str,
        shape=[],
        description='''
        Calculate the third order force constant tensor by finite difference of the selected group,

        where Phi is the third order force constant tensor.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_thirdorder'))

    x_lammps_inout_control_timer = Quantity(
        type=str,
        shape=[],
        description='''
        Select the level of detail at which LAMMPS performs its CPU timings.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_timer'))

    x_lammps_inout_control_timestep = Quantity(
        type=str,
        shape=[],
        description='''
        Set the timestep size for subsequent molecular dynamics simulations.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_timestep'))

    x_lammps_inout_control_uncompute = Quantity(
        type=str,
        shape=[],
        description='''
        Delete a compute that was previously defined with a compute
        command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_uncompute'))

    x_lammps_inout_control_undump = Quantity(
        type=str,
        shape=[],
        description='''
        Turn off a previously defined dump so that it is no longer active.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_undump'))

    x_lammps_inout_control_unfix = Quantity(
        type=str,
        shape=[],
        description='''
        Delete a fix that was previously defined with a fix
        command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_unfix'))

    x_lammps_inout_control_units = Quantity(
        type=str,
        shape=[],
        description='''
        This command sets the style of units used for a simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_units'))

    x_lammps_inout_control_variable = Quantity(
        type=str,
        shape=[],
        description='''
        This command assigns one or more strings to a variable name for
        evaluation later in the input script or during a simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_variable'))

    x_lammps_inout_control_velocity = Quantity(
        type=str,
        shape=[],
        description='''
        Set or change the velocities of a group of atoms in one of several
        styles.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_velocity'))

    x_lammps_inout_control_writecoeff = Quantity(
        type=str,
        shape=[],
        description='''
        Write a text format file with the currently defined force field
        coefficients in a way, that it can be read by LAMMPS with the
        include command.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_writecoeff'))

    x_lammps_inout_control_writedata = Quantity(
        type=str,
        shape=[],
        description='''
        Write a data file in text format of the current state of the
        simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_writedata'))

    x_lammps_inout_control_writedump = Quantity(
        type=str,
        shape=[],
        description='''
        Dump a single snapshot of atom quantities to one or more files for the
        current state of the system.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_writedump'))

    x_lammps_inout_control_writerestart = Quantity(
        type=str,
        shape=[],
        description='''
        Write a binary restart file of the current state of the simulation.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_inout_control_writerestart'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_lammps_section_input_output_files = SubSection(
        sub_section=SectionProxy('x_lammps_section_input_output_files'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_lammps_section_input_output_files'))

    x_lammps_section_control_parameters = SubSection(
        sub_section=SectionProxy('x_lammps_section_control_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_lammps_section_control_parameters'))


class section_topology(common.section_topology):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_topology'))

    x_lammps_xlo_xhi = Quantity(
        type=str,
        shape=[],
        description='''
        test
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_xlo_xhi'))

    x_lammps_data_file_store = Quantity(
        type=str,
        shape=[],
        description='''
        Filename of data file
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_file_store'))

    x_lammps_dummy = Quantity(
        type=str,
        shape=[],
        description='''
        dummy
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_dummy'))

    x_lammps_input_units_store = Quantity(
        type=str,
        shape=[],
        description='''
        It determines the units of all quantities specified in the input script and data
        file, as well as quantities output to the screen, log file, and dump files.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_input_units_store'))

    x_lammps_data_bd_types_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_bd_types_store'))

    x_lammps_data_bd_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_bd_count_store'))

    x_lammps_data_ag_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_ag_count_store'))

    x_lammps_data_at_types_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_at_types_store'))

    x_lammps_data_dh_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_dh_count_store'))

    x_lammps_data_angles_store = Quantity(
        type=str,
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_angles_store'))

    x_lammps_data_angle_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_angle_list_store'))

    x_lammps_data_bond_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_bond_list_store'))

    x_lammps_data_dihedral_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_dihedral_list_store'))

    x_lammps_data_dihedral_coeff_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_dihedral_coeff_list_store'))

    x_lammps_masses_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_masses_store'))

    x_lammps_data_topo_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_data_topo_list_store'))


m_package.__init_metainfo__()
