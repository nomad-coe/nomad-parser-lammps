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


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_lammps_program_version_date = Quantity(
        type=str,
        shape=[],
        description='''
        Program version date.
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_program_version_date'))

    x_lammps_dummy_text = Quantity(
        type=str,
        shape=[],
        description='''
        dummy text
        ''',
        a_legacy=LegacyDefinition(name='x_lammps_dummy_text'))


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
