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
import numpy as np
import os
import logging
from ase import data as asedata
import pint

from .metainfo import m_env
from nomad.parsing.parser import FairdiParser

from nomad.parsing.file_parser import Quantity, TextParser
from nomad.datamodel.metainfo.common_dft import Run, SamplingMethod, System,\
    SingleConfigurationCalculation, EnergyContribution, Workflow, MolecularDynamics
from nomad.datamodel.metainfo.common import section_topology, section_interaction
from .metainfo.lammps import x_lammps_section_input_output_files, x_lammps_section_control_parameters


def get_unit(units_type, property_type=None, dimension=3):
    mole = 6.022140857e+23

    units_type = units_type.lower()
    if units_type == 'real':
        units = dict(
            mass=pint.Quantity(1 / mole, 'g'), distance='angstrom', time='fs',
            energy=pint.Quantity(1 / mole, 'kcal'), velocity='angstrom/fs',
            force=pint.Quantity(1 / mole, 'kcal/angstrom'), torque=pint.Quantity(1 / mole, 'kcal'),
            temperature='K', pressure='atm', dynamic_viscosity='poise', charge='elementary_charge',
            dipole='elementary_charge*angstrom', electric_field='V/angstrom',
            density='g/cm^%d' % dimension)

    elif units_type == 'metal':
        units = dict(
            mass=pint.Quantity(1 / mole, 'g'), distance='angstrom', time='ps',
            energy='eV', velocity='angstrom/ps', force='eV/angstrom', torque='eV',
            temperature='K', pressure='bar', dynamic_viscosity='poise', charge='elementary_charge',
            dipole='elementary_charge*angstrom', electric_field='V/angstrom',
            density='g/cm^%d' % dimension)

    elif units_type == 'si':
        units = dict(
            mass='kg', distance='m', time='s', energy='J', velocity='m/s', force='N',
            torque='N*m', temperature='K', pressure='Pa', dynamic_viscosity='Pa*s',
            charge='C', dipole='C*m', electric_field='V/m', density='kg/m^%d' % dimension)

    elif units_type == 'cgs':
        units = dict(
            mass='g', distance='cm', time='s', energy='erg', velocity='cm/s', force='dyne',
            torque='dyne*cm', temperature='K', pressure='dyne/cm^2', dynamic_viscosity='poise',
            charge='esu', dipole='esu*cm', electric_field='dyne/esu',
            density='g/cm^%d' % dimension)

    elif units_type == 'electron':
        units = dict(
            mass='amu', distance='bohr', time='fs', energy='hartree',
            velocity='bohr/atomic_unit_of_time', force='hartree/bohr', temperature='K',
            pressure='Pa', charge='elementary_charge', dipole='debye', electric_field='V/cm')

    elif units_type == 'micro':
        units = dict(
            mass='pg', distance='microm', time='micros', energy='pg*microm^2/micros^2',
            velocity='microm/micros', force='pg*microm/micros^2', torque='pg*microm^2/micros^2',
            temperature='K', pressure='pg/(microm*micros^2)', dynamic_viscosity='pg/(microm*micros)',
            charge='pC', dipole='pC*microm', electric_field='V/microm',
            density='pg/microm^%d' % dimension)

    elif units_type == 'nano':
        units = dict(
            mass='ag', distance='nm', time='ns', energy='ag*nm^2/ns^2', velocity='nm/ns',
            force='ag*nm/ns^2', torque='ag*nm^2/ns^2', temperature='K', pressure='ag/(nm*ns^2)',
            dynamic_viscosity='ag/(nm*ns)', charge='elementary_charge',
            dipole='elementary_charge*nm', electric_field='V/nm', density='ag/nm^%d' % dimension)

    else:
        units = dict()

    if property_type:
        unit = units.get(property_type, None)
        if isinstance(unit, str):
            unit = pint.Quantity(1, unit)
        return unit
    else:
        for key, val in units.items():
            if isinstance(val, str):
                units[key] = pint.Quantity(1, val)
        return units


class DataParser(TextParser):
    def __init__(self):
        self._headers = [
            'atoms', 'bonds', 'angles', 'dihedrals', 'impropers', 'atom types', 'bond types',
            'angle types', 'dihedral types', 'improper types', 'extra bond per atom',
            'extra/bond/per/atom', 'extra angle per atom', 'extra/angle/per/atom',
            'extra dihedral per atom', 'extra/dihedral/per/atom', 'extra improper per atom',
            'extra/improper/per/atom', 'extra special per atom', 'extra/special/per/atom',
            'ellipsoids', 'lines', 'triangles', 'bodies', 'xlo xhi', 'ylo yhi', 'zlo zhi',
            'xy xz yz']
        self._sections = [
            'Atoms', 'Velocities', 'Masses', 'Ellipsoids', 'Lines', 'Triangles', 'Bodies',
            'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Pair Coeffs', 'PairIJ Coeffs',
            'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs',
            'BondBond Coeffs', 'BondAngle Coeffs', 'MiddleBondTorsion Coeffs',
            'EndBondTorsion Coeffs', 'AngleTorsion Coeffs', 'AngleAngleTorsion Coeffs',
            'BondBond13 Coeffs', 'AngleAngle Coeffs']
        self._interactions = [
            section for section in self._sections if section.endswith('Coeffs')]
        super().__init__(None)

    def init_quantities(self):
        self._quantities = []
        for header in self._headers:
            self._quantities.append(Quantity(
                header, r'\s*([\+\-eE\d\. ]+)\s*%s\s*\n' % header, comment='#', repeats=True))

        def get_section_value(val):
            val = val.split('\n')
            name = None

            if val[0][0] == '#':
                name = val[0][1:].strip()
                val = val[1:]

            value = []
            for i in range(len(val)):
                v = val[i].split('#')[0].split()
                if not v:
                    continue

                try:
                    value.append(np.array(v, dtype=float))
                except Exception:
                    break

            return name, np.array(value)

        for section in self._sections:
            self._quantities.append(
                Quantity(
                    section, r'\s*%s\s*(#*\s*[\s\S]*?\n)\n*([\deE\-\+\.\s]+)\n' % section,
                    str_operation=get_section_value, repeats=True))

    def get_interactions(self):
        styles_coeffs = []
        for interaction in self._interactions:
            coeffs = self.get(interaction, None)
            if coeffs is None:
                continue
            if isinstance(coeffs, tuple):
                coeffs = [coeffs]

            styles_coeffs += coeffs

        return styles_coeffs


class TrajParser(TextParser):
    def __init__(self):
        self._masses = None
        self._reference_masses = dict(
            masses=np.array(asedata.atomic_masses), symbols=asedata.chemical_symbols)
        self._chemical_symbols = None
        self._units = None
        super().__init__(None)

    def init_quantities(self):

        def get_pbc_cell(val):
            val = val.split()
            pbc = [v == 'pp' for v in val[:3]]

            cell = np.zeros((3, 3))
            for i in range(3):
                cell[i][i] = float(val[i * 2 + 4]) - float(val[i * 2 + 3])

            return pbc, cell

        def get_atoms_info(val):
            val = val.split('\n')
            keys = val[0].split()
            values = np.array([v.split() for v in val[1:] if v], dtype=float)
            values = values[values[:, 0].argsort()].T
            return {keys[i]: values[i] for i in range(len(keys))}

        self._quantities = [
            Quantity(
                'time_step', r'\s*ITEM:\s*TIMESTEP\s*\n\s*(\d+)\s*\n', comment='#',
                repeats=True),
            Quantity(
                'n_atoms', r'\s*ITEM:\s*NUMBER OF ATOMS\s*\n\s*(\d+)\s*\n', comment='#',
                repeats=True),
            Quantity(
                'pbc_cell', r'\s*ITEM: BOX BOUNDS\s*([\s\w]+)([\+\-\d\.eE\s]+)\n',
                str_operation=get_pbc_cell, comment='#', repeats=True),
            Quantity(
                'atoms_info', r's*ITEM:\s*ATOMS\s*([ \w]+\n)*?([\+\-eE\d\.\n ]+)\n*I*',
                str_operation=get_atoms_info, comment='#', repeats=True)
        ]

    def with_trajectory(self):
        return self.get('atoms_info') is not None

    @property
    def masses(self):
        return self._masses

    @masses.setter
    def masses(self, val):
        self._masses = val
        if self._masses is None:
            return

        self._masses = val
        if self._chemical_symbols is None:
            masses = self._masses[0][1]
            self._chemical_symbols = {}
            for i in range(len(masses)):
                symbol_idx = np.argmin(abs(self._reference_masses['masses'] - masses[i][1]))
                self._chemical_symbols[masses[i][0]] = self._reference_masses['symbols'][symbol_idx]

    def get_atom_labels(self, idx):
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return

        atoms_type = atoms_info[idx].get('type')
        if atoms_type is None:
            return

        if self._chemical_symbols is None:
            return

        atom_labels = [self._chemical_symbols[atype] for atype in atoms_type]

        return atom_labels

    def get_positions(self, idx):
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return

        atoms_info = atoms_info[idx]

        cell = self.get('pbc_cell')
        if cell is None:
            return

        cell = cell[idx][1]

        if 'xs' in atoms_info and 'ys' in atoms_info and 'zs' in atoms_info:
            positions = np.array([atoms_info['xs'], atoms_info['ys'], atoms_info['zs']]).T

            positions = positions * np.linalg.norm(cell, axis=1) + np.amin(cell, axis=1)

        else:
            positions = np.array([atoms_info['x'], atoms_info['y'], atoms_info['z']]).T

            if 'ix' in atoms_info and 'iy' in atoms_info and 'iz' in atoms_info:
                positions_img = np.array([
                    atoms_info['ix'], atoms_info['iy'], atoms_info['iz']]).T

                positions += positions_img * np.linalg.norm(cell, axis=1)

        return pint.Quantity(positions, self._units.get('distance', None))

    def get_velocities(self, idx):
        atoms_info = self.get('atoms_info')

        if atoms_info is None:
            return

        atoms_info = atoms_info[idx]

        if 'vx' not in atoms_info or 'vy' not in atoms_info or 'vz' not in atoms_info:
            return

        velocities = np.array([atoms_info['vx'], atoms_info['vy'], atoms_info['vz']]).T

        return pint.Quantity(velocities, self._units.get('velocity', None))

    def get_forces(self, idx):
        atoms_info = self.get('atoms_info')

        if atoms_info is None:
            return

        atoms_info = atoms_info[idx]

        if 'fx' not in atoms_info or 'fy' not in atoms_info or 'fz' not in atoms_info:
            return

        forces = np.array([atoms_info['fx'], atoms_info['fy'], atoms_info['fz']]).T

        return pint.Quantity(forces, self._units.get('force', None))


class LogParser(TextParser):
    def __init__(self):
        self._commands = [
            'angle_coeff', 'angle_style', 'atom_modify', 'atom_style', 'balance',
            'bond_coeff', 'bond_style', 'bond_write', 'boundary', 'change_box', 'clear',
            'comm_modify', 'comm_style', 'compute', 'compute_modify', 'create_atoms',
            'create_bonds', 'create_box', 'delete_bonds', 'dielectric', 'dihedral_coeff',
            'dihedral_style', 'dimension', 'displace_atoms', 'dump', 'dump_modify',
            'dynamical_matrix', 'echo', 'fix', 'fix_modify', 'group', 'group2ndx',
            'ndx2group', 'hyper', 'if', 'improper_coeff', 'improper_style', 'include',
            'info', 'jump', 'kim_init', 'kim_interactions', 'kim_query', 'kim_param',
            'kim_property', 'kspace_modify', 'kspace_style', 'label', 'lattice', 'log',
            'mass', 'message', 'min_modify', 'min_style', 'minimize', 'minimize/kk',
            'molecule', 'neb', 'neb/spin', 'neigh_modify', 'neighbor', 'newton', 'next',
            'package', 'pair_coeff', 'pair_modify', 'pair_style', 'pair_write',
            'partition', 'prd', 'print', 'processors', 'quit', 'read_data', 'read_dump',
            'read_restart', 'region', 'replicate', 'rerun', 'reset_atom_ids',
            'reset_mol_ids', 'reset_timestep', 'restart', 'run', 'run_style', 'server',
            'set', 'shell', 'special_bonds', 'suffix', 'tad', 'temper/grem', 'temper/npt',
            'thermo', 'thermo_modify', 'thermo_style', 'third_order', 'timer', 'timestep',
            'uncompute', 'undump', 'unfix', 'units', 'variable', 'velocity', 'write_coeff',
            'write_data', 'write_dump', 'write_restart']
        self._interactions = [
            'atom', 'pair', 'bond', 'angle', 'dihedral', 'improper', 'kspace']
        self._thermo_data = None
        self._units = None
        super().__init__(None)

    def init_quantities(self):
        def str_op(val):
            val = val.split('#')[0]
            val = val.replace('&\n', ' ').split()
            val = val if len(val) > 1 else val[0]
            return val

        self._quantities = [
            Quantity(
                name, r'\n\s*%s\s+([\w\. \/\#\-]+)(\&\n[\w\. \/\#\-]*)*' % name,
                str_operation=str_op, comment='#', repeats=True) for name in self._commands]

        self._quantities.append(Quantity(
            'program_version', r'\s*LAMMPS\s*\(([\w ]+)\)\n', dtype=str, repeats=False,
            flatten=False)
        )

        self._quantities.append(Quantity(
            'finished', r'\s*Dangerous builds\s*=\s*(\d+)', repeats=False)
        )

        def str_to_thermo(val):
            res = {}
            if val.count('Step') > 1:
                val = val.replace('-', '').replace('=', '').replace('(sec)', '').split()
                val = [v.strip() for v in val]

                for i in range(len(val)):
                    if val[i][0].isalpha():
                        res.setdefault(val[i], [])
                        res[val[i]].append(float(val[i + 1]))

            else:
                val = val.split('\n')
                keys = [v.strip() for v in val[0].split()]
                val = np.array([v.split() for v in val[1:] if v], dtype=float).T

                res = {key: [] for key in keys}
                for i in range(len(keys)):
                    res[keys[i]] = val[i]

            return res

        self._quantities.append(Quantity(
            'thermo_data', r'\s*\-*(\s*Step\s*[\-\s\w\.\=\(\)]*[ \-\.\d\n]+)Loop',
            str_operation=str_to_thermo, repeats=False, convert=False)
        )

    @property
    def units(self):
        if self._file_handler is None or self._units is None:
            units_type = self.get('units', ['lj'])[0]
            self._units = get_unit(units_type)
        return self._units

    def get_thermodynamic_data(self):
        self._thermo_data = self.get('thermo_data')

        if self._thermo_data is None:
            return

        for key, val in self._thermo_data.items():
            low_key = key.lower()
            if low_key.startswith('e_') or low_key.endswith('eng'):
                self._thermo_data[key] = val * self.units['energy']
            elif low_key == 'press':
                self._thermo_data[key] = val * self.units['pressure']
            elif low_key == 'temp':
                self._thermo_data[key] = val * self.units['temperature']

        return self._thermo_data

    def get_traj_files(self):
        dump = self.get('dump')
        if dump is None:
            self.logger.warn(
                'Trajectory not specified in directory, will scan.',
                data=dict(directory=self.maindir))
            traj_files = os.listdir(self.maindir)
            traj_files = [f for f in traj_files if f.endswith('trj')]

        else:
            traj_files = []
            if type(dump[0]) in [str, int]:
                dump = [dump]
            traj_files = [d[4] for d in dump]

        return [os.path.join(self.maindir, f) for f in traj_files]

    def get_data_files(self):
        read_data = self.get('read_data')
        if read_data is None:
            self.logger.warn(
                'Data file not specified in directory, will scan.',
                data=dict(directory=self.maindir))
            data_files = os.listdir(self.maindir)
            data_files = [f for f in data_files if f.endswith('data') or f.startswith('data')]

        else:
            data_files = read_data

        return [os.path.join(self.maindir, f) for f in data_files]

    def get_pbc(self):
        pbc = self.get('boundary', ['p', 'p', 'p'])
        return [v == 'p' for v in pbc]

    def get_sampling_method(self):
        fix_style = self.get('fix', [[''] * 3])[0][2]

        sampling_method = 'langevin_dynamics' if 'langevin' in fix_style else 'molecular_dynamics'
        return sampling_method, fix_style

    def get_thermostat_settings(self):
        fix = self.get('fix', [None])[0]
        if fix is None:
            return {}

        try:
            fix_style = fix[2]
        except IndexError:
            return {}

        temp_unit = self.units['temperature']
        press_unit = self.units['pressure']
        time_unit = self.units['time']

        res = dict()
        if fix_style.lower() == 'nvt':
            try:
                res['target_T'] = float(fix[5]) * temp_unit
                res['thermostat_tau'] = float(fix[6]) * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'npt':
            try:
                res['target_T'] = float(fix[5]) * temp_unit
                res['thermostat_tau'] = float(fix[6]) * time_unit
                res['target_P'] = float(fix[9]) * press_unit
                res['barostat_tau'] = float(fix[10]) * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'nph':
            try:
                res['target_P'] = float(fix[5]) * press_unit
                res['barostat_tau'] = float(fix[6]) * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'langevin':
            try:
                res['target_T'] = float(fix[4]) * temp_unit
                res['langevin_gamma'] = float(fix[5]) * time_unit
            except Exception:
                pass

        else:
            self.logger.warn('Fix style not supported', data=dict(style=fix_style))

        return res

    def get_interactions(self):
        styles_coeffs = []
        for interaction in self._interactions:
            styles = self.get('%s_style' % interaction, None)
            if styles is None:
                continue

            if isinstance(styles[0], str):
                styles = [styles]

            for i in range(len(styles)):
                if interaction == 'kspace':
                    coeff = [float(c) for c in styles[i][1:]]
                    style = styles[i][0]

                else:
                    coeff = self.get("%s_coeff" % interaction)
                    style = ' '.join([str(si) for si in styles[i]])

                styles_coeffs.append((style.strip(), coeff))

        return styles_coeffs


class LammpsParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/lammps', code_name='LAMMPS', code_homepage='https://lammps.sandia.gov/',
            domain='dft', mainfile_contents_re=r'^LAMMPS')

        self._metainfo_env = m_env
        self.log_parser = LogParser()
        self.traj_parser = TrajParser()
        self.data_parser = DataParser()

    def parse_thermodynamic_data(self):
        thermo_data = self.log_parser.get_thermodynamic_data()
        if not thermo_data:
            return
        n_evaluations = len(thermo_data['Step'])

        energy_keys_mapping = {
            'e_pair': 'pair', 'e_bond': 'bond', 'e_angle': 'angle', 'e_dihed': 'dihedral',
            'e_impro': 'improper', 'e_coul': 'coulomb', 'e_vdwl': 'van der Waals',
            'e_mol': 'molecular', 'e_long': 'kspace long range',
            'e_tail': 'van der Waals long range', 'kineng': 'kinetic', 'poteng': 'potential',
        }

        sec_run = self.archive.section_run[-1]

        sec_sccs = sec_run.section_single_configuration_calculation

        create_scc = True
        if sec_sccs:
            if len(sec_sccs) != n_evaluations:
                self.logger.warn(
                    '''Mismatch in number of calculations and number of property
                    evaluations!, will create new sections''',
                    data=dict(n_calculations=len(sec_sccs), n_evaluations=n_evaluations))

            else:
                create_scc = False

        for n in range(n_evaluations):
            if create_scc:
                sec_scc = sec_run.m_create(SingleConfigurationCalculation)
            else:
                sec_scc = sec_sccs[n]

            for key, val in thermo_data.items():
                key = key.lower()
                if key in energy_keys_mapping:
                    sec_energy = sec_scc.m_create(EnergyContribution)
                    sec_energy.energy_contibution_kind = energy_keys_mapping[key]
                    sec_energy.energy_contribution_value = val[n]

                elif key == 'toteng':
                    sec_scc.energy_method_current = val[n]

                elif key == 'press':
                    sec_scc.pressure = val[n]

                elif key == 'temp':
                    sec_scc.temperature = val[n]

                elif key == 'step':
                    sec_scc.time_step = int(val[n])

                elif key == 'cpu':
                    sec_scc.time_calculation = float(val[n])

                else:
                    if n == 0:
                        self.logger.warn(
                            'Unsupported property in thermodynamic data',
                            data=dict(property=key))

    def parse_sampling_method(self):
        sec_run = self.archive.section_run[-1]
        sec_sampling_method = sec_run.m_create(SamplingMethod)

        run_style = self.log_parser.get('run_style', ['verlet'])[0]
        run = self.log_parser.get('run', [0])[0]

        units = self.log_parser.get('units', ['lj'])[0]
        time_unit = get_unit(units, 'time')
        timestep = self.log_parser.get('timestep', [0], unit=time_unit)[0]
        sampling_method, ensemble_type = self.log_parser.get_sampling_method()

        sec_sampling_method.x_lammps_integrator_type = run_style
        sec_sampling_method.x_lammps_number_of_steps_requested = run
        sec_sampling_method.x_lammps_integrator_dt = timestep
        sec_sampling_method.sampling_method = sampling_method
        sec_sampling_method.ensemble_type = ensemble_type

        thermo_settings = self.log_parser.get_thermostat_settings()
        target_T = thermo_settings.get('target_T', None)
        if target_T is not None:
            sec_sampling_method.x_lammps_thermostat_target_temperature = target_T
        thermostat_tau = thermo_settings.get('thermostat_tau', None)
        if thermostat_tau is not None:
            sec_sampling_method.x_lammps_thermostat_tau = thermostat_tau
        target_P = thermo_settings.get('target_P', None)
        if target_P is not None:
            sec_sampling_method.x_lammps_barostat_target_pressure = target_P
        barostat_tau = thermo_settings.get('barostat_P', None)
        if barostat_tau is not None:
            sec_sampling_method.x_lammps_barostat_tau = barostat_tau
        langevin_gamma = thermo_settings.get('langevin_gamma', None)
        if langevin_gamma is not None:
            sec_sampling_method.x_lammps_langevin_gamma = langevin_gamma

    def parse_system(self):
        sec_run = self.archive.section_run[-1]

        pbc_cell = self.traj_parser.get('pbc_cell', [])
        n_atoms = self.traj_parser.get('n_atoms', [])

        if len(n_atoms) == 0:
            return

        create_scc = True
        sec_sccs = sec_run.section_single_configuration_calculation
        if sec_sccs:
            if len(sec_sccs) != len(pbc_cell):
                self.logger.warn(
                    '''Mismatch in number of calculations and number of structures,
                    will create new sections''',
                    data=dict(n_calculations=len(sec_sccs), n_structures=len(pbc_cell)))

            else:
                create_scc = False

        distance_unit = self.log_parser.units.get('distance', None)
        self.traj_parser._units = self.log_parser.units

        for i in range(len(pbc_cell)):
            sec_system = sec_run.m_create(System)
            sec_system.number_of_atoms = n_atoms[i]
            sec_system.configuration_periodic_dimensions = pbc_cell[i][0]
            sec_system.simulation_cell = pbc_cell[i][1] * distance_unit
            sec_system.atom_positions = self.traj_parser.get_positions(i)
            atom_labels = self.traj_parser.get_atom_labels(i)
            if atom_labels is None:
                atom_labels = ['X'] * n_atoms[i]
            sec_system.atom_labels = atom_labels

            velocities = self.traj_parser.get_velocities(i)
            if velocities is not None:
                sec_system.atom_velocities = velocities

            forces = self.traj_parser.get_forces(i)
            if forces is not None:
                if create_scc:
                    sec_scc = sec_run.m_create(SingleConfigurationCalculation)
                else:
                    sec_scc = sec_sccs[i]

                sec_scc.atom_forces = forces

    def parse_topology(self):
        sec_run = self.archive.section_run[-1]

        if self.traj_parser.mainfile is None or self.data_parser.mainfile is None:
            return

        masses = self.data_parser.get('Masses', None)

        self.traj_parser.masses = masses

        sec_topology = sec_run.m_create(section_topology)
        sec_topology.number_of_topology_atoms = self.data_parser.get('atoms', [None])[0]

        interactions = self.log_parser.get_interactions()
        if not interactions:
            interactions = self.data_parser.get_interactions()

        for interaction in interactions:
            if not interaction[0] or not interaction[1]:
                continue
            sec_interaction = sec_topology.m_create(section_interaction)
            sec_interaction.interaction_kind = str(interaction[0])
            sec_interaction.interaction_parameters = [list(a) for a in interaction[1]]

    def parse_input(self):
        sec_run = self.archive.section_run[-1]
        sec_input_output_files = sec_run.m_create(x_lammps_section_input_output_files)

        if self.data_parser.mainfile is not None:
            sec_input_output_files.x_lammps_inout_file_data = os.path.basename(
                self.data_parser.mainfile)

        if self.traj_parser.mainfile is not None:
            sec_input_output_files.x_lammps_inout_file_trajectory = os.path.basename(
                self.traj_parser.mainfile)

        sec_control_parameters = sec_run.m_create(x_lammps_section_control_parameters)
        keys = self.log_parser._commands
        for key in keys:
            val = self.log_parser.get(key, None)
            if val is None:
                continue
            val = val[0] if len(val) == 1 else val
            key = 'x_lammps_inout_control_%s' % key.replace('_', '').replace('/', '').lower()
            if hasattr(sec_control_parameters, key):
                if isinstance(val, list):
                    val = ' '.join([str(v) for v in val])
                setattr(sec_control_parameters, key, str(val))

    def init_parser(self):
        self.log_parser.mainfile = self.filepath
        self.log_parser.logger = self.logger
        self.traj_parser.logger = self.logger
        self.data_parser.logger = self.logger

    def reuse_parser(self, parser):
        self.log_parser.quantities = parser.log_parser.quantities
        self.traj_parser.quantities = parser.traj_parser.quantities
        self.data_parser.quantities = parser.data_parser.quantities

    def parse(self, filepath, archive, logger):
        self.filepath = filepath
        self.archive = archive
        self.logger = logger if logger is not None else logging

        self.init_parser()

        sec_run = self.archive.m_create(Run)

        # parse basic
        sec_run.program_name = 'LAMMPS'
        sec_run.program_version = self.log_parser.get('program_version', '')

        # parse method-related
        self.parse_sampling_method()

        # parse all data files associated with calculation, normally only one
        # specified in log file, otherwise will scan directory
        data_files = self.log_parser.get_data_files()
        for data_file in data_files:
            self.data_parser.mainfile = data_file
            self.parse_topology()

        # parse all trajectorty files associated with calculation, normally only one
        # specified in log file, otherwise will scan directory
        traj_files = self.log_parser.get_traj_files()
        for traj_file in traj_files:
            self.traj_parser.mainfile = traj_file
            self.parse_system()

        # parse thermodynamic data from log file
        self.parse_thermodynamic_data()

        # include input controls from log file
        self.parse_input()

        # create workflow
        if sec_run.section_sampling_method[0].sampling_method:
            sec_workflow = self.archive.m_create(Workflow)
            sec_workflow.workflow_type = sec_run.section_sampling_method[0].sampling_method

            if sec_workflow.workflow_type == 'molecular_dynamics':
                sec_md = sec_workflow.m_create(MolecularDynamics)

                sec_md.finished_normally = self.log_parser.get('finished') is not None
                sec_md.with_trajectory = self.traj_parser.with_trajectory()
                sec_md.with_thermodynamics = self.log_parser.get('thermo_data') is not None
