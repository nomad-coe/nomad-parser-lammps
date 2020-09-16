import numpy as np
import os
import logging
from ase import data as asedata

from .LammpsCommon import converter

from nomad.parsing.text_parser import Quantity, UnstructuredTextFileParser
from nomad.datamodel.metainfo.public import section_run, section_sampling_method,\
    section_system, section_single_configuration_calculation, section_energy_contribution,\
    Workflow, MolecularDynamics
from nomad.datamodel.metainfo.common import section_topology, section_interaction


class Parser:
    def __init__(self, mainfile, logger):
        self._mainfile = mainfile
        self._results = None
        self._quantities = None
        self.logger = logger

    def __getitem__(self, key):
        val = self.results.get(key, None)
        return val

    def get(self, key, default=None):
        val = self[key]
        if not val:
            return default
        return val

    @property
    def maindir(self):
        return os.path.dirname(self._mainfile)

    @property
    def mainfile(self):
        return self._mainfile

    @mainfile.setter
    def mainfile(self, value):
        self._mainfile = value
        self._results = None

    @property
    def results(self):
        if self._results is None:
            self.parse()
        if self._results is None:
            self._results = {}

        return self._results

    def parse(self):
        pass


class DataParser(Parser):
    def __init__(self, mainfile, logger):
        super().__init__(mainfile, logger)
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

    def parse(self):
        if self.mainfile is None:
            return

        if not os.path.isfile(self.mainfile):
            return

        if self._quantities is None:
            self._quantities = []
            for header in self._headers:
                self._quantities.append(
                    Quantity(header, r'\s*([\+\-eE\d\. ]+)\s*%s\s*\n' % header, comment='#'))

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
                        str_operation=get_section_value))

        parser = UnstructuredTextFileParser(self.mainfile, self._quantities)

        self._results = {key: parser[key] for key in parser.keys() if parser[key]}

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


class TrajParser(Parser):
    def __init__(self, mainfile, logger):
        super().__init__(mainfile, logger)
        self._masses = None
        self._reference_masses = dict(
            masses=np.array(asedata.atomic_masses), symbols=asedata.chemical_symbols)
        self._chemical_symbols = None

    def parse(self):
        if self.mainfile is None:
            return

        if not os.path.isfile(self.mainfile):
            return

        if self._quantities is None:
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
                    'time_step', r'\s*ITEM:\s*TIMESTEP\s*\n\s*(\d+)\s*\n', comment='#'),
                Quantity(
                    'n_atoms', r'\s*ITEM:\s*NUMBER OF ATOMS\s*\n\s*(\d+)\s*\n', comment='#'),
                Quantity(
                    'pbc_cell', r'\s*ITEM: BOX BOUNDS\s*([\s\w]+)([\+\-\d\.eE\s]+)\n',
                    str_operation=get_pbc_cell, comment='#'),
                Quantity(
                    'atoms_info', r's*ITEM:\s*ATOMS\s*([ \w]+\n)*?([\+\-eE\d\.\n ]+)\n*I*',
                    str_operation=get_atoms_info, comment='#')
            ]

        parser = UnstructuredTextFileParser(self.mainfile, self._quantities)

        self._results = {key: parser[key] for key in parser.keys() if parser[key]}

    def with_trajectory(self):
        if self._results is None:
            return False

        return self._results.get('atoms_info', None) is not None

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
        atoms_info = self.results.get('atoms_info')[idx]
        atoms_type = atoms_info.get('type', None)
        if atoms_type is None:
            return

        if self._chemical_symbols is None:
            return

        atom_labels = [self._chemical_symbols[atype] for atype in atoms_type]

        return atom_labels

    def get_positions(self, idx):
        atoms_info = self.results.get('atoms_info')[idx]

        cell = self.results.get('pbc_cell')[idx][1]

        if 'xs' in atoms_info and 'ys' in atoms_info and 'zs' in atoms_info:
            positions = np.array([atoms_info['xs'], atoms_info['ys'], atoms_info['zs']]).T

            positions = positions * np.linalg.norm(cell, axis=1) + np.amin(cell, axis=1)

        else:
            positions = np.array([atoms_info['x'], atoms_info['y'], atoms_info['z']]).T

            if 'ix' in atoms_info and 'iy' in atoms_info and 'iz' in atoms_info:
                positions_img = np.array([
                    atoms_info['ix'], atoms_info['iy'], atoms_info['iz']]).T

                positions += positions_img * np.linalg.norm(cell, axis=1)

        return positions

    def get_velocities(self, idx):
        atoms_info = self.results.get('atoms_info')[idx]

        if 'vx' not in atoms_info or 'vy' not in atoms_info or 'vz' not in atoms_info:
            return

        return np.array([atoms_info['vx'], atoms_info['vy'], atoms_info['vz']]).T

    def get_forces(self, idx):
        atoms_info = self.results.get('atoms_info')[idx]

        if 'fx' not in atoms_info or 'fy' not in atoms_info or 'fz' not in atoms_info:
            return

        return np.array([atoms_info['fx'], atoms_info['fy'], atoms_info['fz']]).T


class LogParser(Parser):
    def __init__(self, mainfile, logger):
        super().__init__(mainfile, logger)
        self.logger = logger
        self._thermo_data = None
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

    def get_thermodynamic_data(self):
        def str_operation(val):
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

        pattern = r'\s*\-*(\s*Step\s*[\-\s\w\.\=\(\)]*[ \-\.\d\n]+)Loop'

        parser = UnstructuredTextFileParser(
            self.mainfile, [Quantity('thermo_data', pattern, str_operation=str_operation)])

        thermo_data = parser['thermo_data']

        self._thermo_data = list(thermo_data)[0] if thermo_data is not None else thermo_data

        return self._thermo_data

    def get_traj_files(self):
        dump = self.results.get('dump', None)
        if dump is None:
            self.logger.warn(
                'Trajectory not specified in %s, will scan directory' % self.maindir)
            traj_files = os.listdir(self.maindir)
            traj_files = [f for f in traj_files if f.endswith('trj')]

        else:
            traj_files = []
            if type(dump[0]) in [str, int]:
                dump = [dump]
            traj_files = [d[4] for d in dump]

        return [os.path.join(self.maindir, f) for f in traj_files]

    def get_data_files(self):
        read_data = self.results.get('read_data', None)
        if read_data is None:
            self.logger.warn(
                'Data file not specified in %s, will scan directory' % self.maindir)
            data_files = os.listdir(self.maindir)
            data_files = [f for f in data_files if f.endswith('data') or f.startswith('data')]

        else:
            data_files = read_data

        return [os.path.join(self.maindir, f) for f in data_files]

    def get_pbc(self):
        pbc = self.results.get('boundary', ['p', 'p', 'p'])
        return [v == 'p' for v in pbc]

    def get_program_version(self):
        version = self.results.get('program_version', [''])[0]
        return ' '.join(version)

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

        res = dict()
        if fix_style.lower() == 'nvt':
            try:
                res['target_T'] = float(fix[5])
                res['thermostat_tau'] = float(fix[6])
            except Exception:
                pass

        elif fix_style.lower() == 'npt':
            try:
                res['target_T'] = float(fix[5])
                res['thermostat_tau'] = float(fix[6])
                res['target_P'] = float(fix[9])
                res['barostat_tau'] = float(fix[10])
            except Exception:
                pass

        elif fix_style.lower() == 'nph':
            try:
                res['target_P'] = float(fix[5])
                res['barostat_tau'] = float(fix[6])
            except Exception:
                pass

        elif fix_style.lower() == 'langevin':
            try:
                res['target_T'] = float(fix[4])
                res['langevin_gamma'] = float(fix[5])
            except Exception:
                pass

        else:
            self.logger.warn('Fix style %s not supported' % fix_style)

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

    def finished_normally(self):
        return self.get('finished', None) is not None

    def with_thermodynamics(self):
        return self._thermo_data is not None

    def parse(self):
        if self.mainfile is None:
            return

        if not os.path.isfile(self.mainfile):
            return

        if self._quantities is None:
            def str_op(val):
                val = val.split('#')[0]
                val = val.replace('&\n', ' ').split()
                val = val if len(val) > 1 else val[0]
                return val

            self._quantities = [
                Quantity(
                    name, r'\n\s*%s\s+([\w\. \/\#\-]+)(\&\n[\w\. \/\#\-]*)*' % name,
                    str_operation=str_op, comment='#') for name in self._commands]

            self._quantities.append(
                Quantity('program_version', r'\s*LAMMPS\s*\(([\w ]+)\)\n', dtype=str)
            )

            self._quantities.append(
                Quantity('finished', r'\s*Dangerous builds\s*=\s*(\d+)')
            )

        parser = UnstructuredTextFileParser(self.mainfile, self._quantities)

        self._results = {
            key: parser[key] for key in parser.keys() if parser[key]}


class LammpsOutput:
    def __init__(self, filepath, archive, logger=None):
        self.filepath = filepath
        self.archive = archive
        self.logger = logger if logger is not None else logging
        self.log_parser = LogParser(filepath, self.logger)
        self.traj_parser = TrajParser(None, self.logger)
        self.data_parser = DataParser(None, self.logger)
        self._converter = None

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
                    '''Mismatch in number of calculations %d and number of property
                    evaluations %d!, will create new sections''' % (
                        len(sec_sccs), n_evaluations))

            else:
                create_scc = False

        for n in range(n_evaluations):
            if create_scc:
                sec_scc = sec_run.m_create(section_single_configuration_calculation)
            else:
                sec_scc = sec_sccs[n]

            for key, val in thermo_data.items():
                key = key.lower()
                if key in energy_keys_mapping:
                    sec_energy = sec_scc.m_create(section_energy_contribution)
                    sec_energy.energy_contibution_kind = energy_keys_mapping[key]
                    sec_energy.energy_contribution_value = self._converter.Energy(val[n])

                elif key == 'toteng':
                    sec_scc.energy_method_current = self._converter.Energy(val[n])

                elif key == 'press':
                    sec_scc.pressure = self._converter.Press(val[n])

                elif key == 'temp':
                    sec_scc.temperature = self._converter.Temp(val[n])

                elif key == 'step':
                    sec_scc.time_step = int(val[n])

                elif key == 'cpu':
                    sec_scc.time_calculation = float(val[n])

                else:
                    if n == 0:
                        self.logger.warn(
                            'Unsupported property %s in thermodynamic data' % key)

    def parse_sampling_method(self):
        sec_run = self.archive.section_run[-1]
        sec_sampling_method = sec_run.m_create(section_sampling_method)

        run_style = self.log_parser.get('run_style', ['verlet'])[0]
        run = self.log_parser.get('run', [0])[0]
        timestep = self._converter.Time(self.log_parser.get('timestep', [0])[0])
        sampling_method, ensemble_type = self.log_parser.get_sampling_method()

        sec_sampling_method.x_lammps_integrator_type = run_style
        sec_sampling_method.x_lammps_number_of_steps_requested = run
        sec_sampling_method.x_lammps_integrator_dt = timestep
        sec_sampling_method.sampling_method = sampling_method
        sec_sampling_method.ensemble_type = ensemble_type

        thermo_settings = self.log_parser.get_thermostat_settings()
        target_T = thermo_settings.get('target_T', None)
        if target_T is not None:
            target_T = self._converter.Temp(target_T)
            sec_sampling_method.x_lammps_thermostat_target_temperature = target_T
        thermostat_tau = thermo_settings.get('thermostat_tau', None)
        if thermostat_tau is not None:
            thermostat_tau = self._converter.Time(thermostat_tau)
            sec_sampling_method.x_lammps_thermostat_tau = thermostat_tau
        target_P = thermo_settings.get('target_P', None)
        if target_P is not None:
            target_P = self._converter.Press(target_P)
            sec_sampling_method.x_lammps_barostat_target_pressure = target_P
        barostat_tau = thermo_settings.get('barostat_P', None)
        if barostat_tau is not None:
            barostat_tau = self._converter.Time(barostat_tau)
            sec_sampling_method.x_lammps_barostat_tau = barostat_tau
        langevin_gamma = thermo_settings.get('langevin_gamma', None)
        if langevin_gamma is not None:
            langevin_gamma = self._converter.Time(langevin_gamma)
            sec_sampling_method.x_lammps_langevin_gamma = langevin_gamma

    def parse_system(self, trajfile):
        sec_run = self.archive.section_run[-1]

        self.traj_parser.mainfile = trajfile

        pbc_cell = self.traj_parser.get('pbc_cell', [])
        n_atoms = self.traj_parser.get('n_atoms', [])

        create_scc = True
        sec_sccs = sec_run.section_single_configuration_calculation
        if sec_sccs:
            if len(sec_sccs) != len(pbc_cell):
                self.logger.warn(
                    '''Mismatch in number of calculations %d and number of property
                    evaluations %d!, will create new sections''' % (
                        len(sec_sccs), len(pbc_cell)))

            else:
                create_scc = False

        for i in range(len(pbc_cell)):
            sec_system = sec_run.m_create(section_system)
            sec_system.number_of_atoms = n_atoms[i]
            sec_system.configuration_periodic_dimensions = pbc_cell[i][0]
            sec_system.simulation_cell = self._converter.Distance(pbc_cell[i][1])
            sec_system.atom_positions = self.traj_parser.get_positions(i)
            atom_labels = self.traj_parser.get_atom_labels(i)
            if atom_labels is None:
                atom_labels = ['X'] * n_atoms[i]
            sec_system.atom_labels = atom_labels

            velocities = self.traj_parser.get_velocities(i)
            if velocities is not None:
                sec_system.atom_velocities = self._converter.Velocity(velocities)

            forces = self.traj_parser.get_forces(i)
            if forces is not None:
                if create_scc:
                    sec_scc = sec_run.m_create(section_single_configuration_calculation)
                else:
                    sec_scc = sec_sccs[i]

                sec_scc.atom_forces = self._converter.Force(forces)

    def parse_topology(self, datafile):
        sec_run = self.archive.section_run[-1]

        self.data_parser.mainfile = datafile

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

    def parse(self):
        sec_run = self.archive.m_create(section_run)
        # TODO Warning! lj units are not supported in LammpsCommon!
        self._converter = converter(self.log_parser.get('units', ['lj'])[0])

        # parse basic
        sec_run.program_name = 'LAMMPS'
        sec_run.program_version = self.log_parser.get_program_version()

        # parse method-related
        self.parse_sampling_method()

        # parse all data files associated with calculation, normally only one
        # specified in log file, otherwise will scan directory
        data_files = self.log_parser.get_data_files()
        for data_file in data_files:
            self.parse_topology(data_file)

        # parse all trajectorty files associated with calculation, normally only one
        # specified in log file, otherwise will scan directory
        traj_files = self.log_parser.get_traj_files()
        for traj_file in traj_files:
            self.parse_system(traj_file)

        # parse thermodynamic data from log file
        self.parse_thermodynamic_data()

        # create workflow
        if sec_run.section_sampling_method[0].sampling_method:
            sec_workflow = self.archive.m_create(Workflow)
            sec_workflow.workflow_type = sec_run.section_sampling_method[0].sampling_method

            if sec_workflow.workflow_type == 'molecular_dynamics':
                sec_md = sec_workflow.m_create(MolecularDynamics)

                sec_md.finished_normally = self.log_parser.finished_normally()
                sec_md.with_trajectory = self.traj_parser.with_trajectory()
                sec_md.with_thermodynamics = self.log_parser.with_thermodynamics()
