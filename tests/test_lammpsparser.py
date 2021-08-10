#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from lammpsparser import LammpsParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return LammpsParser()


def test_nvt(parser):
    archive = EntryArchive()
    parser.parse('tests/data/hexane_cyclohexane/log.hexane_cyclohexane_nvt', archive, None)

    sec_run = archive.run[0]
    assert sec_run.program.version == '14 May 2016'

    sec_workflow = archive.workflow[0]
    assert sec_workflow.type == 'molecular_dynamics'
    assert sec_workflow.molecular_dynamics.x_lammps_integrator_dt.magnitude == 2.5e-16
    assert sec_workflow.molecular_dynamics.x_lammps_thermostat_target_temperature.magnitude == 300.
    assert sec_workflow.molecular_dynamics.ensemble_type == 'nvt'

    sec_method = sec_run.method[0]
    assert len(sec_method.force_field.model[0].contributions) == 4
    assert sec_method.force_field.model[0].contributions[2].type == 'harmonic'
    assert sec_method.force_field.model[0].contributions[0].parameters[0][2] == 0.066

    sec_system = sec_run.system
    assert len(sec_system) == 201
    assert sec_system[5].atoms.lattice_vectors[1][1].magnitude == approx(2.24235e-09)
    assert False not in sec_system[0].atoms.periodic
    assert sec_system[80].atoms.labels[91:96] == ['H', 'H', 'H', 'C', 'C']

    sec_scc = sec_run.calculation
    assert len(sec_scc) == 201
    assert sec_scc[21].energy.current.value.magnitude == approx(8.86689197e-18)
    assert sec_scc[180].time_calculation.magnitude == 218.5357
    assert sec_scc[56].thermodynamics[0].pressure.magnitude == approx(-77642135.4975)
    assert sec_scc[103].thermodynamics[0].temperature.magnitude == 291.4591
    assert sec_scc[11].thermodynamics[0].time_step == 4400
    assert len(sec_scc[1].energy.contributions) == 9
    assert sec_scc[112].energy.contributions[8].kind == 'kspace long range'
    assert sec_scc[96].energy.contributions[2].value.magnitude == approx(1.19666271e-18)
    assert sec_scc[47].energy.contributions[4].value.magnitude == approx(1.42166035e-18)

    assert sec_run.x_lammps_section_control_parameters[0].x_lammps_inout_control_atomstyle == 'full'


def test_thermo_format(parser):
    archive = EntryArchive()
    parser.parse('tests/data/1_methyl_naphthalene/log.1_methyl_naphthalene', archive, None)

    sec_sccs = archive.run[0].calculation
    assert len(sec_sccs) == 301
    assert sec_sccs[98].energy.total.value.magnitude == approx(1.45322428e-17)

    assert len(archive.run[0].system) == 4


def test_traj_xyz(parser):
    archive = EntryArchive()
    parser.parse('tests/data/methane_xyz/log.methane_nvt_traj_xyz_thermo_style_custom', archive, None)

    sec_systems = archive.run[0].system
    assert len(sec_systems) == 201
    assert sec_systems[13].atoms.positions[7][0].magnitude == approx(-8.00436e-10)


def test_traj_dcd(parser):
    archive = EntryArchive()
    parser.parse('tests/data/methane_dcd/log.methane_nvt_traj_dcd_thermo_style_custom', archive, None)

    assert len(archive.run[0].calculation) == 201
    sec_systems = archive.run[0].system
    assert np.shape(sec_systems[56].atoms.positions) == (320, 3)
    assert len(sec_systems[107].atoms.labels) == 320
