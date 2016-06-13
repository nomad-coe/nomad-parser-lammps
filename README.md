# LAMMPS Parser

This is the parser for [LAMMPS](http://lammps.sandia.gov/).
It is part of the [NOMAD Laboratory](http://nomad-lab.eu).
The official version lives at

    git@gitlab.mpcdf.mpg.de:nomad-lab/parser-lammps.git

you can browse it at

    https://gitlab.mpcdf.mpg.de/nomad-lab/parser-lammps

It relies on having the nomad-meta-info and the python common repositories one level higher.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

then this will be in parsers/lammps.

# Running and Testing the Parser

## Usage
LAMMPS input/output files can be parsed with:

    python LAMMPSParserControl.py ../../test/examples/input.1_methyl_naphthalene

## Test Files
Example input/output files of LAMMPS can be found in the directory test/examples.

    - /methane/                 : NVT simulation on pure CH4 with different logging styles (see input file name)
    - /methane_MD_traj/         : NVT simulation on pure CH4 with different trajectory styles (see input file name)
    - /hexane_cyclohexane/      : NVT simulation on a binary (equimolar mixture)
    - /1_methyl_naphthalene/    : NVT simualtion on pure 1-methylnaphthalene (no trajectory)    


