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

from LAMMPS.BaseClasses import abstract_cmd, abstract_cmd_arg


def read_data(filename, *args):
    """Wrapper for read_data command"""


    class add(abstract_cmd_arg):
        """add arg = append or Nstart or merge
                append = add new atoms with IDs appended to current IDs
                Nstart = add new atoms with IDs starting with Nstart
                merge = add new atoms with their IDs unchanged"""

        def __init__(self, arg):
            self.arg = arg

    class offset(abstract_cmd_arg):
        """offset args = toff boff aoff doff ioff
                toff = offset to add to atom types
                boff = offset to add to bond types
                aoff = offset to add to angle types
                doff = offset to add to dihedral types
                ioff = offset to add to improper types"""

        def __init__(self, toff, boff, aoff, doff, ioff):
            self.toff = toff
            self.boff = boff
            self.aoff = aoff
            self.doff = doff
            self.ioff = ioff

    class shift(abstract_cmd_arg):
        """shift args = Sx Sy Sz
                Sx,Sy,Sz = distance to shift atoms when adding to system (distance units)"""

        def __init__(self, Sx, Sy, Sz):
            self.Sx = Sx
            self.Sy = Sy
            self.Sz = Sz

    class group(abstract_cmd_arg):
        """group args = groupID
                groupID = add atoms in data file to this group"""
        def __init__(self, group_id):
            self.group_id = group_id

    class nocoeff(abstract_cmd_arg):
        """nocoeff = ignore force field parameters"""

    class fix(abstract_cmd_arg):
        """fix args = fix-ID header-string section-string
                fix-ID = ID of fix to process header lines and sections of data file
                header-string = header lines containing this string will be passed to fix
                section-string = section names with this string will be passed to fix"""

        def __init__(self, id, header, section):
            self.id = id
            self.header = header
            self.section = section


    class extra_atom_types(abstract_cmd_arg):
        """extra/atom/types arg = # of extra atom types"""
        pass
    
    class extra_bond_types(abstract_cmd_arg):
        """extra/bond/types arg = # of extra bond types"""
        pass
    
    class extra_angle_types(abstract_cmd_arg):
        """extra/angle/types arg = # of extra angle types"""
        pass
    
    class extra_dihedral_types(abstract_cmd_arg):
        """extra/dihedral/types arg = # of extra dihedral types"""
        pass
    
    class extra_improper_types(abstract_cmd_arg):
        """extra/improper/types arg = # of extra improper types"""
        pass
    

if __name__ == '__main__':

    read_data('data.lj')
    read_data('../run7/data.polymer.gz')
    read_data('data.protein', 'fix', 'mycmap', 'crossterm' 'CMAP')
    read_data('data.water', 'add', 'append', 'offset', 3, 1, 1, 1, 1, 'shift', 0.0, 0.0, 50.0)
    # data = read_data('data.water')
    # data.add('append')
    # data.offset(3, 1, 1, 1, 1)
    # data.shift( 0.0, 0.0, 50.0)

    read_data('data.water', 'add', 'merge', 1, 'group', 'solvent')
