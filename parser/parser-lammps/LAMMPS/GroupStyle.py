from LAMMPS.BaseClasses import Group

# Default
#
# All atoms belong to the “all” group.


class Delete(Group):
    name = 'delete'
    pass


class Region(Group):
    name = 'region'
    pass


class Type(Group):
    name = 'type'
    pass


class Id(Group):
    name = 'id'
    pass


class Molecule(Group):
    name = 'molecule'
    pass


class Variable(Group):
    name = 'variable'
    pass


class Include(Group):
    name = 'include'
    pass


class Subtract(Group):
    name = 'subtract'
    pass


class Union(Group):
    name = 'union'
    pass


class Intersect(Group):
    name = 'intersect'
    pass


class Dynamic(Group):
    name = 'dynamic'
    pass


class Static(Group):
    name = 'static'
    pass

