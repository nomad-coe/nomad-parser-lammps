from LAMMPS.BaseClasses import Region


class Block(Region):
    name = 'block'

    def __init__(self, *args):
        pass


class Cone(Region):
    name = 'cone'

    def __init__(self, *args):
        pass


class Cylinder(Region):
    name = 'cylinder'

    def __init__(self, *args):
        pass


class Plane(Region):
    name = 'plane'

    def __init__(self, *args):
        pass


class Prism(Region):
    name = 'prism'

    def __init__(self, *args):
        pass


class Sphere(Region):
    name = 'sphere'

    def __init__(self, *args):
        pass


class Union(Region):
    name = 'union'

    def __init__(self, *args):
        pass


class Intersect(Region):
    name = 'intersect'

    def __init__(self, *args):
        pass


if __name__ == '__main__':
    pass
