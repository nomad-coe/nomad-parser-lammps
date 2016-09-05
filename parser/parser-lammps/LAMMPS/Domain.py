from LAMMPS.BaseClasses import Domain as AbsDomain


class Domain(AbsDomain):

    def __init__(self, dimension=3):
        self.dimension = dimension
        self.periodicity = None
        self.boundary = None

        self.minxlo = 0.0
        self.minxhi = 0.0
        self.minylo = 0.0
        self.minyhi = 0.0
        self.minzlo = 0.0
        self.minzhi = 0.0

if __name__ == '__main__':
    pass
