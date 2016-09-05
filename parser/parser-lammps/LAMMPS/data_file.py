import logging
from collections import UserDict, OrderedDict
from abc import ABCMeta
import numpy as np

# logging.basicConfig()
logger = logging.getLogger(__name__)


class AtomType(OrderedDict):
    pass


class Atoms:
    def __init__(self):
        self.data = []
        self.atom_id = {}
        self.molecule_id = {}

    def add(self, id, values):
        self.data.append(values)
        self.atom_id[id] = values


    # def __hash__(self):
    #     return hash(self.atom_id)


class AtomStyle(metaclass=ABCMeta):
    pass


class AtomStyleFull(AtomStyle):
    name = 'full'

    # atom-ID molecule-ID atom-type q x y z
    def __init__(self, tom_id, molecule_id, atom_type, q, x, y, z):
        super(AtomStyleFull, self).__init__()

        self.tom_id = tom_id
        self.molecule_id = molecule_id
        self.atom_type = atom_type
        self.q = q
        self.x = x
        self.y = y
        self.z = z


    # @property
    # def x(self):
    #     return self.x


"""TODO
filereader
data numpy? dict? class?

"""

class DataFile(object):
    def __init__(self, filename):
        self.filereader = filename


    @property
    def filereader(self):
        return self._filereader

    @filereader.setter
    def filereader(self, filename):
        self._filereader = filename



if __name__ == '__main__':

    logger.info("Start Datafile module")
    print("hello")

    data = DataFile("data.methene")

    a = Atoms()
    a.add(1,[1,2,3])
    a.add(2,[11,22,33])
    print(a.data)
    print(a.atom_id)

    a.atom_id[1][2]=11111
    print(a.data)
    print(a.atom_id)


    atom_type = AtomType()
    atom_type[2] = {'mass': 1}
    atom_type[1] = {'mass': 11}

    print(atom_type)
























    pass
