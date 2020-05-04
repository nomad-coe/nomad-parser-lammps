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

from LAMMPS.BaseClasses import AbsVariable
import logging

logger = logging.getLogger(__name__)


class Variable(AbsVariable):
    """Variables: This command assigns one or more strings to a variable name for evaluation later in the input script or during a simulation.
    - name (string) can only contain alphanumeric characters and underscores
    - value (“string”) can be simple text, it can contain other variables, or it can be a formula.

    equal-style variable VS internal-style variable

    The variables can NOT be re-defined in an input script. There are two exceptions to this rule:
      - First, variables of style string, getenv, internal, equal, vector, atom, and python ARE redefined each time the command is encountered.
      - Second, as described below, if a variable is iterated on to the end of its list of strings via the next command, it is removed from the list of active variables, and is thus available to be re-defined in a subsequent variable command. The delete style does the same thing.

    - the next command can be used to increment which string is assigned to the variable
    - equal: store a formula which when evaluated produces a single numeric value
    - vector: store a formula which produces a vector of such values
    - atom: store a formula which when evaluated produces one numeric value per atom
    - atomfile: can be used anywhere in an input script that atom-style variables are used; they get their per-atom values from a file rather than from a formula
    - python: can be hooked to Python functions using code you provide, so that the variable gets its value from the evaluation of the Python code
    - internal: are used by a few commands which set their value directly
    """
    names = {}

    def __init__(self, *args, **kwargs):
        pass

    def parse(self, name, style, *args):
        pass

    def delete(self, name):
        """delete = no args
        :param name: name of variable to define
        :return:
        """
        if name in self.names:
            del self.names[name]


    def index(self, name, args):
        """index
        :param name: name of variable to define
        :param args: one or more strings
        :return:
        """
        self.names[name] = args

    def loop(self, name, *args):
        """ loop args = N
              N = integer size of loop, loop from 1 to N inclusive
            loop args = N pad
              N = integer size of loop, loop from 1 to N inclusive
              pad = all values will be same length, e.g. 001, 002, ..., 100
            loop args = N1 N2
              N1,N2 = loop from N1 to N2 inclusive
            loop args = N1 N2 pad
              N1,N2 = loop from N1 to N2 inclusive
              pad = all values will be same length, e.g. 050, 051, ..., 100
        :param name: name of variable to define
        :param N: integer size of loop, loop from 1 to N inclusive
        """

        if len(args) == 1:
            self.names[name] = [x for x in range(1, N+1)]
        elif len(args) == 2:
            pass
        elif len(args) == 3:
            pass
        else:
            pass


    def world(self, *args, **kwargs):
        logger.warn("Variable style of {} is not supported yet.".format("word"))
        pass

    def universe(self, *args, **kwargs):
        pass

    def uloop(self, *args, **kwargs):
        pass

    def string(self, *args, **kwargs):
        pass

    def format(self, *args, **kwargs):
        pass

    def getenv(self, *args, **kwargs):
        pass

    def file(self, *args, **kwargs):
        pass

    def atomfile(self, *args, **kwargs):
        pass

    def python(self, *args, **kwargs):
        pass

    def equal(self, *args, **kwargs):
        pass

    def vector(self, *args, **kwargs):
        pass

    def atom(self, *args, **kwargs):
        pass





if __name__ == '__main__':
    pass
    v = Variable()

    str="a b c"
    v.index('a', str.split())

    v.delete('a')

    v.world('a')

    # v.loop('b', 3)


    logger.warn("Not supported yet.")
    print(v.names)
