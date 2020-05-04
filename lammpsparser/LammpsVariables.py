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

class LammpsVariables:

    def __init__(self):
        self.vars = {}

    def eval(self, exp):
        pass

    def parse(self, input):
        """parse
        :param input: "variable name style args"
        :return:
        """
        line = input.split(maxsplit=3)
        try:
            getattr(self, line[1])(line[0], line[2:])
        except AttributeError:
            print(line[1] + ": currently unsupported")

            # width self.read() as name,style,args:

    def delete(self, name, *args):
        """delete
        delete = no args
        """
        del self.vars[name]

    def index(self, name, *args):
        """index
        index args = one or more strings
        """
        self.vars[name] = args

    def loop(self, name, *args):
        """loop
        loop args = N
          N = integer size of loop, loop from 1 to N inclusive
        loop args = N pad
          N = integer size of loop, loop from 1 to N inclusive
          pad = all values will be same length, e.g. 001, 002, ..., 100
        loop args = N1 N2
          N1,N2 = loop from N1 to N2 inclusive
        loop args = N1 N2 pad
          N1,N2 = loop from N1 to N2 inclusive
          pad = all values will be same length, e.g. 050, 051, ..., 100
        """
        self.vars[name] = args

    def world(self, name, *args):
        """word
        world args = one string for each partition of processors
        """
        self.vars[name] = args

    def universe(self, name, *args):
        """universe
        universe args = one or more strings
        """
        self.vars[name] = args

    def uloop(self, name, *args):
        """uloop
        uloop args = N
          N = integer size of loop
        uloop args = N pad
          N = integer size of loop
          pad = all values will be same length, e.g. 001, 002, ..., 100
        """
        self.vars[name] = args

    def string(self, name, *args):
        """String vaiable
        string arg = one string
        """
        self.vars[name] = args

    def format(self, name, *args):
        """format
        format args = vname fstr
          vname = name of equal-style variable to evaluate
          fstr = C-style format string
        """
        self.vars[name] = args

    def getenv(self, name, *args):
        """getenv
        getenv arg = one string
        """
        self.vars[name] = args

    def file(self, name, *args):
        """file
        file arg = filename
        """
        self.vars[name] = args

    def atomfile(self, name, *args):
        """atomfile
        atomfile arg = filename
        """
        self.vars[name] = args

    def python(self, name, *args):
        """python
        python arg = function
        """
        self.vars[name] = args

    """equal
    equal or vector or atom args = one formula containing numbers, thermo keywords, math operations, group functions, atom values and vectors, compute/fix/variable references
    numbers = 0.0, 100, -5.4, 2.8e-4, etc
    constants = PI, version, on, off, true, false, yes, no
    thermo keywords = vol, ke, press, etc from thermo_style
    math operators = (), -x, x+y, x-y, x*y, x/y, x^y, x%y,
    x == y, x != y, x < y, x <= y, x > y, x >= y, x && y, x || y, !x
    math functions = sqrt(x), exp(x), ln(x), log(x), abs(x),
    sin(x), cos(x), tan(x), asin(x), acos(x), atan(x), atan2(y,x),
    random(x,y,z), normal(x,y,z), ceil(x), floor(x), round(x)
    ramp(x,y), stagger(x,y), logfreq(x,y,z), logfreq2(x,y,z),
    stride(x,y,z), stride2(x,y,z,a,b,c),
    vdisplace(x,y), swiggle(x,y,z), cwiggle(x,y,z)
    group functions = count(group), mass(group), charge(group),
    xcm(group,dim), vcm(group,dim), fcm(group,dim),
    bound(group,dir), gyration(group), ke(group),
    angmom(group,dim), torque(group,dim),
    inertia(group,dimdim), omega(group,dim)
    region functions = count(group,region), mass(group,region), charge(group,region),
    xcm(group,dim,region), vcm(group,dim,region), fcm(group,dim,region),
    bound(group,dir,region), gyration(group,region), ke(group,reigon),
    angmom(group,dim,region), torque(group,dim,region),
    inertia(group,dimdim,region), omega(group,dim,region)
    special functions = sum(x), min(x), max(x), ave(x), trap(x), slope(x), gmask(x), rmask(x), grmask(x,y), next(x)
    feature functions = is_active(category,feature,exact), is_defined(category,id,exact)
    atom value = id[i], mass[i], type[i], mol[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], fx[i], fy[i], fz[i], q[i]
    atom vector = id, mass, type, mol, x, y, z, vx, vy, vz, fx, fy, fz, q
    compute references = c_ID, c_ID[i], c_ID[i][j], C_ID, C_ID[i]
    fix references = f_ID, f_ID[i], f_ID[i][j], F_ID, F_ID[i]
    variable references = v_name, v_name[i]
    """


    def equal(self, name, *args):
        self.vars[name] = args

    def vector(self, name, *args):
        self.vars[name] = args

    def atom(self, name, *args):
        self.vars[name] = args







if __name__ == '__main__':

    vars = LammpsVariables()

    examples = [
        'variable x index run1 run2 run3 run4 run5 run6 run7 run8',
        'variable LoopVar loop $n',
        'variable beta equal temp/3.0',
        'variable b1 equal x[234]+0.5*vol',
        'variable b1 equal "x[234] + 0.5*vol"',
        'variable b equal xcm(mol1,x)/2.0',
        'variable b equal c_myTemp',
        'variable b atom x*y/vol',
        'variable foo string myfile',
        'variable myPy python increase',
        'variable f file values.txt',
        'variable temp world 300.0 310.0 320.0 ${Tfinal}',
        'variable x universe 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15',
        'variable x uloop 15 pad',
        'variable str format x %.6g',
        'variable x delete',
    ]

    for line in examples:
        exp = line[len('variable'):]
        print(exp)
        vars.parse(exp)

    pass


