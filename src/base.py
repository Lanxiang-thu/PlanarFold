import os
import sys
import operator
import numpy as np
import dsflyext as ext
from scipy.io.netcdf import netcdf_file

#===============================================================================
class AnsiColor:
    # normal; k:black;  lr:light red
    # to be added
    n = '\x1b[0;0m'
    k = '\x1b[0;30m'
    r = '\x1b[0;31m'
    g = '\x1b[0;32m'
    y = '\x1b[0;33m'
    b = '\x1b[0;34m'
    p = '\x1b[0;35m'
    c = '\x1b[0;36m'
    lk = '\x1b[1;30m'
    lr = '\x1b[1;31m'
    lg = '\x1b[1;32m'
    ly = '\x1b[1;33m'
    lb = '\x1b[1;34m'
    lp = '\x1b[1;35m'
    lc = '\x1b[1;36m'
    # background color
    _k = '\x1b[0;40m'
    _r = '\x1b[0;41m'
    _g = '\x1b[0;42m'
    _y = '\x1b[0;43m'
    _b = '\x1b[0;44m'
    _p = '\x1b[0;45m'
    _c = '\x1b[0;46m'


cl = AnsiColor()


#===============================================================================
def pager(texts):
    try:
        fp = os.popen('less -R','w')
        fp.writelines(texts)
        fp.close()
    except IOError:
        pass

#===============================================================================
# divide a list of text into blocks
# according to user-defined getkey()
# return a list containing subgroups
# apply for AAABBBBCC... key sequence
def divide(lines, getkey):
    keys = [getkey(a) for a in lines]
    mydiff = [b!=a for a,b in zip(keys[:-1],keys[1:])]
    idx = [i+1 for i,a in enumerate(mydiff) if a==1]
    start = [0]+idx
    end = idx+[len(lines)]
    return [lines[i:j] for i,j in zip(start,end)]

#===============================================================================
# divide a list of text into blocks
# according to user-definded isborder()
# return a list containing subgroups
# include: can be 'header', 'tailer' or None
# include='header': keep border as 1st line of each block; if 1st line of text
#                   is not border, then no border in this block; if last line
#                   of text is border, then this border is discarded
# include='tailer': keep border as the last line of each block; if the last line
#                   of text is not border, then no border in this block; if 1st
#                   line of text is border, then this border is discarded
# include = None:   discard border in each block
def partition(lines, isborder, include=None):
    if not lines:
        return []
    idx = [i for i,a in enumerate(lines) if isborder(a)]
    if not idx:
        return [lines]

    if include == 'header':
        if idx[0] != 0:
            idx = [0]+idx
        if idx[-1] != len(lines)-1:
            idx = idx+[len(lines)]
        start = idx[:-1]
        end = idx[1:]
        return [lines[i:j] for i,j in zip(start,end)]
    elif include == 'tailer':
        if idx[0] != 0:
            idx = [-1]+idx
        if idx[-1] != len(lines)-1:
            idx = idx+[len(lines)-1]
        start = idx[:-1]
        end = idx[1:]
        return [lines[i+1:j+1] for i,j in zip(start,end)]
    else:  # by default
        if idx[0] != 0:
            idx = [-1]+idx
        if idx[-1] != len(lines)-1:
            idx = idx+[len(lines)]
        start = idx[:-1]
        end = idx[1:]
        return [lines[i+1:j] for i,j in zip(start,end)]

#===============================================================================
def findcommon(a, b):
    if not isinstance(a, list):
        a = list(a)
    if not isinstance(b, list):
        b = list(b)
    c = list(set(a) & set(b))
    c.sort()
    idx1 = [a.index(item) for item in c]
    idx2 = [b.index(item) for item in c]
    return idx1,idx2

#===============================================================================
def range2list(mystr):
    fds = [map(int, fd.split('-')) for fd in mystr.split(',')]
    mylist = reduce(operator.add, [range(fd[0], fd[-1]+1) for fd in fds])
    return mylist

#===============================================================================
def alignstr(strs):
    str0 = strs[0]
    print('%s'%str0)
    for str1 in strs[1:]:
        cs = [(c1 if c1!=c0 else '-') for c0,c1 in zip(str0, str1)]
        cs = ''.join(cs)
        print('%s'%cs)


#===============================================================================
def dot2helix(filename):
    lines = open(filename).readlines()
    mark0 = lines[1].strip()
    #print("Restaint:\n", mark0)
    nRes = len(mark0)
    lbrk = 0
    rbrk = 0
    #--------------- check whether the '(' and ')' match----------------------
    for i in range(nRes):
        if mark0[i]=='(':
            lbrk += 1
        elif mark0[i]==')':
            rbrk += 1
    if lbrk!=rbrk:
        print("'('[{0:d}] not equil to ')'[{1:d})]".format(lbrk, rbrk))
        sys.exit(0)
    #else:
    #    print("RNA has {0:d} basepairs".format(lbrk))
    #------------------------------------------------------
    mark = '.'
    for i in range(nRes):
        if mark0[i]=='[' or mark0[i]==']':
            mark += '.'
        elif mark0[i]=='(' or mark0[i]==')' or mark0[i]=='.':
            mark += mark0[i]
        else:
            print('this is neither a dot nor a bracket:', mark0[i])
            sys.exit(0)
    #print("Restaint trimmed:\n", mark[1:])
    #------------------------------------------------------
    rBplist = []
    left = []
    slope = '+'
    stem = []
    notLp = []
    for i in range(nRes):
        if mark[i+1]=='(':
            slope += '+'
            if mark[i]=='(' or mark[i+2]=='(':
                left.append(i)
                notLp.append(1)
            else:   # this is a lone-pair
                notLp.append(0)
        elif mark[i+1]==')':
            slope += '-'
            if notLp[-1]:   # Exclude lone-pair
                stem.append(left[-1])
                stem.append(i)
                left.pop()
            notLp.pop()
        else:
            slope += slope[-1]
        if slope[i]== '+' and slope[i+1]=='-':
            rBplist.append(stem)
        elif slope[i]== '-' and slope[i+1]=='+':
            stem = []
    #print("slope:\n", slope[1:])
    i = 0
    nrBp = []
    maxHelix = 0
    for c in rBplist:
        length = int(len(c)/2)
        if (length>maxHelix):
            maxHelix = length
        nrBp.append(length)
        i += 1
        #print('chain-{0:d}({1:d}bp):'.format(i, length),  c)
    nchain = i
    for i in range(len(rBplist)):
        rBplist[i] += [0] * (maxHelix*2 - len(rBplist[i]))
    neorBplist = []
    for c in rBplist:
        for cc in c:
            neorBplist.append(cc)
    return [nchain, maxHelix, nrBp, neorBplist]


#===============================================================================
def dot2bps(filename, minloop):
    lines = open(filename).readlines()
    dots = lines[2].strip()
    nRes = len(dots)
    print('2D restraint:', dots)
    print('There are {} residues'.format(nRes))
    ddd = dots.replace('.', '').replace('(', '').replace(')', '')
    if len(ddd)>0:
        print("Error !!!")
        print("Restraint:", dots)
        print("Un-recognizable format:", ddd)
        print("Pseudo-knot is currently not allowed, please remove them")
        sys.exit(1)
    if len(dots.replace('(', '')) != len(dots.replace(')', '')):
        print("Error !!!")
        print("Unmatching '(' and ')' !")
        print("Please check your 2D structure restraint!")
        sys.exit(1)
    if minloop == 3:
        dots = dots.replace('(..)', '....')
    elif minloop == 4:
        dots = dots.replace('(..)', '....')
        dots = dots.replace('(...)', '.....')
    elif minloop == 5:
        dots = dots.replace('(..)', '....')
        dots = dots.replace('(...)', '.....')
        dots = dots.replace('(....)', '......')
    else:
        print("Error !!!")
        print("Currently, we don't suggest seting: min_loop>5")
        print("Please reconsider change 'min_loop' parameter to a smaller value!")
        sys.exit(1)
    print('Note: min_loop check will remove the unsatisfied basepairs!')
    print('2D restraint:', dots)
    bp5 = []
    bplist = -1 * np.ones(nRes, dtype='int32')
    for i in range(nRes):
        if dots[i] == '(':
            bp5.append(i)
        elif dots[i] == ')':
            #bps.append(bp5[-1])
            #bps.append(i)
            #
            #print(i, bp5[-1])
            bplist[bp5[-1]] = i
            bplist[i] = bp5[-1]
            #
            bp5.pop()
    #-------------------------------------------------
    #print('bplist:', bplist)
    bps = []
    for i in range(nRes):
        if bplist[i] > -1:
            bps.append(i)
            n = bplist[i]
            bps.append(n)
            bplist[i] = -1
            bplist[n] = -1
            #print(i, 'bplist:', bplist)
    nBP = int(len(bps) / 2)
    print('There are {} basepairs'.format(nBP))
    #print('bps:', bps)
    return [nBP, bps]

#====== Prmtop ================================================================
class Prmtop:
    pass

    #--------------------------------------------------------------------------
    def __init__(self, inp=None):
        if isinstance(inp, str):
            if inp[-7:] == '.prmtop':
                self.read(inp)
            else:
                print('Unsupported filename suffix for Prmtop:', inp)
                print('Standard filename suffix for Prmtop: xxxx.prmtop')
                sys.exit(0)
        else:
            print('wrong file name format', inp)
            sys.exit(0)

    #--------------------------------------------------------------------------
    def read(self, inp):
        #----------------------------------------------------------------------
        # READ_LIST: raed different section with a given format into a List
        #
        # Arguments:
        #   sect:    section to be read
        #   rdt:     read times of each raw
        #            (except the last raw becasut it might be incomplete)
        #   typ:     type of the data(a\I\E are optinal)
        #   rdl:     read length of each element
        #----------------------------------------------------------------------
        def read_list(sect, rdt, typ, rdl):
            # Read data from complete raws
            array = []
            # Skip the first two mark raws
            for nraw in range(2,len(sect)-1):
                i = 0
                j = rdl
                for t in range(0,rdt):
                    array.append(sect[nraw][i:j])
                    i = j
                    j = j+rdl
            # Calculate how many element and read thm in the last raw
            # -1 because of the '\n' at the end
            remain = int((len(sect[-1])-1)/rdl)
            i = 0
            j = rdl
            for t in range(0,remain):
                array.append(sect[-1][i:j])
                i = j
                j = j+rdl
            # Convert the data to desired type(str or int or float )
            if (typ == 'A' or typ == 'a'):
                return array
            elif (typ == 'I' or typ == 'i'):
                arr_int = list(map(int,array))
                return arr_int
            elif (typ == 'E' or typ == 'e'):
                arr_double = list(map(float,array))
                return arr_double

        #----------------------------------------------------------------------
        # READ_ARRAY: raed bonds/angles/dihedrals into numpy array
        #               (the default formats are 10I8 )
        #
        # Arguments:
        #   sect:   section to be read
        #   m:      row of numpy array (3/4/5 for bnd/angl/dihed,respectively)
        #   n:      column of numpy array (number of bnd/angl/dihed)
        #----------------------------------------------------------------------
        def read_array(sect, m, n):
            data = []
            arrays = ''
            # Skip the first two raws with % mark
            for nraw in range(2,len(sect)):
                # remove the '\n' at the end of each raw
                arrays = arrays+sect[nraw].rstrip('\n')
            a = 0
            b = 8 # 10I8 fortran format
            for i in range(m):
                for j in range(n):
                    data.append(int(arrays[a:b]))
                    a = a+8
                    b = b+8
            return np.array(data).reshape((n,m))

        #--------------------------------------------------------------------------
        lines = open(inp).readlines()

        # Info of first 4 raws are insignificant, pop them put
        for i in range(0,4):
            lines.pop(0)
        
        # Partition the different sections
        sections = partition(lines, lambda x:x[1:5]=='FLAG', include='header')
        
        # Read each section into corresponding 'self.variables'
        for section in sections:
            # First, the pointer information
            if section[0][6:14] == 'POINTERS':
                self.pointers = read_list(section, 10,'I', 8)
                nbonh = self.pointers[2]
                mbona = self.pointers[3]
                ntheth = self.pointers[4]
                mtheta = self.pointers[5]
                nphih = self.pointers[6]
                mphia = self.pointers[7]
                
            # Now, read atomname, charge, mass, etc.
            elif section[0][6:15] == 'ATOM_NAME':                   # string
                self.igraph = read_list(section, 20, 'a', 4)
            elif section[0][6:12] == 'CHARGE':
                self.charge = read_list(section, 5, 'E', 16)
            elif section[0][6:19] == 'ATOMIC_NUMBER':               # useless
                self.atnum = read_list(section, 10, 'I', 8)
            elif section[0][6:10] == 'MASS':
                self.amass = read_list(section, 5, 'E', 16)
            elif section[0][6:21] == 'ATOM_TYPE_INDEX':
                self.iac = read_list(section, 10,'I', 8)
            elif section[0][6:26] == 'NONBONDED_PARM_INDEX':
                self.ico = read_list(section, 10,'I', 8)
            elif section[0][6:25] == 'EXCLUDED_ATOMS_LIST':
                self.inb = read_list(section, 10,'I', 8)
                self.nNBres = len(self.inb)
                #print('EXCLUDED_ATOMS_LIST({0}): {1}'.format(self.nNBres, self.inb))

#====== Restart ===============================================================
class Restart:
    pass

    #--------------------------------------------------------------------------
    def __init__(self, inp=None):
        if isinstance(inp, str):
            if inp[-4:] == '.rst':
                self.read_rst(inp)
            elif inp[-7:] == '.inpcrd':
                self.read_inpcrd(inp)
            else:
                print('Unsupported filename suffix for Restart:', inp)
                print('Standard filename suffix for Restart: xxxx.rst or xxxx.inpcrd')
                sys.exit(0)
        else:
            print('wrong file name format', inp)
            sys.exit(0)
    
    #---------- Inpcrd --------------------------------------------------------
    def read_inpcrd(self, inp):
        lines = open(inp).readlines()
        
        # Skip the title (usually as 'default_name')
        self.title = lines[0].rstrip()

        # Only if the second raw contains only natom but not time info
        self.natom = int(lines[1][0:6])
        self.box = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        nline = int((self.natom+1) / 2)

        crd = []
        data = ''
        # Read coordiante lines
        for i in range(2, 2+nline):
            data += lines[i].rstrip('\n')
        # print('coordinates:',data)
        # Pass crd string to crd float
        k = 0
        for i in range(self.natom * 3):
            crd.append(float(data[k:k+12]))
            k += 12

        #print('len(lines)=', len(lines), 'self.natom=', self.natom)
        if len(lines) > self.natom:   # which means this inpcrd has no vel info
            vel = []
            data = ''
            # Read velocity lines
            for i in range(2+nline, 2+2*nline):
                data += lines[i].rstrip('\n')
            # Pass vel string to vel float
            k = 0
            for i in range(self.natom * 3):
                vel.append(float(data[k:k+12]))
                k += 12

        # Read the last raw (box_lengths and box_angles)
        #self.box[0] = float(lines[-1][0:12])
        #self.box[1] = float(lines[-1][12:24])
        #self.box[2] = float(lines[-1][24:36])
        #self.box[3] = float(lines[-1][36:48])
        #self.box[4] = float(lines[-1][48:60])
        #self.box[5] = float(lines[-1][60:72])
        # Read the last raw
        self.crd = np.array(crd, dtype='float64').flatten()
        self.vel = np.array(vel, dtype='float64').flatten()
        # Initailize velocities
        

    #----------- Rst ----------------------------------------------------------
    def read_rst(self, inp):
        nc = netcdf_file(inp, 'r', mmap = False)
        #print('\n read from{0}'.format(inp))
        vars = list(nc.variables.keys())
        # use 'cell' to store cell_lengths[0:3], cell_angles[3:6], and time[6]
        self.box = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        for i in range(0, len(vars)):
            var = vars[i]
            data = nc.variables[var]

            if var == 'spatial':                # string
                self.spatial = data[:]
            #elif var == 'time':                 # not use yet
            #    self.time = data[:]

            # Crd and vel info
            elif var == 'coordinates':
                # store coordinates as numpy array
                natom = data.shape[0]
                self.crd = np.array(data[:], dtype='float64').flatten()
                self.crd0 = data[:]
            elif var == 'velocities':
                # store velocities as numpy array
                natom = data.shape[0]
                self.vel = np.array(data[:], dtype='float64').flatten()
                self.vel0 = data[:]

            # Cell info
            elif var == 'cell_spatial':         # string
                self.cell_spatial = data[:]
            elif var == 'cell_lengths':
                self.box[0:3] = data[:]
            elif var == 'cell_angles':
                self.box[3:6] = data[:]
            elif var == 'cell_angular':         # string
                self.cell_angular = data[:]
        
        #print('the dtype of crd read in python', self.crd.dtype, 'crd[0]=', self.crd[0])
        #print('the dtype of vel read in python', self.vel.dtype, 'vel[0]=', self.vel[0])
        nc.close()

#====== Controls ===============================================================
class Controls:
    pass
    #--------------------------------------------------------------------------
    def __init__(self, inp=None):
        if isinstance(inp, str):
            if inp[-3:] == '.in':
                self.read_controls(inp)
            else:
                print('Unsupported filename suffix for Controls:', inp)
                print('Standard filename suffix for Controls: xxxx.in')
                sys.exit(0)
        else:
            print('wrong file name format', inp)
            sys.exit(0)
    #---------- Controls --------------------------------------------------------
    def read_controls(self, inp):
        #--- assign the default values ---------------
        self.ntpr = 1000
        self.ntwx = 1000
        self.ntr = 0
        self.ismd = 0
        self.dumpfreq = 1000
        self.stretch_force = 1.0
        self.stretch_spring = 0.01
        self.dist_init = 10.0
        self.dist_final = 110.0
        self.restraint_wt = 0.1
        self.restraint_eq = 10.0
        self.rcut = 11.5
        self.nscm = 1000
        self.nstlim = 1000
        self.dt = 0.004
        self.ntt = 1
        self.temp0 = 283
        self.tempi = 283
        self.duration = 10000
        self.ntemp = 1
        self.tautp = 0.005
        self.iCircular = 0
        self.verlet_cut = 30.0
        self.nVerlet = 500
        self.dpr = 0
        self.iCollide = 0
        #---- Read and undate the values ------------
        lines = open(inp).readlines()
        for i in range(len(lines)):
            if lines[i].strip() == '&cntl':
                #print(lines[i])
                for j in range(i+1, len(lines)):
                    if lines[j].strip() == '/':
                        break
                    variable, value = lines[j].split('=')
                    if value[-2] != ',':
                        print('No blank space is allow after the value:', lines[j])
                        sys.exit(0)
                    variable = variable.strip()
                    value = value.strip()
                    if variable == 'ntpr':
                        #print('#ntpr=', int(value[:-1]))
                        self.ntpr = int(value[:-1])
                    elif variable == 'ntwx':
                        self.ntwx = int(value[:-1])
                    elif variable == 'ntr':
                        self.ntr = int(value[:-1])
                    elif variable == 'ismd':
                        self.ismd = int(value[:-1])
                    elif variable == 'dumpfreq':
                        self.dumpfreq = int(value[:-1])
                    elif variable == 'stretch_force':
                        self.stretch_force = float(value[:-1])
                    elif variable == 'stretch_spring':
                        self.stretch_spring = float(value[:-1])
                    elif variable == 'dist_init':
                        self.dist_init = float(value[:-1])
                    elif variable == 'dist_final':
                        self.dist_final = float(value[:-1])
                    elif variable == 'restraint_wt':
                        self.restraint_wt = float(value[:-1])
                    elif variable == 'restraint_eq':
                        self.restraint_eq = float(value[:-1])
                    elif variable == 'rcut':
                        self.rcut = float(value[:-1])
                    elif variable == 'nscm':
                        self.nscm = int(value[:-1])
                    elif variable == 'nstlim':
                        self.nstlim = int(value[:-1])
                    elif variable == 'dt':
                        self.dt = float(value[:-1])
                    elif variable == 'ntt':
                        self.ntt = int(value[:-1])
                    elif variable == 'ntemp':
                        self.ntemp = int(value[:-1])
                    elif variable == 'temp0':   # multiple temperatures
                        self.temp0 = np.array(list(map(float,value.split(',')[:-1])), dtype='float64')
                    elif variable == 'tempi':   # multiple temperatures
                        self.tempi = np.array(list(map(float,value.split(',')[:-1])), dtype='float64')
                    elif variable == 'duration':    # multiple temperatures
                        self.duration = np.array(list(map(float,value.split(',')[:-1])), dtype='int32')
                    elif variable == 'tautp':
                        self.tautp = float(value[:-1])
                    elif variable == 'iCircular':
                        self.iCircular = int(value[:-1])
                    elif variable == 'verlet_cut':
                        self.verlet_cut = float(value[:-1])
                    elif variable == 'nVerlet':
                        self.nVerlet = int(value[:-1])
                    elif variable == 'dpr':
                        self.dpr = int(value[:-1])
                    elif variable == 'iCollide':
                        self.iCollide = int(value[:-1])
#====== ForceParms ===============================================================
class ForceParms:
    pass
    #--------------------------------------------------------------------------
    def __init__(self, inp=None):
        if isinstance(inp, str):
            if inp[-3:] == '.in':
                self.read_fparm(inp)
            else:
                print('Unsupported filename suffix for ForceParms:', inp)
                print('Standard filename suffix for ForceParms: xxxx.in')
                sys.exit(0)
        else:
            print('wrong file name format', inp)
            sys.exit(0)
    
    #---------- Force parameters --------------------------------------------------------
    def read_fparm(self, inp):
        #--- assign the default values ---------------
        self.sigma = 9.8
        self.epsilon = 0.9
        self.vcut = 13.5
        self.bphm_eq = 10.0
        self.bphm_k = 0.3
        self.bphm_p = -2.0
        self.angle_eq = 3.1415926535898
        self.angle_k = 4.0
        self.angle_p = -3.7
        self.bond_eq = 5.5
        self.bond_k = 20.0
        self.charge = -2.93
        self.ecut = 20.0
        self.len_Debye = 10.0
        self.min_loop = 4
        self.minbpFLAG = 2
        self.minStembp = 2
        self.bulge_p = 0.02
        self.bcut = 14.0
        self.stk_eq = 1.57079633
        self.stk_k = [0.9,  1.1,  2.1,  2.2,  0.6,  1.4, 1.3,  0.9,  2.1,  2.4,  1.0,  1.3, 2.4,  2.2,  3.3,  3.4,  1.5,  2.5, 2.1,  2.1,  2.4,  3.3,  1.4,  2.1, 1.3,  1.4,  2.1,  2.5,  0.5, -1.3, 1.0,  0.6,  1.4,  1.5, -0.3,  0.5]
        self.stk_ks = 1.5
        self.stk_p = -3.0
        self.bpFLAG = [1, 4, 1, 2, 4, 1, 1, 1, 6, 1]
        self.nHB = [0, 2, 0, 0, 2, 0 ,0 ,0, 3, 0]
        self.f_scale = 1.0
        #---- Read and undate the values ------------
        lines = open(inp).readlines()
        for i in range(len(lines)):
            if lines[i].strip() == '&parm':
                #print(lines[i])
                for j in range(i+1, len(lines)):
                    if lines[j].strip() == '/':
                        break
                    variable, value = lines[j].split('=')           
                    if value[-2] != ',':
                        print('No blank space is allow after the value:', lines[j])
                        sys.exit(0)
                    variable = variable.strip()
                    value = value.strip()     
                    if variable == 'sigma':
                        self.sigma = float(value[:-1])
                    elif variable == 'epsilon':
                        self.epsilon = float(value[:-1])
                    elif variable == 'vcut':
                        self.vcut = float(value[:-1])
                    elif variable == 'bphm_eq':
                        self.bphm_eq = float(value[:-1])
                    elif variable == 'bphm_k':
                        self.bphm_k = float(value[:-1])
                    elif variable == 'bphm_p':
                        self.bphm_p = float(value[:-1])
                    elif variable == 'angle_eq':
                        self.angle_eq = float(value[:-1])
                    elif variable == 'angle_k':
                        self.angle_k = float(value[:-1])
                    elif variable == 'angle_p':
                        self.angle_p = float(value[:-1])
                    elif variable == 'bond_eq':
                        self.bond_eq = float(value[:-1])
                    elif variable == 'bond_k':   # multiple temperatures
                        self.bond_k = float(value[:-1])
                    elif variable == 'charge':   # multiple temperatures
                        self.charge = float(value[:-1])
                    elif variable == 'ecut':    # multiple temperatures
                        self.ecut = float(value[:-1])
                    elif variable == 'len_Debye':
                        self.len_Debye = float(value[:-1])
                    elif variable == 'min_loop':
                        self.min_loop = int(value[:-1])
                    elif variable == 'minbpFLAG':
                        self.minbpFLAG = int(value[:-1])
                    elif variable == 'minStembp':
                        self.minStembp = int(value[:-1])
                    elif variable == 'bulge_p':
                        self.bulge_p = float(value[:-1])
                    elif variable == 'bcut':
                        self.bcut = float(value[:-1])
                    elif variable == 'stk_eq':
                        self.stk_eq = float(value[:-1])
                    elif variable == 'stk_k':     # Turner's parameters
                        self.stk_k = np.array(list(map(float,value.split(',')[:-1])), dtype='float64')
                    elif variable == 'stk_p':
                        self.stk_p = float(value[:-1])
                    elif variable == 'stk_ks':
                        self.stk_ks = float(value[:-1])
                    elif variable == 'bpFLAG':      # The number of HBs for each BP type
                        self.bpFLAG = np.array(list(map(int,value.split(',')[:-1])), dtype='int32')
                    elif variable == 'nHB':      # The number of HBs for each BP type
                        self.nHB = np.array(list(map(int,value.split(',')[:-1])), dtype='int32')
                    elif variable == 'f_scale':
                        self.f_scale = float(value[:-1])

#====== Restraint ===============================================================
class Restraint:
    rBplist = []
    nrBp = []

    #--------------------------------------------------------------------------
    def __init__(self, inp, minloop):
        if isinstance(inp, str):
            if inp[-3:] == '.in':
                self.read_rsn(inp, minloop)
            else:
                print('Unsupported filename suffix for Restraint:', inp)
                print('Standard filename suffix for Restraint: restraint.in')
                sys.exit(0)
        else:
            print('wrong file name format', inp)
            sys.exit(0)
    def read_rsn(self, inp, minloop):
        self.nrBp, self.rBplist = dot2bps(inp, minloop)

#====== Systen ================================================================
class System:

    #--------------------------------------------------------------------------
    def __init__(self, fntop, fncrd, fninp, fnprm, fnrsn_ntr, fnrsn_dpr, fnout, fntrj, fnrst):
        self.read_top(fntop)
        self.read_crd(fncrd)
        self.read_ctrl(fninp)
        self.read_parm(fnprm)
        self.read_ntr(fnrsn_ntr, self.fparm.min_loop)
        self.read_dpr(fnrsn_dpr, self.fparm.min_loop)
        self.fnout = fnout
        self.fntrj = fntrj
        self.fnrst = fnrst

    #----- Topology --------------------------------------------------------------
    def read_top(self, fn):
        self.top = Prmtop(fn)

    #----- Coordinates and velocities --------------------------------------------
    def read_crd(self, fn):
        self.crd = Restart(fn)
        
    #----- Control parameters --------------------------------------------------
    def read_ctrl(self, fn):
        self.ctrl = Controls(fn)
    
    #----- Force fiele parameters ----------------------------------------------
    def read_parm(self, fn):
        self.fparm = ForceParms(fn)

    #------- Restraints (ntr)----------------------------------------------------------
    def read_ntr(self, fn, minloop):
        self.ntr = Restraint(fn, minloop)   # use Prmtop class or include them directly?
    #------- Restraints (dpr)----------------------------------------------------------
    def read_dpr(self, fn, minloop):
        self.dpr = Restraint(fn, minloop)   # use Prmtop class or include them directly?

    def run(self):
        # Call the C/C++ func produced by boost library
        toplist = [np.array([self.top.pointers],dtype='int32'),
                   np.array([self.top.charge]),
                   np.array([self.top.amass]),
                   np.array([self.top.iac], dtype='int32'),
                   np.array([self.top.ico], dtype='int32'),
                   np.array([self.top.inb], dtype='int32'),
                   self.top.nNBres]
        rstlist = [self.crd.box[0:6], self.crd.crd, self.crd.vel]
        #print('coordinates:\n', self.crd.crd)
        #print('velocities:\n', self.crd.vel)
        ctrlist = [self.ctrl.ntpr, self.ctrl.ntwx,
                   self.ctrl.ntr, self.ctrl.ismd,
                   self.ctrl.dumpfreq, self.ctrl.stretch_force, self.ctrl.stretch_spring,
                   self.ctrl.dist_init, self.ctrl.dist_final,
                   self.ctrl.restraint_wt, self.ctrl.restraint_eq, self.ctrl.rcut,
                   self.ctrl.nscm, self.ctrl.nstlim, self.ctrl.dt,
                   self.ctrl.ntt, self.ctrl.temp0, self.ctrl.tempi,
                   self.ctrl.duration, self.ctrl.ntemp, self.ctrl.tautp,
                   self.ctrl.iCircular,
                   self.ctrl.verlet_cut, self.ctrl.nVerlet,
                   self.ctrl.dpr,
                   self.ctrl.iCollide,
                   ]
        if self.ctrl.ntemp != len(self.ctrl.duration):
            print("Error !!!")
            print("ntemp != len(duration)")
            print("Please change the temperature controlling parameters !!!")
            sys.exit(1)
        if self.ctrl.ntemp != len(self.ctrl.tempi):
            print("Error !!!")
            print("ntemp != len(tempi)")
            print("Please change the temperature controlling parameters !!!")
            sys.exit(1)
        if self.ctrl.ntemp != len(self.ctrl.temp0):
            print("Error !!!")
            print("ntemp != len(temp0)")
            print("Please change the temperature controlling parameters !!!")
            sys.exit(1)
        #print('ctrlist:', ctrlist)
        prmlist = [self.fparm.sigma, self.fparm.epsilon, self.fparm.vcut,
                   self.fparm.bphm_eq, self.fparm.bphm_k, self.fparm.bphm_p,
                   self.fparm.angle_eq, self.fparm.angle_k, self.fparm.angle_p,
                   self.fparm.bond_eq, self.fparm.bond_k,
                   self.fparm.charge, self.fparm.ecut, self.fparm.len_Debye,
                   self.fparm.min_loop, self.fparm.minbpFLAG, self.fparm.minStembp,
                   self.fparm.bulge_p, self.fparm.bcut,
                   self.fparm.stk_eq, self.fparm.stk_k, self.fparm.stk_ks,
                   self.fparm.stk_p, self.fparm.bpFLAG, self.fparm.nHB, self.fparm.f_scale]
        # check the cutoffs
        if self.fparm.ecut < self.fparm.vcut or self.fparm.ecut < self.fparm.bcut:
            print("Error !!!")
            print('ecut{0} must be larger than vcut{1} and bcut{2}'.format(self.fparm.ecut, self.fparm.vcut, self.fparm.bcut))
            print("Please change the cutoff parameters !!!")
            sys.exit(1)
        #print('nchain:\n', self.rsn.nchain)
        #print('maxHelix:\n', self.rsn.maxHelix)
        #print('ntr_nBp:\n', self.ntr.nrBp)
        #print('ntrBplist:\n', np.array(self.ntr.rBplist))
        #print('dpr_nBp:\n', self.dpr.nrBp)
        #print('dprBplist:\n', np.array(self.dpr.rBplist))
        ntrlist = [self.ntr.nrBp, np.array([self.ntr.rBplist], dtype='int32')]
        dprlist = [self.dpr.nrBp, np.array([self.dpr.rBplist], dtype='int32')]
        ext.run_md(toplist, rstlist, ctrlist, prmlist, ntrlist, dprlist, self.fnout, self.fntrj, self.fnrst)
