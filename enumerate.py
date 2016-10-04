import itertools as it
import os

class Structure:
    def __init__(self):
        self.system = None
        self.lat = None
        self.elems = None
        self.natoms = None
        self.totatoms = None
        self.coord_system = None
        self.coords = None

def read_poscar(fposcar):
    fr = open(fposcar, 'r')
    poscar = [line.strip() for line in fr.readlines()]
    s = Structure()
    s.system = poscar[0]
    s.lat = [[float(i) for i in poscar[j].split()] for j in range(2,5)]
    s.elems = poscar[5].split()
    s.elems.pop(-1)
    s.natoms = [int(n) for n in poscar[6].split()]
    s.totatoms = sum(s.natoms)
    s.coord_system = poscar[7]
    s.coords = [[float(i) for i in poscar[j].split()] for j in range(8,8+s.totatoms)]
    seqs = generate_combinations(s.natoms[-1])
    print seqs
    write_poscar('pos', s)
    return s

def generate_combinations(n):
    return ["".join(seq) for seq in it.product("01", repeat=n)]

def write_poscar(fposcar, struct):
    ##fw = open(fposcar, 'w')
    ##fw.write(struct.system + '\n')
    print "writing:", struct.system

def enum(fposcar):
    struct = read_poscar(fposcar)
    print struct

if __name__=='__main__':
    enum('template.vasp')

