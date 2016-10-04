import os
import shutil
import sys
import subprocess
import glob
import copy
import collections
from collections import OrderedDict
import itertools


class Structure:
    def __init__(self):
        self.system = None
        self.scale = None
        self.lat = []
        self.elems = []
        self.natoms = []
        self.totatoms = None
        self.comp_dict = OrderedDict()
        self.atoms = []
        self.coord_system = None
        self.coords = []


def read_poscar(fposcar):
    fr = open(fposcar, 'r')
    poscar = [line.strip() for line in fr.readlines()]
    s = Structure()
    s.system = poscar[0]
    s.scale = float(poscar[1])
    s.lat = [[float(i) for i in poscar[j].split()] for j in range(2,5)]
    s.elems = poscar[5].split()
    s.natoms = [int(n) for n in poscar[6].split()]
    s.totatoms = sum(s.natoms)
    for k, v in zip(s.elems, s.natoms):
        s.comp_dict[k] = v
        for i in range(v):
            s.atoms.append(k)
    s.coord_system = poscar[7]
    s.coords = [[float(i) for i in poscar[j].split()] for j in range(8,8+s.totatoms)]
    return s


def generate_binary_strings(N):
    return ["".join(seq) for seq in itertools.product("01", repeat=N)]


def generate_mixing_configurations(struct, mixing_atom, seqs):
    mixing_sites = []
    for atom, coord in zip(struct.atoms, struct.coords):
        if atom != mixing_atom:
            continue
        mixing_sites.append(coord)
    configs = []
    for seq in seqs:
        config = collections.deque()
        for i, char in enumerate(seq):
            if char == '0':
                config.appendleft(mixing_sites[i])
            elif char == '1':
                config.append(mixing_sites[i])
        natoms = [seq.count('0'), seq.count('1')]
        configs.append((seq, natoms, config))
    return configs


def write_poscar(fposcar, s):
    fw = open(fposcar, 'w')
    fw.write("%s\n" %s.system)
    fw.write("%s\n" %s.scale)
    for a in s.lat:
        fw.write("%.15f  %.15f  %.15f\n" %(a[0], a[1], a[2]))
    elem_string = " ".join(s.elems)
    fw.write("%s\n" %elem_string)
    natoms_string = " ".join(map(str, s.natoms))
    fw.write("%s\n" %natoms_string)
    fw.write("%s\n" %s.coord_system)
    for c in s.coords:
        fw.write("%.15f  %.15f  %.15f\n" %(c[0], c[1], c[2]))


def remove_duplicate_structures():
    poscars = glob.glob('*/POSCAR')
    unique_structs = [poscars[0]]
    for p in poscars:
        unique = True
        for u in unique_structs:
            mint_call = ['mint', '%s'%p, '%s'%u, '-compare']
            mint_stdout = subprocess.check_output(mint_call)
            if "are the same" in mint_stdout:
                unique = False
                break
        if unique:
            unique_structs.append(p)

    for p in poscars:
        if p in unique_structs:
            continue
        folder = p.split('/')[0]
        shutil.rmtree(folder)


def enum(fposcar, mixing_atom='X', species=['E1', 'E2'], comp_e2=0.5, thresh=0.05):
    sys.stdout.write('reading structure from %s ... ' %fposcar)
    struct = read_poscar(fposcar)
    sys.stdout.write('done.\n')

    N = struct.comp_dict[mixing_atom]
    sys.stdout.write('number of mixing sites = %s \n' %N)
    sys.stdout.write('generating binary strings of length %s ... ' %N)
    seqs = generate_binary_strings(N)
    sys.stdout.write('done.\n')

    comps = [s.count('1')/float(N) for s in seqs]
    comp_seqs = []
    for seq, comp in zip(seqs, comps):
        if abs(comp - comp_e2) > thresh:
            continue
        comp_seqs.append(seq)
    sys.stdout.write("arrangements with given composition (%s) = %s\n" %(comp_e2, comp_seqs))

    sys.stdout.write('generating mixing configurations ... ')
    configs = generate_mixing_configurations(struct, mixing_atom, comp_seqs)
    sys.stdout.write('done.\n')

    sys.stdout.write('writing structures into POSCAR files ... ')
    for i, (seq, natoms, config) in enumerate(configs):
        struct2 = Structure()
        struct2.system = struct.system
        struct2.scale = struct.scale
        struct2.lat = struct.lat
        for e, n in zip(struct.elems, struct.natoms):
            if e == mixing_atom:
                struct2.elems.extend(species)
                struct2.natoms.extend(natoms)
            else:
                struct2.elems.append(e)
                struct2.natoms.append(n)
        struct2.coord_system = struct.coord_system

        mixing_site_count = 0
        for atom, coord in zip(struct.atoms, struct.coords):
            if atom == mixing_atom:
                struct2.coords.append(config[mixing_site_count])
                mixing_site_count += 1
            else:
                struct2.coords.append(coord)

        folder = "_".join(map(str, [i, seq]))
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)
        config_poscar = os.path.join(folder, 'POSCAR')
        write_poscar(config_poscar, struct2)
    sys.stdout.write('done.\n')
    
    sys.stdout.write('removing duplicates ... ')
    remove_duplicate_structures()
    sys.stdout.write('done.\n')
    

