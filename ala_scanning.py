#!/usr/bin/env python3
"""Script to perform in silico Ala scanning
   Usage: ala_scanning.py pdb_file > log
   Uses Amber based Residue Library and forcefields
"""
import sys
import math
from Bio.PDB.PDBParser import PDBParser
from forcefield import VdwParamset
from residue_library import ResiduesDataLib

backbone_atoms=['N','C','CA','O','HA','NH']

def interaction_energy(r1,r2, backbone=False, vdw=True, elec=True):
    energy=0.
    for at1 in r1.get_atoms():
        for at2 in r2.get_atoms():
            if not backbone and (at1.id in backbone_atoms or at2.id in backbone_atoms):
                continue
            energy += int_at_energy(at1,at2, vdw=vdw, elec=elec)
    return energy

def int_at_energy(at1,at2, vdw=True, elec=True):
    e = 0
    if vdw:
        e += vdw_at_energy(at1,at2)
    if elec:
        e += elec_at_energy(at1,at2)
    return e

def vdw_at_energy(at1,at2):
    d = at1-at2
    return math.sqrt(at1.vdw.eps * at1.vdw.eps) * ((at1.vdw.sig/d)**12 -(at1.vdw.sig/d)**6)

def elec_at_energy(at1,at2):
    return 332.6 * at1.charge * at2.charge / (at1-at2)



# MAIN

# Check that the number of paramters, and send warning in none
if not len(sys.argv) > 1:
    print ("Usage: ala_scanning.py pdb_file")
    sys.exit(1)
# set the pdb_path and load the structure
pdb_path = sys.argv[1]
# Setting the Bio.PDB.Parser object
parser = PDBParser(PERMISSIVE=1)
# Loading structure
try:
    st = parser.get_structure('st', pdb_path)
except OSError:
    print ("Error: file not found ")
    sys.exit(1)
# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib('data/aaLib.lib')
# loading VdW parameters
ff_params = VdwParamset('data/vdwprm')
# assign data types, and charges
total_charge=0
for at in st.get_atoms():
    resname = at.get_parent().get_resname()
    params = residue_library.get_params(resname,at.id)
    if not params:
        sys.exit(1)
    at.atom_type = params.at_type
    at.charge = params.charge
    at.vdw = ff_params.at_types[at.atom_type]
    total_charge += at.charge
print ("Total Charge: ",total_charge)
# Calculating res-res interaction energy
# First add an index, 
i=0
for r in st.get_residues():
    r.index = i
    i+=1

vdw_energies=i*[0.]
elec_energies=i*[0.]
total_vdw_energy = 0.
total_elec_energy = 0.
for r1 in st.get_residues():
    vdw_energies[r1.index]=i*[0.]
    elec_energies[r1.index]=i*[0.]
    for r2 in st.get_residues():
        if r2.index<=r1.index:
            continue
        vdw_energies[r1.index][r2.index] =  interaction_energy(r1, r2, elec=False)
        elec_energies[r1.index][r2.index] =  interaction_energy(r1, r2, vdw=False)
        total_vdw_energy += vdw_energies[r1.index][r2.index]
        total_elec_energy += elec_energies[r1.index][r2.index]
# Reference values for WT
print (total_vdw_energy, total_elec_energy)
# Ala scanning DDG is the negative of interaction energies of side chain atoms except CB
ddg_energies=i*[0.]
for r1 in st.get_residues():
    for r2 in st.get_residues():
        if r1.get_resname() in ['GLY','ALA'] or r2.index <= r1.index:
            continue
        for at1 in r1.get_atoms():
            if at1.id in backbone_atoms or at1.id in ['CB', 'HB1', 'HB2', 'HB3']:
                continue
            for at2 in r2.get_atoms():
                ddg_energies[r1.index] += vdw_at_energy(at1,at2) + elec_at_energy(at1,at2)
    print (r1,ddg_energies[r1.index])
