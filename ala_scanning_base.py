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
print ("Total Charge: ",total_charge)

# Calculating res-res interaction energy
# Reference values for WT
# Ala scanning DDG is the negative of interaction energies of side chain atoms except CB
