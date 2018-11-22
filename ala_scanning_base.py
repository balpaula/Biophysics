#!/usr/bin/env python3
"""Script to perform in silico Ala scanning
   Usage: ala_scanning.py pdb_file > log
   Uses Amber based Residue Library and forcefields
"""
import sys
import math
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset
from residue_library import ResiduesDataLib

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
for at in st.get_atoms():
    resname = at.get_parent().get_resname()
    params = residue_library.get_params(resname,at.id)
    if not params:
        sys.exit(1)
    at.xtra['atom_type'] = params.at_type
    at.xtra['charge'] = params.charge
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

# Calculating surfaces
srf = NACCESS_atomic(st[0],naccess_binary ='PATH_TO_NACCESS' )

