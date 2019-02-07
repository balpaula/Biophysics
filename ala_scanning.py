#!/usr/bin/env python3
''' Paula Balcells Falgueras
    Biophysics, Bachelor's Degree in Bioinformatics
    December 2018
'''

''' Original script: https://github.com/jlgelpi/Biophysics
    Original author: @jgelpi
'''

'''Script to perform in silico Ala scanning
   Usage: ala_scanning.py pdb_file > log
   Uses Amber based Residue Library and forcefields
'''

import sys
import math
from split_residues import *
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset
from residue_library import ResiduesDataLib

### SET UP

# Check that the number of paramters, and send warning in none
if not len(sys.argv) > 1:
    print ("Usage: ala_scanning.py pdb_file")
    sys.exit(1)

def set_structure(pdb_path):
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
        #print(at.id)
        if not params:
            sys.exit(1)
        at.xtra['atom_type'] = params.at_type
        at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

    # Calculating surfaces
    srf = NACCESS_atomic(st[0],naccess_binary ='PATH TO NACCESS' )

    return st

# set the pdb_path and load the structure
pdb_path = sys.argv[1]
st = set_structure(pdb_path)

##########################################

### FUNCTIONS DECLARATION

# Mehler-Solmajer dielectric
def ms_dielectric(r):
    return (86.9525/(1-7.7839*math.exp(-0.3153*r)))-8.5525

# Solvation energy
def calc_solvation(at):
    if 'EXP_NACCESS' not in at.xtra:
        return 0.
    else:
        return at.xtra['vdw'].fsrf * float(at.xtra['EXP_NACCESS'])

# Electrostatic interaction energy
def calc_elec(at, atoms, dielectric = None):
    total = 0.
    for at2 in atoms:
        # not considering intra residue interactions
        if at.get_parent() == at2.get_parent():
            continue
        else:
            r = at - at2
            if dielectric:
                d = dielectric
            else:
                d = ms_dielectric(r)
            total += 332.16*at.xtra['charge']*at2.xtra['charge']/(d*r)
    return total

# Van der Waals interaction energy
def calc_vdw(at1, atoms):
    total = 0.
    for at2 in atoms:
        # not considering intra residue interactions
        if at1.get_parent() == at2.get_parent():
            continue
        
        # Distance between atoms
        distance = at1 - at2

        # Avoiding covalent bonds
        if distance > 2.0:
            SIG = math.sqrt(at1.xtra['vdw'].sig * at2.xtra['vdw'].sig)
            EPS = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)

            f = SIG / distance
            total += 4. * EPS * (pow(f, 12)-pow(f, 6))
    return total

# Differences in total protein energy after Ala mutations


#################################

### WORKFLOW

residues = list(st.get_residues())

# Atoms composing the Ala residue
ala_atoms = ["N","H","CA","HA","CB","HB1","HB2","HB3","C","O"]

# Total energies
total_solv = 0.
total_elec = 0.
total_vdw = 0.

# Storing differences in total energy after mutation for the folded state(ordered by mutation position)
diff = []

# In case we wanted to use a fixed dielectric
diel = 80.

# Printing on the terminal
print("MUTATIONS IN FOLDED")
print()

# Iterating through the residues of the structure
for r in residues:

    # Residue energies
    res_solv_energy = 0.
    res_elec_energy = 0.
    res_vdw_energy = 0.

    # Iterating through the atoms of the residue
    for at in r.get_atoms():

        # Solvation energy for atom
        at_solv_energy = calc_solvation(at)
        res_solv_energy += at_solv_energy

        # Electrostatic interaction energy for atom
        at_elec_energy = calc_elec(at, st.get_atoms())
        # To use a fixed dielectric: calc_elec(at, st, dielectric = diel)
        res_elec_energy += at_elec_energy

        # Vdw interaction energy for atom
        at_vdw_energy = calc_vdw(at, st.get_atoms())
        res_vdw_energy += at_vdw_energy

    print(r)
    print("WT")
    print("Elec: "+str(res_elec_energy))
    print("Solv: "+str(res_solv_energy))
    print("VdW: "+str(res_vdw_energy))
    print()


    # Update total energies
    total_solv += res_solv_energy
    total_elec += res_elec_energy
    total_vdw += res_vdw_energy


    # Select those atoms from the residue that are different from Alanine
    different_atoms = []
    for at in r.get_atoms():
        if at.name not in ala_atoms:
            different_atoms.append(at)

    
    # To get the energies of the mutated residue, substract those interactions implying the different atoms
    # Initialized to the wild type residue
    res_solv_energy_mut = res_solv_energy
    res_elec_energy_mut = res_elec_energy
    res_vdw_energy_mut = res_vdw_energy

    # Interactions from the atoms of the mutated residue
    for at in different_atoms:
        at_solv_energy = calc_solvation(at)
        at_elec_energy = calc_elec(at, st.get_atoms())
        at_vdw_energy = calc_vdw(at, st.get_atoms())
        
        res_solv_energy_mut -= at_solv_energy
        res_elec_energy_mut -= at_elec_energy
        res_vdw_energy_mut -= at_vdw_energy

    # Interactions from the atoms of the other residues to the different atoms of the mutated residue
    for r2 in residues:
        for at in r2.get_atoms():
            at_elec_energy = calc_elec(at, different_atoms)
            res_elec_energy_mut -= at_elec_energy

            at_vdw_energy = calc_vdw(at, different_atoms)
            res_vdw_energy_mut -= at_vdw_energy

    print("MUT")
    print("Elec: "+str(res_elec_energy_mut))
    print("Solv: "+str(res_solv_energy_mut))
    print("VdW: "+str(res_vdw_energy_mut))
    print()    

    # Total energy of the wild type residue
    res_total = res_solv_energy + res_elec_energy + res_vdw_energy
    # Total energy of the mutated residue
    res_total_mut = res_solv_energy_mut + res_elec_energy_mut + res_vdw_energy_mut
    # Difference in energies
    res_diff = res_total_mut - res_total

    #print(r)
    print("Difference: "+str(res_diff))
    print()
    diff.append(res_diff)


# Total energy of the wild type protein
total_energy_protein = total_solv + total_elec + total_vdw

# Delta delta G for X -> Ala in the folded state
DDGFala = diff

# Split the complete protein PDB file into separated files according to residues
# Using the function declared in split_residues.py file
split_residues(pdb_path,st)

# Storing differences in total energy after mutation for the unfolded state(ordered by mutation position)
diff2 = []

print()
print("MUTATIONS IN FOLDED")
print()

for i in range(1,len(residues)+1):

        # Get the PDB file corresponding to the residue position
        aa_path = pdb_path[:-4]+"_unfolded/"+str(i)+".pdb"

        # Set as the structure the individual residue
        st_aa = set_structure(aa_path)

        res_solv_energy_unfolded = 0.
        res_solv_energy_unfolded_mut = 0.

        # Calculate the unfolded energy for the residue
        for at in st_aa.get_atoms():

            # Only solvation energy is considered
            solv_at_unfolded = calc_solvation(at)
            res_solv_energy_unfolded += solv_at_unfolded

            # For the mutated residue, only the solvation energy for common atoms between the residue and alanine are considered
            if at.name in ala_atoms:
                res_solv_energy_unfolded_mut += solv_at_unfolded

        # Difference in energies
        res_diff_unfolded = res_solv_energy_unfolded_mut - res_solv_energy_unfolded

        r = list(st_aa.get_residues())

        print(r[0])
        print("Difference: "+str(res_diff_unfolded))
        print()

        diff2.append(res_diff_unfolded)

# Delta delta G for X -> Ala in the unfolded state
DDGUala = diff2

# Calculate differences on proteing folding energy
DDGf = []

for i in range(len(DDGFala)):
    DDGf.append(DDGFala[i]-DDGUala[i])

# Print a summary of the total energies for the wild type protein
print("WILD TYPE PROTEIN ENERGIES")
print()
print("Solvation energy: "+str(total_solv))
print("Electrostatic energy: "+str(total_elec))
print("Vdw energy: "+str(total_vdw))
print("----------")
print("Total energy: "+str(total_energy_protein))
print()

# Plot changes in energy of folded state after Ala mutation
plt.plot(DDGFala)
plt.title("Change in energy of the folded state after Ala mutation")
plt.xlabel("Position of the residue")
plt.ylabel("DDG")
plt.show()

# Plot changes in energy of unfolded state after Ala mutation
plt.plot(DDGUala)
plt.title("Change in energy of the unfolded state after Ala mutation")
plt.xlabel("Position of the residue")
plt.ylabel("DDG")
plt.show()

# Calculate differences on proteing folding energy
DDGf = []

for i in range(len(DDGFala)):
    DDGf.append(DDGFala[i]-DDGUala[i])

# Plot changes in protein folding energy after Ala mutation
plt.plot(DDGf)
plt.title("Change in protein folding energy after Ala mutation")
plt.xlabel("Position of the residue")
plt.ylabel("DDG folding")
plt.show()

