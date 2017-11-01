# BioPhysics
Source codes and data for BioPhysics course on Bioinformatics

## polarcontacts.py
Code for listing polar contacts ans interaction energies from a PDB file

## StructureWrapper.py
Convenient wrapper classes for Bio.PDB classes

## ForceField.py
Module for vdw parameters management

## ResLib.py
Module for Residue Library management

## data/vdwprm
Simple vdw parameters (based on AMBER ff)

## data/aaLib.lib
Library for obtaining amino acid atom types and partial charges (based on AMBER ff)


Usage:

    python3 polarContacts.py [-h] [--backonly] [--nowats] [--diel E] --vdw vdwprm_file --rlib resLib_file pdb_path


External dependencies

    Bio.PDB.NeighborSearch (BioPython)
    Bio.PDB.PDBParser (Biopython)

