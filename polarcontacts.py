#! /usr/bin/python
#
# Simple parser to extact topology of NA from simulation snapshot
# Chain ids or TER not required 
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

import os
import sys
import argparse


COVLNK = 2.0
HBLNK  = 3.5

all_polars = [
    'N','ND1','ND2','NE','NE1','NE2','NH1','NH2', 'NZ', 
    'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH',
    'S', 'SD', 'SG']
backbone_polars =  ['N','O']
waternames = ['WAT','HOH']

#Appropriate wrapper over BioPython classes

class Residue():
    oneLetterResidueCode = {
        'DA' :'A', 'DC' :'C', 'DG' :'G', 'DT' :'T',
        'A'  :'A', 'C'  :'C', 'G'  :'G', 'U'  :'U',
        'DA3':'A', 'DC3':'C', 'DG3':'G', 'DT3':'T',
        'A3' :'A', 'C3' :'C', 'G3' :'G', 'U3' :'U',
        'DA5':'A', 'DC5':'C', 'DG5':'G', 'DT5':'T',
        'A5' :'A', 'C5' :'C', 'G5' :'G', 'U5' :'U',
        'MRA':'A',
        'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 
        'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 
        'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 
        'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
    }  

    def __init__(self, r, useChains=False):
        self.residue = r
        self.useChains = useChains
    
    def resid(self, compact=False):
        if self.useChains:
            ch = ":"+self.residue.get_parent().id
        else:
            ch = ''
        if compact:
            return self._getOneLetterResidueCode() + ch + str(self.residue.id[1])
        else:
            return self.residue.get_resname() + ch + ':'+ str(self.residue.id[1])       

    def resNum(self):
        return self.residue.id[1]

    def _getOneLetterResidueCode(self):
        id = self.residue.get_resname().rstrip().lstrip()
        if not id in Residue.oneLetterResidueCode:
            return 'X'
        else:
            return Residue.oneLetterResidueCode[id]
    
    def __hash__(self):
        return hash(self.resid())
    
    def __eq__(self, other):        
        if self.useChains:
            id = self.residue.get_parent().id + str(self.residue.id[1])
            otherid = other.residue.get_parent().id + str(other.residue.id[1])
        else:
            id = self.residue.id[1]
            otherid = other.residue.id[1]
        return id == otherid

    def __lt__(self, other):        
        if self.useChains:
            id = self.residue.get_parent().id + str(self.residue.id[1])
            otherid = other.residue.get_parent().id + str(other.residue.id[1])
        else:
            id = self.residue.id[1]
            otherid = other.residue.id[1]
        return id < otherid
    
    def __str__(self):
        return self.resid()
            

class Atom():
    def __init__ (self,at,useChains=False):
         self.at=at
         self.useChains=useChains
    
    def atid(self, compact=False):
        return self.resid(compact)+"."+self.at.id
    
    def resid(self, compact=False):
        return Residue(self.at.get_parent(),self.useChains).resid(compact)
    
    def resNum(self):
        return Residue(self.at.get_parent(),self.useChains).resNum()
  
    def attype(self):
        return Residue(self.at.get_parent(),self.useChains)._getOneLetterResidueCode()+'.'+self.at.id
    
    def __lt__(self,other):
        return self.at.get_serial_number()
    
    def __str__(self):
        return self.atid()

    
def main():
  
    parser = argparse.ArgumentParser(
                prog='polarContacts', 
                description='Polar contacts detector'
            )

    parser.add_argument(
        '--debug', '-d', 
        action='store_true', 
        dest='debug',
        help='Produce DEBUG output'
    )
    
    parser.add_argument(
        '--backonly', 
        action='store_true', 
        dest='backonly',
        help='Restrict to backbone'
    )
    parser.add_argument(
        '--nowats', 
        action='store_true', 
        dest='nowats',
        help='Exclude water molecules'
    )

    parser.add_argument('pdb_path')
    
    args = parser.parse_args()
    
    debug = args.debug
    backonly = args.backonly
    nowats =args.nowats
    pdb_path = args.pdb_path
    
    if not pdb_path:
        parser.print_help()
        sys.exit(2)        

    parser = PDBParser(PERMISSIVE=1)
    
    try:
        st = parser.get_structure('st', pdb_path)
    except OSError:
        print ("#ERROR: loading PDB")
        sys.exit(2)

# Checking for models
    if len(st) > 1:
        print ("#WARNING: Several Models found, using only first")

# Using Model 0 any way
    st = st[0]   

# Making a list of polar atoms
    polats = []
    if backonly:
        selected_atoms = backbone_polars
    else:
        selected_atoms = all_polars
        
    for at in st.get_atoms():
        if at.id in selected_atoms:
            polats.append(at)
#Searching for contacts under HNLNK on diferent residues            
    nbsearch = NeighborSearch(polats)  
    hblist = []
    for at1, at2 in nbsearch.search_all(HBLNK):
        if at1.get_parent() == at2.get_parent():
            continue
        if at1.get_serial_number() > at2.get_serial_number():
            continue
#Discard covalents and neighbours
        if (at1-at2) < COVLNK:
            continue
        if at2.get_parent().id[1] - at1.get_parent().id[1] == 1:
            continue
# remove waters
        if nowats:
            if at1.get_parent().get_resname() in waternames \
                or at2.get_parent().get_resname() in waternames:
                continue
        hblist.append([Atom(at1,1),Atom(at2,1)])
    for hb in sorted (hblist,key=lambda i: i[0].at.get_serial_number()):
        print ('{:14} {:14} {:6.2f}'.format(hb[0].atid(),hb[1].atid(),(hb[0].at-hb[1].at)))
        
    
    
    
if __name__ == "__main__":
    main()
