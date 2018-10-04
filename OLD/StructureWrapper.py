#
# Convenient wrappers over BioPython classes
# Hint: Include here interaction energies
#

__author__ = "gelpi"
__date__ = "$01-nov-2017 7:17:33$"

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
        if self.useChains:
            resnum = self.residue.get_parent().id + str(self.residue.id[1])
        else:
            resnum = self.residue.id[1]
        return resnum

    def _getOneLetterResidueCode(self):
        resid = self.residue.get_resname().rstrip().lstrip()
        if not resid in Residue.oneLetterResidueCode:
            return 'X'
        else:
            return Residue.oneLetterResidueCode[resid]

    def __hash__(self):
        return hash(self.resid())

    def __eq__(self, other):
        return self.resid() == other.resid()

    def __lt__(self, other):
        return self.resnum() < self.resnum()

    def __str__(self):
        return self.resid()

class Atom():
    def __init__ (self, at, useChains=False):
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
        return self.at.get_serial_number()<other.at.get_serial_number()

    def __str__(self):
        return self.atid()
