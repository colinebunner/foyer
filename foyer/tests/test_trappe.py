import os
import parmed as pmd
from foyer import Forcefield
from utils import atomtype

TRAPPE_UA = Forcefield(forcefield_files='/Users/cbunner/MoSDeF-Soft/foyer/foyer/forcefields/trappe-ua-tmp.xml')
TRAPPE_TESTFILES_DIR = "/Users/cbunner/MoSDeF-Soft/foyer/foyer/trappe_validation"

test_mols = ["1-octanol","1-pentanol","2-methyl-2-butanol","propane","cyclopentane"]

mol2_files = ["{}/{}/{}.mol2".format(TRAPPE_TESTFILES_DIR,mol2,mol2) for mol2 in test_mols]

def get_known_atom_types(molfile):
   
    # In
    # -----------------------------------------------
    # molfile: filename for TraPPE-UA .mol2 test file
    #
    # Out
    # -----------------------------------------------
    # known_atom_types: A dictionary with keys from 0 to nbeads-1,
    #  where nbeads is the number of beads in the test compound, and
    #  values of strings corresponding to the correct atom type as
    #  specified in the TraPPE-UA .xml force field.

    known_atom_types = {}
   
    with open(molfile,"r") as f:
        f = [line.strip() for line in f.readlines()]
        i = 0
        # Controls reading of mol2 file
        atom = False

        # read .mol2 file
        while (i < len(f)):
            # Find atoms section
            if f[i].startswith("@<TRIPOS>ATOM"):
                atom = True
                i += 1
                while atom == True:
                    if f[i].startswith("@<TRIPOS>"):
                        i += 1
                        atom = False
                    else:
                        # 1st column is bead number, 6th column is explicitly defined atom type
                        known_atom_types[int(f[i].split()[0])-1] = f[i].split()[5]
                        i += 1
            else:
                i += 1

    return known_atom_types

if __name__ == '__main__':
    for molfile in mol2_files:
        try:
           # .mol2 to ParmEd Structure instance
           structure = pmd.load_file(molfile, structure=True)
	   # apply forcefield to structure
           structure = TRAPPE_UA.apply(structure,assert_angle_params=False,assert_dihedral_params=False,
			                assert_improper_params=False) 

           # Get explicitly coded atom types from .mol2 files
           known_atom_types = get_known_atom_types(molfile)

	   # Check foyer atom types from .xml file and bead assignment algorithm
           # against the known atom type
           print("Atoms:")
           for i,atom in enumerate(structure.atoms):
               print('Atom {} is typed as {} (true={})'.format(atom, atom.type,known_atom_types[i]))
        # Check atom types against correct atom types as listed in .mol2
        except Exception as e:
            print(e)

