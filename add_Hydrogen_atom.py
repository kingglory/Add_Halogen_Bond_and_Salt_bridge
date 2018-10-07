#wensong
#2018.10.7
#define a function to add hydrogen atoms for pdb file in that to find Hydrogen bond
import phenix


def add_hydrogen_atoms(pdb_file, pdb_add_H_file=None):
    phenix.reduce (pdb_file, pdb_add_H_file)
    return pdb_add_H_file


#test
pdb_file = "1b20.pdb"
add_hydrogen_atoms(pdb_file)