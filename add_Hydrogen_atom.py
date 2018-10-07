#wensong
#2018.10.7
#define a function to add hydrogen atoms for pdb file in that to find Hydrogen bond
import phenix
import os

def add_hydrogen_atoms():

   os.system("phenix.reduce %s > %s " %('1b20.pdb' ,'1b20'+'h'+'.pdb'))


#test
pdb_file = "1b20.pdb"

add_hydrogen_atoms()