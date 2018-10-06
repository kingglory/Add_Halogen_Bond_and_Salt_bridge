#wensong
#unfinish!!!!!!
#2018.10.6
#this is a function that can take the atom informention out from the pdb file into a pandas dataframe
from docutils.parsers import null
from Find_Halogen_atom import find_halogen_atoms
from is_number import is_number
import iotbx.pdb
import pandas as pd
def change_the_pdb_file_into_dataframe(pdb_file):
    pdb_inp = open(pdb_file)
    print pdb_inp
    Halogen_bond_atom_label = pd.DataFrame(index=['Atom_ID', 'Chain_label', 'Nucleic_acid_ID', 'X', 'Y', 'Z', 'Atomic_symbol'])
    if (find_halogen_atoms(pdb_file)==True):
        for line in pdb_inp.readlines():
            print line
            L = line.split(' ')
            L = [item for item in filter(lambda x: x != '', L)]
            if (is_number(L[1])) == True:
                if len(L)>=6:
                    if (is_number(L[5])) == True:
                        print (L[1])
                        Halogen_bond_atom_label[0] = [L[1], L[4], L[5], L[6], L[7], L[8], L[11]]



#    awl = pdb_inp.atoms_with_labels()
#    for atom in awl:
#        print atom.xyz,atom.element



# test the change_the_pdb_file_into_dadaframe(pdb_file)

change_the_pdb_file_into_dataframe("5v7d.pdb")
