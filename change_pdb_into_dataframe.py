#wensong
#2018.10.6
#this is a function that can take the atom informention out from the pdb file into a pandas dataframe
from Find_Halogen_atom import find_halogen_atoms
from is_number import is_number
import pandas as pd
def change_the_pdb_file_into_dataframe(pdb_file):
    i = 0
    pdb_inp = open(pdb_file)
    Halogen_bond_pdb_file_atom_label = pd.DataFrame(index=['Atom_ID','Chain_label', 'Nucleic_acid_ID', 'X', 'Y', 'Z', 'Atomic_symbol'])
    if (find_halogen_atoms(pdb_file)==True):
        for line in pdb_inp.readlines():
            L = line.split(' ')
            L = [item for item in filter(lambda x: x != '', L)]
            if (is_number(L[1])) == True:
                if len(L)>=13:
                    if (is_number(L[5])) == True:
                        if (is_number(L[6])) == True:
                            if (is_number(L[7])) == True:
                                Halogen_bond_pdb_file_atom_label[i] = [L[1], L[4], L[5], L[6], L[7], L[8], L[11]]
                                i = i + 1
        print Halogen_bond_pdb_file_atom_label

    return Halogen_bond_pdb_file_atom_label





# test the change_the_pdb_file_into_dadaframe(pdb_file)
pdb_file='5v7d.pdb'
change_the_pdb_file_into_dataframe(pdb_file)


