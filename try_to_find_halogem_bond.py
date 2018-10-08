#wensong
#2018.10.7

from bond_distance import distance
import bond_angle
from change_pdb_into_dataframe import change_the_pdb_file_into_dataframe


def try_to_find_halogen_bond(pdb_file):
    i = 0
    j = 0
    Halogen_bond_pdb_file_atom_label =change_the_pdb_file_into_dataframe(pdb_file)
    atom1 = []
    halogen_atom = [[],[]]
    for Atomic_symbol in Halogen_bond_pdb_file_atom_label.loc['Atomic_symbol',:]:
        if Atomic_symbol == 'BR':
             atom1 = Halogen_bond_pdb_file_atom_label.loc[['X','Y','Z'],i].astype('float16')
             atom1 = list(atom1.values)
             halogen_atom[j] = atom1
             j = j + 1
        i = i+1

    k = 0
    h = 0
    for t in range(len(halogen_atom)):
        for k in Halogen_bond_pdb_file_atom_label:
            atom2 = Halogen_bond_pdb_file_atom_label.loc[['X','Y','Z'],k].astype('float16')
            k = k+1
            atom2 = list(atom2.values)
            bond_lengh = distance(halogen_atom[t],atom2)
            if 1<=bond_lengh<=5:
                for h in Halogen_bond_pdb_file_atom_label:
                    atom3 = Halogen_bond_pdb_file_atom_label.loc[['X','Y','Z'],k].astype('float16')
                    atom3 = list(atom3.values)
                    bond_lengh = distance(atom2,atom3)
                    if 1<bond_lengh<5:
                         from bond_angle import angle
                         bond_angle = angle(halogen_atom[t],atom2,atom3)
                         if 90<bond_angle<180:
                             print 'find halogen bond !'



        
















#test

try_to_find_halogen_bond('5v7d.pdb')











