#wensong
#2018.10.4
# write a function to find the Halogen atoms in one file
import iotbx.pdb
from iotbx.pdb import hierarchy
print "hello phenix"




# write a function to find the halogen atoms in pdb file
def find_halogen_atoms(pdb_file):
    pdb_in = hierarchy.input(file_name=pdb_file)
    for chain in pdb_in.hierarchy.only_model().chains() :
        for residue_group in chain.residue_groups() :
          for atom_group in residue_group.atom_groups() :
            for atom in atom_group.atoms() :
              if (atom.element.strip().upper() == "BR") :
                print('Here is a halogen atom !')
                return True
            if (atom_group.atoms_size() == 0) :
              residue_group.remove_atom_group(atom_group)
          if (residue_group.atom_groups_size() == 0) :
            chain.remove_residue_group(residue_group)

#test the find_halogen_atoms()

pdb_file = "5v7d.pdb"
find_halogen_atoms(pdb_file)