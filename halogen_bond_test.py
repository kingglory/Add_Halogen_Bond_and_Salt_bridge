from __future__ import division
from libtbx import easy_run
import os
from run import hierarchy_No_cif_model
from run import hierarchy_cif_model
from run import find_the_atoms_makeing_up_halogen_bond
# prepare the cif file if the pdb file needs

def add_H_atoms_into_pad_files():
  X_bonds_file = ["5v7d.pdb","2h79.pdb", "2ito.pdb", "2oxy.pdb","2vag.pdb",
                  "2yj8.pdb", "3v04.pdb", "4e7r.pdb"]
  for pdb_file in X_bonds_file:
      add_h_pdb_file = pdb_file[0:4] + 'h.pdb'
      easy_run.call("phenix.reduce %s > %s " % (pdb_file,  add_h_pdb_file))
      
def prepare_cif_for_pdb_file():
  X_bonds_file = ["5v7dh.pdb","2h79h.pdb", "2itoh.pdb", "2oxyh.pdb","2vagh.pdb",
                  "2yj8h.pdb", "3v04h.pdb", "4e7rh.pdb"]
  for pdb_file in X_bonds_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 X_bonds_file = ["5v7dh.pdb","2h79h.pdb", "2itoh.pdb", "2oxyh.pdb","2vagh.pdb",
                  "2yj8h.pdb", "3v04h.pdb", "4e7rh.pdb"]

 i = 0
 for pdb_file in X_bonds_file:
  pdb_file = pdb_file[0:5] + ".updated.pdb"
  pdb_cif = pdb_file[0:5] + ".ligands.cif"
  if os.path.exists(pdb_cif):
    X_bonds_file[i] = [pdb_file, pdb_cif]
  else:
    X_bonds_file[i] = [pdb_file,None]
  i = i + 1
 return  X_bonds_file


def find_the_atoms_makeing_up_halogen_bond_test():
    X_bonds_file = list_cif_and_pdb_file()
    for i in range(len(X_bonds_file)) :
      pdb_file = X_bonds_file[i][0]
      if X_bonds_file[i][1] == None:
        hierarchy, vdwr,model=hierarchy_No_cif_model(pdb_file)
        find_the_atoms_makeing_up_halogen_bond(hierarchy, vdwr,model)
        print (pdb_file[0:4] + " this pdb_file is ok")
      else:
        hierarchy, vdwr,model = hierarchy_cif_model(pdb_file)
        find_the_atoms_makeing_up_halogen_bond(hierarchy, vdwr,model)
        print (pdb_file[0:4] + " this pdb_file_cif is ok")




if __name__ == '__main__':
    #add_H_atoms_into_pad_files()
    #prepare_cif_for_pdb_file()
    find_the_atoms_makeing_up_halogen_bond_test()