from __future__ import division
import iotbx.pdb
import iotbx.cif
from iotbx.pdb import hierarchy
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
from libtbx import easy_run
import os
from run import halogen_find_test_No_cif_model
from run import halogen_find_test_cif_model
# prepare the cif file if the pdb file needs
def prepare_cif_for_pdb_file():
  X_bonds_file = ["5v7d.pdb","2h79.pdb", "2ito.pdb", "2oxy.pdb","2vag.pdb",
                  "2yj8.pdb", "3v04.pdb", "4e7r.pdb"]
  for pdb_file in X_bonds_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 X_bonds_file = ["5v7d.pdb","2h79.pdb", "2ito.pdb", "2oxy.pdb","2vag.pdb",
             "2yj8.pdb", "3v04.pdb", "4e7r.pdb"]

 i = 0
 for pdb_file in X_bonds_file:
  pdb_file = pdb_file[0:4] + ".updated.pdb"
  pdb_cif = pdb_file[0:4] + ".ligands.cif"
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
        halogen_find_test_No_cif_model(pdb_file)
      else:
        halogen_find_test_cif_model(pdb_file)





if __name__ == '__main__':
    find_the_atoms_makeing_up_halogen_bond_test()