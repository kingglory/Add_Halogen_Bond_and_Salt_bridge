from __future__ import division
from libtbx import easy_run
import os
from run import hierarchy_No_cif_model
from run import hierarchy_cif_model
from run import get_salt_bridge


# prepare the cif file if the pdb file needs
def prepare_cif_for_pdb_file():
  X_bonds_file = []
  for pdb_file in X_bonds_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 X_bonds_file = []

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
        hierarchy, vdwr=hierarchy_No_cif_model(pdb_file)
        get_salt_bridge(hierarchy)
        print (pdb_file[0:4] + " this pdb_file is ok")
      else:
        hierarchy, vdwr = hierarchy_cif_model(pdb_file)
        get_salt_bridge(hierarchy)
        print (pdb_file[0:4] + " this pdb_file_cif is ok")
