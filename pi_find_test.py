from __future__ import division
from libtbx import easy_run
import os
from run import hierarchy_No_cif_model
from run import hierarchy_cif_model
from run import define_pi_system
# prepare the cif file if the pdb file needs
def prepare_cif_for_pdb_file():
  pi_file = ["1c14.pdb","3Az9.pdb"]
  for pdb_file in pi_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 pi_file = ["1c14.pdb","3Az9.pdb"]
 i = 0
 for pdb_file in pi_file:
  pdb_file = pdb_file[0:4] + ".updated.pdb"
  pdb_cif = pdb_file[0:4] + ".ligands.cif"
  if os.path.exists(pdb_cif):
    pi_file[i] = [pdb_file, pdb_cif]
  else:
    pi_file[i] = [pdb_file,None]
  i = i + 1
 return  pi_file

def find_pi_test():
    pi_file = list_cif_and_pdb_file()
    for i in range(len(pi_file)) :
      pdb_file = pi_file[i][0]
      if pi_file[i][1] == None:
        hierarchy, vdwr,model=hierarchy_No_cif_model(pdb_file)
        define_pi_system(hierarchy, vdwr,model)
        print (pdb_file[0:4] + " this pdb_file is ok")
      else:
        hierarchy, vdwr,model = hierarchy_cif_model(pdb_file)
        define_pi_system(hierarchy, vdwr,model)
        print (pdb_file[0:4] + " this pdb_file_cif is ok")

if __name__ == '__main__':
 #   prepare_cif_for_pdb_file()
    find_pi_test()