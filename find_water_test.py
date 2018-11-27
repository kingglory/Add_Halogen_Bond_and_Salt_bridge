from __future__ import division
from libtbx import easy_run
import os
from run import hierarchy_No_cif_model
from run import hierarchy_cif_model
from run import find_water

# didn't find water containing pdb files,so just useing the follow files to try,then didn't find water,work hard!!!

def add_H_atoms_into_pad_files():
  W_bonds_file = ["5v7d.pdb","2h79.pdb", "2ito.pdb", "2oxy.pdb","2vag.pdb"]
  for pdb_file in W_bonds_file:
      add_h_pdb_file = pdb_file[0:4] + 'h.pdb'
      easy_run.call("phenix.reduce %s > %s " % (pdb_file,  add_h_pdb_file))

# prepare the cif file if the pdb file needs,
#but here it isn't sucessful for 1mio.pdb
def prepare_cif_for_pdb_file():
  W_bonds_file = ["5v7dh.pdb","2h79h.pdb", "2itoh.pdb", "2oxyh.pdb","2vagh.pdb"]
  for pdb_file in W_bonds_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 W_bonds_file = ["5v7dh.pdb","2h79h.pdb", "2itoh.pdb", "2oxyh.pdb","2vagh.pdb"]
 i = 0
 for pdb_file in W_bonds_file:
  pdb_update_file = pdb_file[0:5] + ".updated.pdb"
  pdb_cif = pdb_file[0:5] + ".ligands.cif"
  if os.path.exists(pdb_cif):
    W_bonds_file[i] = [pdb_update_file, pdb_cif]
  else:
    W_bonds_file[i] = [pdb_file,None]
  i = i + 1
 return  W_bonds_file

def find_water_test():
    W_bonds_file = list_cif_and_pdb_file()
    print (W_bonds_file)
    for i in range(len(W_bonds_file)) :
      pdb_file = W_bonds_file[i][0]
      pdb_cif  = W_bonds_file[i][1]
      print (pdb_file,pdb_cif,"#"*500)
      if pdb_cif  == None:
        hierarchy, vdwr,model=hierarchy_No_cif_model(pdb_file)
        find_water(hierarchy)
        print (pdb_file[0:4] + " this pdb_file is ok")
      else:
        hierarchy, vdwr,model = hierarchy_cif_model(pdb_file)
        find_water(hierarchy)
        print (pdb_file[0:4] + " this pdb_file_cif is ok")



if __name__ == '__main__':
   # add_H_atoms_into_pad_files()
   #prepare_cif_for_pdb_file()
   find_water_test()