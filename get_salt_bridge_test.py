from __future__ import division
from libtbx import easy_run
import os
from run import hierarchy_No_cif_model
from run import hierarchy_cif_model
from run import get_salt_bridge




def add_H_atoms_into_pad_files():
  SB_bonds_file = ["1ifr.pdb", "1cbr.pdb", "1pga.pdb","1eq5.pdb", "2qmt.pdb", "1mio.pdb", "2on8.pdb", "2onq.pdb"]
  for pdb_file in SB_bonds_file:
      add_h_pdb_file = pdb_file[0:4] + 'h.pdb'
      easy_run.call("phenix.reduce %s > %s " % (pdb_file,  add_h_pdb_file))


# prepare the cif file if the pdb file needs,
#but here it isn't sucessful for 1mio.pdb
def prepare_cif_for_pdb_file():
  SB_bonds_file = ["1ifrh.pdb","1cbrh.pdb","1pgah.pdb","1eq5.pdb","2qmth.pdb","1mioh.pdb","2on8h.pdb","2onqh.pdb"]
  for pdb_file in SB_bonds_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 SB_bonds_file = ["1ifrh.pdb","1cbrh.pdb","1pgah.pdb","1eq5.pdb","2qmth.pdb","1mioh.pdb","2on8h.pdb","2onqh.pdb"]

 i = 0
 for pdb_file in SB_bonds_file:
  pdb_update_file = pdb_file[0:4] + ".updated.pdb"
  pdb_cif = pdb_file[0:4] + ".ligands.cif"
  if os.path.exists(pdb_cif):
    SB_bonds_file[i] = [pdb_update_file, pdb_cif]
  else:
    SB_bonds_file[i] = [pdb_file,None]
  i = i + 1
 return  SB_bonds_file

def find_the_atoms_makeing_up_halogen_bond_test():
    SB_bonds_file = list_cif_and_pdb_file()
    print SB_bonds_file
    for i in range(len(SB_bonds_file)) :
      pdb_file = SB_bonds_file[i][0]
      pdb_cif  = SB_bonds_file[i][1]
      print pdb_file,pdb_cif,"#"*500
      if pdb_cif  == None:
        hierarchy, vdwr=hierarchy_No_cif_model(pdb_file)
        get_salt_bridge(hierarchy,vdwr,eps = 0.3)
        print (pdb_file[0:4] + " this pdb_file is ok")
      else:
        hierarchy, vdwr = hierarchy_cif_model(pdb_file)
        get_salt_bridge(hierarchy,vdwr,eps = 0.3)
        print (pdb_file[0:4] + " this pdb_file_cif is ok")



if __name__ == '__main__':
   # add_H_atoms_into_pad_files()
   #prepare_cif_for_pdb_file()
   find_the_atoms_makeing_up_halogen_bond_test()