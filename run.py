from __future__ import division
import iotbx.pdb
import iotbx.cif
from iotbx.pdb import hierarchy
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
from libtbx import easy_run
import os
# first step write codes to find halogen bond in one pdb file

# because the unit and it's copy can't make halogen bond toghter,
# so here is a function that can get rid of this situation
def get_rid_of_no_bonding_situations(copy_ID_1,copy_ID_2,resid_1,resid_2):
    if copy_ID_1  == copy_ID_2 :
      if resid_1 != resid_2 :
        return 1
    if copy_ID_1 != copy_ID_2:
      if copy_ID_1 == " ":
        return 1
      if copy_ID_2 == " ":
        return 1

def find_atom_information(atom):
    return ( (atom.id_str()[9]),
             (atom.parent().resname),
             (atom.parent().parent().resid())
             )
    #print (atom.i_seq)

# define a function to try finding the halogen bond pairs
def get_halogen_bond_pairs(hierarchy, vdwr):
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ['S','O', "N","CL","BR","I"]
  #print vdwr.keys()
  for atom_1 in hierarchy.atoms():
    atom_e1 = filter(str.isalpha,atom_1.name.upper() )
    e1 = atom_1.name.strip().upper()
    if (atom_e1 in halogens):
      print (e1+"find"+"*"*50)
      (copy_ID_1,resname_1,resid_1) = find_atom_information(atom_1)
      for atom_2 in hierarchy.atoms():
        atom_e2 = filter(str.isalpha,atom_2.name.upper() )
        e2 = atom_2.name.strip().upper()
        if (atom_e2 in halogen_bond_pairs_atom):
          (copy_ID_2, resname_2, resid_2) = find_atom_information(atom_2)
          result = get_rid_of_no_bonding_situations(copy_ID_1,copy_ID_2,resid_1,resid_2)
          if  ((result) ==1 ):
            d = atom_1.distance(atom_2)
            # O2' in pdb file will recognized as O2* ,so replace it
            e1 = e1.replace("'","*")
            e2 = e2.replace("'","*")
            sum_vdwr = vdwr[e1] + vdwr[e2]
            sum_vdwr_min = sum_vdwr*0.6
            if (sum_vdwr_min < d < sum_vdwr):
              d_x_p = d/sum_vdwr

              return (d_x_p,atom_1,atom_2, copy_ID_1, resname_1, resid_1,
                      copy_ID_2, resname_2, resid_2)

# define a function trying to find the third atoms that can make up the angles
#the third atoms make up covalent bond with the knowed atoms
#the distance of covalent bond is near from 1 angstrom to 3 angstrom!!!
#one of the angle (that the halogen atom is at side)is near 120 degrees
# another angle (that the halogen atom is in the middle)is near 180 degrees

def find_the_atoms_makeing_up_halogen_bond(hierarchy,vdwr):
  result = get_halogen_bond_pairs(hierarchy,vdwr)
  if (result is not None):
   (d_x_p,atom_1, atom_2,copy_ID_1, resname_1, resid_1,
       copy_ID_2, resname_2, resid_2) = result
   for atom_3 in hierarchy.atoms():
    (copy_ID_3, resname_3, resid_3) = find_atom_information(atom_3)
    if (copy_ID_3 == copy_ID_1):
     if (resid_3 == resid_1 ):
      if (1.3< atom_1.distance(atom_3) <3):
        angle_1 = (atom_1.angle(atom_2,atom_3,deg = True))
        if (140 < angle_1 ):

         for atom_4 in hierarchy.atoms():
          (copy_ID_4, resname_4, resid_4) = find_atom_information(atom_4)
          if (copy_ID_4 == copy_ID_2):
            if ( resid_4 == resid_2):
             if (1.3 < atom_2.distance(atom_4) < 3):
              angle_2 = (atom_2.angle(atom_1,atom_4,deg = True))
              if (90 < angle_2 < 140):
                print ("**********************************************************************************************")
                print ("find halogen bond," 
                        "the information of the four atoms is :")
                print (atom_1.id_str(),atom_2.id_str(),atom_1.distance(atom_2),
                        angle_1,d_x_p,atom_3.id_str(),atom_4.id_str())
                print ("**********************************************************************************************")

Amino_Acids = ["ARG","HIS","LYS","ASP","GLU","SER","THR","ASN","GLU","CYS","SEC",
               "GLY","PRO","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP"]



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





def halogen_find_test_cif_model(pdb_file):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file,
                            source_info=None)
  pdb_cif = pdb_file[0:4] + ".ligands.cif"
  cif_object = iotbx.cif.reader(pdb_cif).model()
  cif_objects = [(pdb_cif, cif_object)]
  model = mmtbx.model.manager(model_input=pdb_inp,
                              build_grm=True,
                              restraint_objects=cif_objects,
                              log=null_out())
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  print (pdb_file[0:4] +" this pdb_file_cif is ok")
  find_the_atoms_makeing_up_halogen_bond(hierarchy,vdwr)
  print (pdb_file[0:4] +" this pdb_file_cif is ok")



def halogen_find_test_No_cif_model(pdb_file):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  model = mmtbx.model.manager(model_input=pdb_inp,
                                  process_input=True,
                                  log=null_out())
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  print (pdb_file[0:4] +" this pdb_file is ok")
  find_the_atoms_makeing_up_halogen_bond(hierarchy,vdwr)
  print (pdb_file[0:4] +" this pdb_file is ok")


def find_the_atoms_makeing_up_halogen_bond_test():
    X_bonds_file = list_cif_and_pdb_file()
    for i in range(len(X_bonds_file)) :
      pdb_file = X_bonds_file[i][0]
      if X_bonds_file[i][1] == None:
        halogen_find_test_No_cif_model(pdb_file)
      else:
        halogen_find_test_cif_model(pdb_file)

#Second step,find salt bridge in one pdb filess

def creat_new_filename(pdb_file):
  new_pdb_file = pdb_file[0:4] + 'h.pdb'
  return new_pdb_file


def add_H_atoms_into_pad_files(pdb_file):
  new_pdb_file = creat_new_filename(pdb_file)
  easy_run.call("phenix.reduce %s > %s " % (pdb_file,  new_pdb_file))

# in the twenty one Amino Acids,[arg,his,lys] are positive Amino Acids
#  with electrically charged side chains,the positive charged atom is N
# there are some H atoms around this N atom
# [asp,glu] are negative,the negative charged atom is O
# the N atom and the O atom will make up the iron bond,the o atom and
# one of the H atom make up the H bond
def get_salt_bridge_atom_pairs(hierarchy):
  positive_acide = ["ARG" , "HIS", "LYS"]
  negative_acids = ["ASP" , "GLU"]
  for atom_1 in hierarchy.atoms():
   for atom_3 in hierarchy.atoms():
    if (atom_1.parent().resname in positive_acide):
     if (atom_3.parent().resname in positive_acide):
      (copy_ID_1, resname_1, resid_1) = find_atom_information(atom_1)
      (copy_ID_3, resname_3, resid_3) = find_atom_information(atom_3)
      if copy_ID_1 == copy_ID_3 :
       if resid_1 == resid_3 :
        e1 = filter(str.isalpha,atom_1.name.upper() )
        e3 = filter(str.isalpha,atom_3.name.upper() )
        if e1 == "N" :
         if e3 == "H":
          if (0.95 < atom_1.distance(atom_3) <1.05):
           for atom_2 in hierarchy.atoms():
            for atom_4 in hierarchy.atoms():
             if (atom_2.parent().resname in negative_acids):
              if (atom_4.parent().resname in negative_acids):
               (copy_ID_2, resname_2, resid_2) = find_atom_information(atom_2)
               (copy_ID_4, resname_4, resid_4) = find_atom_information(atom_4)
               if copy_ID_2 == copy_ID_4 :
                if resid_2 == resid_4 :
                 e2 = filter(str.isalpha,atom_2.name.upper() )
                 e4 = filter(str.isalpha,atom_4.name.upper() )
                 if e2 == "O" :
                  if e4 == "C" :
                   if(1.3 < atom_1.distance(atom_2) < 2.3):
                    if (1.3 < atom_2.distance(atom_3) <2.4):
                     if (170 < atom_2.angle(atom_3 ,atom_4) ):
                      print ("find the salt bridge ")
                      print (atom_1.id_str(),atom_2.id_str(),
                             atom_3.id_str(),atom_4.id_str())



if __name__ == '__main__':
 #   path = "/home/pdb/mirror/pub/pdb/data/structures/divided/pdb"
 #   of   = open("".join([path,"INDEX"]),"r")
 #   files= ["".join([path,f]).strip() for f in of.readlines()]
 #   of.close()
 #   list = []
 #   for f in files:
 #     pdb_code = os.path.basename(f)[3:7]
 #     pdb_file = str(pdb_code) + ".pdb"


 """   pdb_file = "5v7d.pdb"
    pdb_inp = iotbx.pdb.input(file_name=pdb_file)
    model = mmtbx.model.manager(
        model_input=pdb_inp,
        process_input=True,
       log=null_out())
    vdwr = model.get_vdw_radii()
    print vdwr.keys()
    get_halogen_bond_pairs(
       hierarchy=model.get_hierarchy(),
       vdwr=vdwr)

    find_the_atoms_makeing_up_halogen_bond(
       hierarchy=model.get_hierarchy(),
       vdwr=vdwr)
    prepare_cif_for_pdb_file()
    
       """

 find_the_atoms_makeing_up_halogen_bond_test()

#    add_H_atoms_into_pad_files(pdb_file)

#    get_salt_bridge_atom_pairs(
#       hierarchy=model.get_hierarchy())