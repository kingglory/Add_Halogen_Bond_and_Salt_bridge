from __future__ import division
import iotbx.pdb
from iotbx.pdb import hierarchy
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
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
  halogen_bond_pairs_atom = ["S", 'O', "N"]
  for atom_1 in hierarchy.atoms():
    e1 = atom_1.name.strip().upper()
    if (e1 in halogens):
      (copy_ID_1,resname_1,resid_1) = find_atom_information(atom_1)
      for atom_2 in hierarchy.atoms():
        e2 = atom_2.name.strip().upper()
        if (e2 in halogen_bond_pairs_atom):
          (copy_ID_2, resname_2, resid_2) = find_atom_information(atom_2)
          result = get_rid_of_no_bonding_situations(copy_ID_1,copy_ID_2,resid_1,resid_2)
          if  ((result) ==1 ):
            d = atom_1.distance(atom_2)
            sum_vdwr = vdwr[e1] + vdwr[e2]
            sum_vdwr_min = sum_vdwr*0.6
            if (sum_vdwr_min < d < sum_vdwr):
              print  "hello"
              return (atom_1,atom_2, copy_ID_1, resname_1, resid_1,
                      copy_ID_2, resname_2, resid_2)

# define a function trying to find the third atoms that can make up the angles
#the third atoms make up covalent bond with the knowed atoms
#the distance of covalent bond is near from 1 angstrom to 3 angstrom!!!
#one of the angle (that the halogen atom is at side)is near 120 degrees
# another angle (that the halogen atom is in the middle)is near 180 degrees

def find_the_atoms_makeing_up_halogen_bond(hierarchy,vdwr):
    (atom_1,atom_2,
     copy_ID_1, resname_1, resid_1,
     copy_ID_2, resname_2, resid_2) = get_halogen_bond_pairs(
                                          hierarchy=model.get_hierarchy(),
                                          vdwr=vdwr)
    for atom_3 in hierarchy.atoms():
      (copy_ID_3, resname_3, resid_3) = find_atom_information(atom_3)
      if (copy_ID_3 == copy_ID_1):
        if (resid_3 == resid_1 ):
          if (1< atom_1.distance(atom_3) <3):
            angle_1 = (atom_1.angle(atom_2,atom_3,deg = True))
            if (140 < angle_1 < 180):
              for atom_4 in hierarchy.atoms():
                (copy_ID_4, resname_4, resid_4) = find_atom_information(atom_4)
                if (copy_ID_4 == copy_ID_2):
                  if ( resid_4 == resid_2):
                    if (1 < atom_2.distance(atom_4) < 3):
                      angle_2 = (atom_2.angle(atom_1,atom_4,deg = True))
                      if (90 < angle_2 < 140):
                        print "find halogen bond," \
                              "the information of the four atoms is :"
                        print atom_1.id_str(),atom_2.id_str(),\
                            atom_3.id_str(),atom_4.id_str()



# Second step,find salt bridge in one pdb filess'3
def creat_new_filename(pdb_file):
  new_pdb_file = pdb_file[0:4] + 'h.pdb'
  return new_pdb_file


def add_H_atoms_into_pad_files(pdb_file):
  new_pdb_file = creat_new_filename(pdb_file)
  os.system("phenix.reduce %s > %s " % (pdb_file,  new_pdb_file))


def get_salt_bridge_atom_pairs(hierarchy):
  positive_acide = ["ARG" , "HIS", "LYS"]
  negative_acids = ["ASP" , "GLU"]
  for atom_1 in hierarchy.atoms():
   for atom_3 in hierarchy.atoms():
    if (atom_1.parent().resname in positive_acide):
     if (atom_3.parent().resname in positive_acide):
      (copy_ID_1, resname_1, resid_1) = find_atom_information(atom_1)
      (copy_ID_3, resname_3, resid_3) = find_atom_information(atom_3)
      e1 = atom_1.name.strip().upper()
      e3 = atom_3.name.strip().upper()
      if e1 == "N" :
       if e3 == "H":
        if (0.95 < atom_1.distance(atom_3) <1.05):
         for atom_2 in hierarchy.atoms():
          for atom_4 in hierarchy.atoms():
           if (atom_2.parent().resname in negative_acids):
            if (atom_4.parent().resname in negative_acids):
             (copy_ID_2, resname_2, resid_2) = find_atom_information(atom_2)
             (copy_ID_4, resname_4, resid_4) = find_atom_information(atom_4)
             e2 = atom_2.name.strip().upper()
             e4 = atom_2.name.strip().upper()
             if e2 == "O" :
              if e4 == "C" :
               if(1.3 < atom_1.distance(atom_2) < 2.3):
                if (1.3 < atom_2.distance(atom_3) <2.4):
                 if (170 < atom_2.angle(atom_3 ,atom_4) ):
                  print atom_1.id_str(),atom_2.id_str(),atom_3.id_str(),atom_4.id_str()




























if __name__ == '__main__':
    pdb_file = "5v7d.pdb"
    pdb_inp = iotbx.pdb.input(file_name=pdb_file)
    model = mmtbx.model.manager(
        model_input=pdb_inp,
        process_input=True,
        log=null_out())
    vdwr = model.get_vdw_radii()
    #print vdwr.keys()
    get_halogen_bond_pairs(
       hierarchy=model.get_hierarchy(),
       vdwr=vdwr)

    find_the_atoms_makeing_up_halogen_bond(
       hierarchy=model.get_hierarchy(),
       vdwr=vdwr)

    add_H_atoms_into_pad_files(pdb_file)

    get_salt_bridge_atom_pairs(
        hierarchy=model.get_hierarchy())