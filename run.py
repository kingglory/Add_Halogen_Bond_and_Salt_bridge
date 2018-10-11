from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex
import math
import mmtbx.model
from libtbx.utils import null_out
import os
# first step write codes to find halogen bond in one pdb file
def distance(s1, s2):
    return math.sqrt(
        (s1[0] - s2[0]) ** 2 +
        (s1[1] - s2[1]) ** 2 +
        (s1[2] - s2[2]) ** 2)
# because the unit and it's copy can't make halogen bond toghter,
# so here is a function that can get rid of this situation
def get_rid_of_no_bonding_situations(atom1,atom2):
    if (len(atom1[1]) == 4):
        if (len(atom2[1]) == 4):
            if (((atom1[1])[0]) != ((atom2[1])[0])):
                #f (((atom1[1])[1:4]) == ((atom2[1])[1:4])):
              return 0
    #if it can't make the halogen bond between the same residue,
            elif((atom1[2:4]==atom2[2:4])and(((atom1[1])[0]) == ((atom2[1])[0]))):
              return 0
    if (len(atom1[1])==3):
      if(len(atom2[1])==3):
         if(atom1[1:4]==atom2[1:4]):
             return 0
    # The cut-off line for elif statement above
    elif(len(atom1[1])!=len(atom2[1])):
        return 1
    else:
        return 1


# define a function that can change from the atom got to a list of atom.id_str()
def atom_to_list_of_atom_id_str(atom):
    ai = atom.id_str().replace("pdb=","")
    ai = ai.replace('\"',"")
    list_ai = ai.split(' ')
    list_ai_no_space = [item for item in filter(lambda x: x != '', list_ai)]
    #print list_ai_no_space
    return list_ai_no_space





# define a function to try finding the halogen bond pairs
def get_halogen_bond_pairs(hierarchy, vdwr):
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ["s", 'O', "N"]
  for atom_1 in hierarchy.atoms():
    e1 = atom_1.element.strip().upper()
    if (e1 in halogens):
      for atom_2 in hierarchy.atoms():
        e2 = atom_2.element.strip().upper()
        if (e2 in halogen_bond_pairs_atom):
          #print atom_1.id_str()
          #ai1 = atom_1.id_str().replace("pdb=","")
          #ai2 = atom_2.id_str().replace("pdb=","")
          #ai1 = ai1.replace('\"',"")
          #ai2 = ai2.replace('\"',"")
          #list_ai1 = ai1.split(' ')
          #list_ai2 = ai2.split(' ')
          #list_ai1_no_space = [item for item in filter(lambda x: x != '', list_ai1)]
          #list_ai2_no_space = [item for item in filter(lambda x: x != '', list_ai2)]
          #print  list_ai1_no_space[1],list_ai2_no_space[1]
          # d = distance(atom_1.xyz, atom_2.xyz)
          list_ai1_no_space = atom_to_list_of_atom_id_str(atom_1)
          list_ai2_no_space = atom_to_list_of_atom_id_str(atom_2)
          (result1) = (get_rid_of_no_bonding_situations(list_ai1_no_space, list_ai2_no_space))
          if  ((result1)!=(0)):
            d = atom_1.distance(atom_2)
            sum_vdwr = vdwr[e1] + vdwr[e2]
            sum_vdwr_min = sum_vdwr*0.6
            if (sum_vdwr_min< d <sum_vdwr):
              #print (d, e1, e2,list_ai1_no_space,list_ai2_no_space,vdwr[e2],vdwr[e1])
              return (atom_1,atom_2)
             #ai1 = atom_1.id_str().replace("pdb=","")
             #ai2 = atom_2.id_str().replace("pdb=","")
             # a = atom_1.angle(atom_1=b, atom_3=b, deg=True)
            #print d
# define a function trying to find the third atoms that can make up the angles
#the third atoms make up covalent bond with the knowed atoms
#the distance of covalent bond is near from 1 angstrom to 3 angstrom!!!
#one of the angle (that the halogen atom is at side)is near 120 degrees
# another angle (that the halogen atom is in the middle)is near 180 degrees

def find_the_atoms_makeing_up_halogen_bond(hierarchy):
    (atom_1,atom_2) = get_halogen_bond_pairs(
                          hierarchy=model.get_hierarchy(),
                          vdwr=vdwr)
    list_ai1_no_space = atom_to_list_of_atom_id_str(atom_1)
    list_ai2_no_space = atom_to_list_of_atom_id_str(atom_2)
    for atom_3 in hierarchy.atoms():
      list_ai3_no_space = atom_to_list_of_atom_id_str(atom_3)
      if (list_ai3_no_space[2:4]==list_ai1_no_space[2:4]):
        if (list_ai3_no_space[2:4]!=list_ai2_no_space[2:4]):
          if (1< atom_1.distance(atom_3) <3):
            angle_1 = (atom_1.angle(atom_2,atom_3,deg = True))
            if (90 < angle_1 < 180):
              for atom_4 in hierarchy.atoms():
                list_ai4_no_space = atom_to_list_of_atom_id_str(atom_4)
                if (list_ai4_no_space[2:4] == list_ai2_no_space[2:4]):
                  if (list_ai4_no_space[2:4] != list_ai1_no_space[2:4]):
                    if (1 < atom_2.distance(atom_4) < 3):
                      angle_2 = (atom_2.angle(atom_1,atom_4,deg = True))
                      if (90 < angle_2 < 140):
                        print "find halogen bond,the information of the four atoms is :"
                        print (list_ai3_no_space ,list_ai1_no_space,list_ai2_no_space,list_ai4_no_space)
                        return list_ai3_no_space ,list_ai1_no_space,list_ai2_no_space,list_ai4_no_space
                  #print (atom_1,atom_2,atom_3,atom_4)





# Second step,find salt bridge in one pdb file
def add_H_atoms_into_pad_files(pdbfile):
    os.system("phenix.reduce %s > %s " % (pdb_file,+ 'h' + '.pdb'))


























if __name__ == '__main__':
    pdb_file = "5v7d.pdb"
    pdb_inp = iotbx.pdb.input(file_name=pdb_file)
    model = mmtbx.model.manager(
        model_input=pdb_inp,
        process_input=True,
        log=null_out())
    vdwr = model.get_vdw_radii()
    #print vdwr.keys()
    #get_halogen_bond_pairs(
    #   hierarchy=model.get_hierarchy(),
    #   vdwr=vdwr)

    find_the_atoms_makeing_up_halogen_bond(
        hierarchy=model.get_hierarchy())
