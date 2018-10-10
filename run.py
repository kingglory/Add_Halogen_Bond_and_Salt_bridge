from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex
import math
import mmtbx.model
from libtbx.utils import null_out


def distance(s1, s2):
    return math.sqrt(
        (s1[0] - s2[0]) ** 2 +
        (s1[1] - s2[1]) ** 2 +
        (s1[2] - s2[2]) ** 2)

def get_rid_of_no_bonding_situations(atom1,atom2):
    if (len(atom1[1]) == 4):
        if (len(atom2[1]) == 4):
            if (((atom1[1])[0]) != ((atom2[1])[0])):
                if (((atom1[1])[1:4]) == ((atom2[1])[1:4])):
                    return 0
    elif(len(atom1[1])!=len(atom2[1])):
        return 1
    else:
        return 1

def get_halogen_bond_pairs(hierarchy, vdwr):
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ["O","N"]
  for atom_1 in hierarchy.atoms():
    e1 = atom_1.element.strip().upper()
    if (e1 in halogens):
      #print e1, vdwr[e1]
      for atom_2 in hierarchy.atoms():
        e2 = atom_2.element.strip().upper()
        if (e2 in halogen_bond_pairs_atom):
          #print atom_1.id_str()
          ai1 = atom_1.id_str().replace("pdb=","")
          ai2 = atom_2.id_str().replace("pdb=","")
          ai1 = ai1.replace('\"',"")
          ai2 = ai2.replace('\"',"")
          list_ai1 = ai1.split(' ')
          list_ai2 = ai2.split(' ')
          list_ai1_no_space = [item for item in filter(lambda x: x != '', list_ai1)]
          list_ai2_no_space = [item for item in filter(lambda x: x != '', list_ai2)]
          #print  list_ai1_no_space[1],list_ai2_no_space[1]
          # d = distance(atom_1.xyz, atom_2.xyz)
          (result1) = (get_rid_of_no_bonding_situations(list_ai1_no_space, list_ai2_no_space))
          if  ((result1)!=(0)):
            d = atom_1.distance(atom_2)
            sum_vdwr = vdwr[e1] + vdwr[e2]
            sum_vdwr_min = sum_vdwr*0.6
            if (sum_vdwr_min<d<sum_vdwr):
             #ai1 = atom_1.id_str().replace("pdb=","")
             #ai2 = atom_2.id_str().replace("pdb=","")
              print (d,e1,e2,ai1,ai2)
             # a = atom_1.angle(atom_1=b, atom_3=b, deg=True)
            #print d


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


