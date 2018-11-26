from __future__ import division
import iotbx.pdb
import iotbx.cif
from iotbx.pdb import hierarchy
import mmtbx.model
from libtbx.utils import null_out
from libtbx import easy_run
# first step write codes to find halogen bond in one pdb file
# define a function to try finding the halogen bond pairs
def get_halogen_bond_pairs(hierarchy, vdwr,eps = 0.3):
  halogens = ["CL", "BR", "I", "F"]
  # less pi in  halogen_bond_pairs_atom!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  halogen_bond_pairs_atom = ["S","O", "N"]
  result = []
  for atom_1 in hierarchy.atoms():
    atom_e1 = filter(str.isalpha,atom_1.element.upper() )
    e1 = atom_1.name.strip().upper()
    if (atom_e1 in halogens):
      for atom_2 in hierarchy.atoms():
       atom_e2 = filter(str.isalpha,atom_2.element.upper() )
       e2 = atom_2.name.strip().upper()
       if (not atom_1.is_in_same_conformer_as(atom_2)): continue
       if (atom_1.parent().parent().resseq == atom_2.parent().parent().resseq): continue
       if (atom_e2[0] in halogen_bond_pairs_atom):
        # O2' in 3v04.pdb file will recognized as O2* ,so replace it
        e1 = e1.replace("'","*")
        e2 = e2.replace("'","*")
        # 2yj8.pdb  vdwr can't recognize 'OXT'
        e2 = e2.replace("XT","")
        d = atom_1.distance(atom_2)
        sum_vdwr = vdwr[e1] + vdwr[e2]
        sum_vdwr_min = sum_vdwr*0.6
        if (sum_vdwr_min-eps < d < sum_vdwr+eps):
          d_x_p = d/sum_vdwr
          result.append([d,sum_vdwr,d_x_p,atom_1,atom_2])
  return result
# test。。。。。。。。。。

"""define a function trying to find the third atoms that can make up the angles
   the third atoms make up covalent bond with the knowed atoms
   the distance of covalent bond is near from 1 angstrom to 3 angstrom!!!
   one of the angle (that the halogen atom is at side)is near 120 degrees
   another angle (that the halogen atom is in the middle)is near 180 degrees
"""

def find_the_atoms_makeing_up_halogen_bond(hierarchy,vdwr):
  result = get_halogen_bond_pairs(hierarchy, vdwr,eps = 0.3)
  d_list = []
  info_result = []
  if (result is not None):
   for i in range(len(result)):
    (d,sum_vdwr,d_x_p,atom_1,atom_2) = result[i]
    for atom_3 in hierarchy.atoms():
     e3 = filter(str.isalpha, atom_3.element.upper())
     if e3[0] == "C":
      e1 = atom_1.name.strip().upper()
      sum_vdwr1 = vdwr[e1] + vdwr[e3]
      if (not atom_1.is_in_same_conformer_as(atom_3)): continue
      if (atom_1.parent().parent().resseq==atom_3.parent().parent().resseq):
       if (1.3< atom_1.distance(atom_3) <3):
        angle_1 = (atom_1.angle(atom_2,atom_3,deg = True))
        if (140 < angle_1 ):
         for atom_4 in hierarchy.atoms():
          e4 = filter(str.isalpha, atom_4.element.upper())
          if (not atom_2.is_in_same_conformer_as(atom_4)): continue
          if (atom_2.parent().parent().resseq==atom_4.parent().parent().resseq):
           if e4[0] == "C":
            if (atom_3.distance(atom_2)<1.9):continue
            if (atom_4.distance(atom_1) < sum_vdwr1): continue
            if (1 < atom_2.distance(atom_4) < 3):
             angle_2 = (atom_2.angle(atom_1,atom_4,deg = True))
             if (90 < angle_2 < 160):
              d_list.append(result[i][0])
              info_result.append((atom_3.id_str(),atom_1.id_str(),
                                  atom_2.id_str(), atom_4.id_str(),
                                  d,sum_vdwr,d_x_p,angle_1,angle_2
                                  ))
   d_min = min(d_list)
   print (d_min)
   for i in range(len(result)):
    if (info_result[i][4])==d_min:
     print ("find halogen bond," "the information of the four atoms is :")
     print ("**" * 50)
     print (info_result[i])
     print info_result[i][4]



def hierarchy_cif_model(pdb_file):
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
  return hierarchy, vdwr


def hierarchy_No_cif_model(pdb_file):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  model = mmtbx.model.manager(model_input=pdb_inp,
                                  process_input=True,
                                  log=null_out())
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  return hierarchy,vdwr


#Second step,find salt bridge in one pdb filess

Amino_Acids = ["ARG","HIS","LYS","ASP","GLU","SER","THR","ASN","GLU","CYS","SEC",
               "GLY","PRO","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP"]


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
def get_salt_bridge(hierarchy,vdwr,eps = 0.3):
  positive_acide = ["ARG" , "HIS", "LYS"]
  negative_acids = ["ASP" , "GLU"]
  for atom_1 in hierarchy.atoms():
   for atom_3 in hierarchy.atoms():
    if (atom_1.parent().resname in positive_acide):
     if (atom_3.parent().resname in positive_acide):
      if (not atom_1.is_in_same_conformer_as(atom_3)): continue
      if(atom_1.parent().parent().resseq==atom_3.parent().parent().resseq) :
        e1 = filter(str.isalpha,atom_1.name.upper() )
        e3 = filter(str.isalpha,atom_3.name.upper() )
        if e1[0] == "N" :
         if e3[0] == "H":
          Length_n_h =1.01
          if ( Length_n_h - eps < atom_1.distance(atom_3) < Length_n_h + eps):
           for atom_2 in hierarchy.atoms():
            if (atom_2.parent().resname in negative_acids):
             if (not atom_1.is_in_same_conformer_as(atom_2)): continue
             if (atom_1.parent().parent().resseq != atom_2.parent().parent().resseq):
              e2 = filter(str.isalpha, atom_2.name.upper())
              dict_h_bond_lengh = { "O" : 0.98,"N" : 1.45,"F" : 0.92}
              dict_ionic_bond_lengh={}
              if e2[0] in dict_h_bond_lengh.keys():
               Length_o_n = 1.46
               h_d = dict_h_bond_lengh[e2[0]]
               sum_vdwr = vdwr[e2[0]] + vdwr[e3]
               if(Length_o_n - 0.1 < atom_1.distance(atom_2) < 4):
                 if (h_d-eps< atom_2.distance(atom_3) < sum_vdwr ):
                  angle_1 = (atom_3.angle(atom_1, atom_2, deg=True))
                  if(90< angle_1):
                   print ("find the salt bridge ")
                   print (atom_1.id_str(),atom_2.id_str(),
                             atom_3.id_str())




""" path = "/home/pdb/mirror/pub/pdb/data/structures/divided/pdb"
    of  = open("".join([path,"INDEX"]),"r")
    files= ["".join([path,f]).strip() for f in of.readlines()]
    of.close()
    list = []
    for f in files:
      pdb_code = os.path.basename(f)[3:7]
      pdb_file = str(pdb_code) + ".pdb"



    add_H_atoms_into_pad_files(pdb_file)

    get_salt_bridge_atom_pairs(
       hierarchy=model.get_hierarchy())
"""