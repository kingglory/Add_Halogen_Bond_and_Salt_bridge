from __future__ import division
import iotbx.pdb
import iotbx.cif
from iotbx.pdb import hierarchy
import mmtbx.model
from libtbx.utils import null_out
from libtbx import easy_run
from libtbx import group_args

def ideal_distance(atom_1,atom_2,model):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
        sites_cart=model.get_sites_cart())
  for proxy in bond_proxies_simple:
    i_seq, j_seq = proxy.i_seqs
    if (atom_1.i_seq in proxy.i_seqs):
     if (atom_1.i_seq == i_seq):
      if (atom_2.i_seq == j_seq):
        d = atom_1.distance(atom_2)
        return atom_1 ,atom_2,d

# first step write codes to find halogen bond in one pdb file
# define a function to try finding the halogen bond pairs
def find_water(hierarchy):
    get_class = iotbx.pdb.common_residue_names_get_class
    for model in hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                for ag in rg.atom_groups():
                    if (get_class(ag.resname)=="water"):
                        print  ("water here :",ag)
                        return ag

def find_halogen_bonds(model, eps = 0.3, emp_scale = 0.6):
  hierarchy = model.get_hierarchy()
  vdwr      = model.get_vdw_radii()
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ["S", "O", "N"] #XXX Add halogens here
  result = []
  for a1 in hierarchy.atoms():
    e1 = a1.element.upper()
    n1 = a1.name.strip().upper()
    if(e1 in halogens):
      for a2 in hierarchy.atoms():
        e2 = a2.element.upper()
        n2 = a2.name.strip().upper()
        if(not a1.is_in_same_conformer_as(a2)): continue
        if(a1.parent().parent().resseq == a2.parent().parent().resseq): continue
        if(e2 in halogen_bond_pairs_atom):
          # O2' in 3v04.pdb file will recognized as O2* ,so replace it
          n1 = n1.replace("'","*")
          n2 = n2.replace("'","*")
          # 2yj8.pdb  vdwr can't recognize 'OXT'
          n2 = n2.replace("XT","")
          d = a1.distance(a2)
          sum_vdwr = vdwr[n1] + vdwr[n2]
          sum_vdwr_min = sum_vdwr*emp_scale
          if(sum_vdwr_min-eps < d < sum_vdwr+eps): # found HB pairs-candidates
            d_x_p = d/sum_vdwr
            for a3 in hierarchy.atoms():
              e3 = a3.element.upper()
              if(e3 in ["C","P","S"]): # XXX Cite paper
                e3 = a3.name.strip().upper() # XXX Fix naming
                e1 = a1.name.strip().upper()
                e3 = e3.replace("'", "*") # XXX Fix naming
                sum_vdwr1 = vdwr[e1] + vdwr[e3] # XXX Fix naming; use pair_proxy instead!
                if(not a1.is_in_same_conformer_as(a3)): continue
                if(a1.parent().parent().resseq==a3.parent().parent().resseq):
                  if(1.3 < a1.distance(a3) < 3): # XXX use pair_proxy instead, don't calculate distance!
                    angle_123 = (a1.angle(a2, a3, deg = True)) # theta_1 angle in PAPER
                    if(140 < angle_123): # See figure XXX in paper XXX
                      for a4 in hierarchy.atoms():
                        e4 = a4.element.upper()
                        if(not a2.is_in_same_conformer_as(a4)): continue
                        if(a2.parent().parent().resseq == 
                           a4.parent().parent().resseq):
                          if(a3.distance(a2) < 1.9): continue # XXX use pair_proxy instead
                          if(a4.distance(a1) < sum_vdwr1): continue # XXX use pair_proxy instead
                          if(1 < a2.distance(a4) < 3): # XXX use pair_proxy instead
                            angle_214 = (a2.angle(a1, a4, deg = True)) # theta_2 angle in PAPER
                            if(90 < angle_214 < 160): # theta_2 angle in PAPER
                              result.append(group_args(
                                atom_1    = a1,
                                atom_2    = a2,
                                atom_3    = a3,
                                atom_4    = a4,
                                d12       = d,
                                sum_vdwr  = sum_vdwr,
                                d_x_p     = d_x_p,
                                angle_123 = angle_123,
                                angle_214 = angle_214))
  return result

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
  return hierarchy, vdwr,model


def hierarchy_No_cif_model(pdb_file):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  model = mmtbx.model.manager(model_input=pdb_inp,
                                  process_input=True,
                                  log=null_out())
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  return hierarchy,vdwr,model


#Second step,find salt bridge in one pdb filess

Amino_Acids = ["ARG","HIS","LYS","ASP","GLU","SER","THR",
               "ASN","GLU","CYS","SEC","GLY","PRO","ALA",
               "VAL","ILE","LEU","MET","PHE","TYR","TRP"]


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
def get_salt_bridge(hierarchy,vdwr,model,eps = 0.3):
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
          if ( Length_n_h - eps <
                   atom_1.distance(atom_3)
                   < Length_n_h + eps):
           geometry = model.get_restraints_manager()
           bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                  sites_cart=model.get_sites_cart())
           for proxy in bond_proxies_simple:
            i_seq, j_seq = proxy.i_seqs
            if (atom_1.i_seq in proxy.i_seqs):
             if (atom_1.i_seq == i_seq):
              if (atom_3.i_seq == j_seq):
               for atom_2 in hierarchy.atoms():
                if (atom_2.parent().resname in negative_acids):
                 if (not atom_1.is_in_same_conformer_as(atom_2)): continue
                 if (atom_1.parent().parent().resseq !=
                         atom_2.parent().parent().resseq):
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

def define_pi_system(hierarchy,vdwr,model,eps = 5):
 pi_amino_acids = ["HIS","PRO","PHE","TYR","TYP"]
 for atom_1 in hierarchy.atoms():
  for atom_2 in hierarchy.atoms():
   for atom_3 in hierarchy.atoms():
    for atom_4 in hierarchy.atoms():
     for atom_5 in hierarchy.atoms():
      for atom_6 in hierarchy.atoms():
       if (atom_1.parent().resname==
           atom_2.parent().resname==
           atom_3.parent().resname==
           atom_4.parent().resname==
           atom_5.parent().resname==
           atom_6.parent().resname):
        if (not atom_1.parent().resname in pi_amino_acids): continue
        if (not atom_1.is_in_same_conformer_as(atom_2)): continue
        if (not atom_1.is_in_same_conformer_as(atom_3)): continue
        if (not atom_1.is_in_same_conformer_as(atom_4)): continue
        if (not atom_1.is_in_same_conformer_as(atom_5)): continue
        if (not atom_1.is_in_same_conformer_as(atom_6)): continue
        result_12 = ideal_distance(atom_1, atom_2, model)
        if result_12 is None:continue
        result_23 = ideal_distance(atom_2, atom_3, model)
        if result_23 is None: continue
        result_34 = ideal_distance(atom_3, atom_4, model)
        if result_34 is None: continue
        result_45 = ideal_distance(atom_4, atom_5, model)
        if result_45 is None: continue
        result_56 = ideal_distance(atom_5, atom_6, model)
        if result_56 is None: continue
        angle1 = atom_1.angle(atom_2 , atom_6, degree = True)
        angle2 = atom_2.angle(atom_1 , atom_3, degree = True)
        angle3 = atom_3.angle(atom_2 , atom_4, degree = True)
        angle4 = atom_4.angle(atom_3 , atom_5, degree = True)
        angle5 = atom_5.angle(atom_4 , atom_6, degree = True)
        angle6 = atom_6.angle(atom_5 , atom_1, degree = True)
        if (abs(angle1 - 120) > eps): continue
        if (abs(angle2 - 120) > eps): continue
        if (abs(angle3 - 120) > eps): continue
        if (abs(angle4 - 120) > eps): continue
        if (abs(angle5 - 120) > eps): continue
        if (abs(angle6 - 120) > eps): continue
        print (atom_1 ,atom_2 ,atom_3 ,atom_4 ,atom_5 ,atom_6 )






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