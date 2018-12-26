from __future__ import division
import iotbx.pdb
import iotbx.cif
from libtbx import group_args
from libtbx import easy_run

def is_bonded(atom_1, atom_2, bond_proxies_simple):
  result = False
  for proxy in bond_proxies_simple:
    i_seq, j_seq = proxy.i_seqs
    if(atom_1.i_seq in proxy.i_seqs and atom_2.i_seq in proxy.i_seqs):
      result = True
      break
  return result

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
"""
1,when one halogen atoms make halogen bond when other atom,when the distance and angle
limations all fited,the bond distance more shorter ,more possible;
2, when chooseing the third atoms to make up the angle1,the angle1 is more near 180,more possible;
3,when chooseing the fourth atoms to make up the angle2,the angle1 is more near 120,more possible;
4,see Figure 4 in paper "Halogen bonding(X-bonding):  A biological perspective"
5,theta_2 angle :Geometry of X-bonds in paper "Halogen bond in biological molecules"
"""

def find_atom_3(d, sum_vdwr, hierarchy, a1, a2, bond_proxies_simple):
    """
    this function is help to short and complete the find X or H bonds function
    :param d: the distance between a1 and a2
    :param sum_vdwr: the sum of (a1 and a2)'s vdwr
    :param bond_proxies_simple: a object that will help to define if it will make covalent bond between
    :param result_123: the informations of a1 and a2 an a3
    """
    result_123 = []
    d_x_p = d / sum_vdwr
    for a3 in hierarchy.atoms():
      if (not is_bonded(a1, a3, bond_proxies_simple)): continue
      # theta_1 angle in paper "Halogen bond in biological molecules"
      angle_312 = (a1.angle(a2, a3, deg=True))
      # See Fig.1 in paper "Halogen bond in biological molecules"
      if (130 < angle_312):
        result_123.append(group_args(
                a1        = a1,
                a2        = a2,
                a3        = a3,
                d_12      = d,
                d_x_p     = d_x_p,
                sum_vdwr  = sum_vdwr,
                angle_312 = angle_312
            ))
    if result_123 is None:
      for a3 in hierarchy.atoms():
        if (not is_bonded(a1, a3, bond_proxies_simple)): continue
        # theta_1 angle in paper "Halogen bond in biological molecules"
        angle_312 = (a1.angle(a2, a3, deg=True))
        # See Fig.1 in paper "Halogen bond in biological molecules"
        if (90 < angle_312):
          result_123.append(group_args(
                    a1        = a1,
                    a2        = a2,
                    a3        = a3,
                    d_12      = d,
                    d_x_p     = d_x_p,
                    sum_vdwr  = sum_vdwr,
                    angle_312 = angle_312
                ))
    return result_123

def find_atom4(hierarchy,result_123, bond_proxies_simple, angle_eps):
  result = []
  result4 = None
  for r_123 in result_123 :
    a1        = r_123.a1
    a2        = r_123.a2
    a3        = r_123.a3
    d         = r_123.d_12
    d_x_p     = r_123.d_x_p
    sum_vdwr  = r_123.sum_vdwr
    angle_312 = r_123.angle_312
    for a4 in hierarchy.atoms():
      e4 = a4.element.upper()
      # See Fig.1 in paper "Halogen bond in biological molecules"
      if (e4[0] in ["C", "P", "S"]):
        if (not is_bonded(a2, a4, bond_proxies_simple)): continue
        # theta_2 angle in paper "Halogen bond in biological molecules"
        angle_214 = (a2.angle(a1, a4, deg=True))
        if (120 - angle_eps < angle_214 < 120 + angle_eps):  # 5,line32
          result.append(group_args(
                    atom_1    = a1,
                    atom_2    = a2,
                    atom_3    = a3,
                    atom_4    = a4,
                    d_12      = d,
                    sum_vdwr  = sum_vdwr,
                    d_x_p     = d_x_p,
                    angle_312 = angle_312,
                    angle_214 = angle_214))
    diff_best = 1.e+9
    for r in result:
      diff = abs(120 - r.angle_214)
      if(diff < diff_best):
        diff_best = diff
        result4 = r
  return result4


def find_halogen_bonds(model, eps = 0.15, emp_scale1 = 0.6,emp_scale2 = 0.75, angle_eps=40):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                sites_cart = model.get_sites_cart())
  hierarchy = model.get_hierarchy()
  vdwr      = model.get_vdw_radii()
  halogens  = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ["S", "O", "N","F","CL","BR","I"]
  atom2_list   = []
  final_result = []
  for a1 in hierarchy.atoms():
    e1 = a1.element.upper()
    n1 = a1.name.strip().upper()
    if(e1 in halogens):
      for a2 in hierarchy.atoms():
        e2 = a2.element.upper()
        n2 = a2.name.strip().upper()
        if(not a1.is_in_same_conformer_as(a2)): continue
        if(is_bonded(a1,a2, bond_proxies_simple)): continue
        if(a1.parent().parent().resseq == a2.parent().parent().resseq): continue
        if(e2 in halogen_bond_pairs_atom):
          # O2' in 3v04.pdb file will recognized as O2* ,so replace it
          n1            = n1.replace("'","*")
          n2            = n2.replace("'","*")
          n2            = n2.replace("XT","")# 2yj8.pdb  vdwr can't recognize 'OXT'
          d             = a1.distance(a2)
          if n1 not in vdwr.keys(): continue
          if n2 not in vdwr.keys(): continue
          sum_vdwr      = vdwr[n1] + vdwr[n2]
          sum_vdwr_min1 = sum_vdwr * emp_scale1
          sum_vdwr_min2 = sum_vdwr*emp_scale2 #4 line 31
          if(sum_vdwr_min2-eps < d < sum_vdwr+eps):# found HB pairs-candidates
            result_123 = find_atom_3(d, sum_vdwr, hierarchy, a1, a2, bond_proxies_simple)
            if result_123 is None:continue
            result4 = find_atom4(hierarchy, result_123, bond_proxies_simple, angle_eps)
            if result4 is None:continue
            final_result.append(result4)
            atom2_list.append(a2)
          if atom2_list is None:
           if(sum_vdwr_min1-eps < d < sum_vdwr_min2-eps ):
            result_123 = find_atom_3(d, sum_vdwr, hierarchy, a1, a2, bond_proxies_simple)
            if result_123 is None: continue
            result4    = find_atom4(hierarchy,result_123, bond_proxies_simple, angle_eps)
            if result4 is None: continue
            final_result.append(result4)
  return final_result

#Second step,find salt bridge in one pdb filess

Amino_Acids = ["ARG","HIS","LYS","ASP","GLU","SER","THR",
               "ASN","GLU","CYS","SEC","GLY","PRO","ALA",
               "VAL","ILE","LEU","MET","PHE","TYR","TRP"]

def find_hydrogen_bonds(model, eps = 0.15,emp_scale = 0.75):
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                     sites_cart=model.get_sites_cart())
    hierarchy = model.get_hierarchy()
    vdwr = model.get_vdw_radii()
    final_result = []
    atom2_list = []
    for a1 in hierarchy.atoms():
      e1 = a1.element.upper()
      n1 = a1.name.strip().upper()
      if (e1[0]=="H"):
        for a2 in hierarchy.atoms():
          n2 = a2.name.strip().upper()
          if (not a1.is_in_same_conformer_as(a2)): continue
          if (is_bonded(a1, a2, bond_proxies_simple)): continue
          dict_h_bond_lengh = {"O": 0.98, "N": 1.01, "F": 0.92}
          if n2[0] in dict_h_bond_lengh.keys():
            if n1 not in vdwr.keys(): continue
            if n2 not in vdwr.keys(): continue
            d_12 = a1.distance(a2)
            sum_vdwr = vdwr[n1] + vdwr[n2]
            sum_vdwr_min1 = dict_h_bond_lengh[n2[0]]
            sum_vdwr_min2 = sum_vdwr * emp_scale  # 4 line 31
            if (sum_vdwr_min2 - eps < d_12 < sum_vdwr + eps):  # found HB pairs-candidates
              result_123 = find_atom_3(d_12, sum_vdwr, hierarchy, a1, a2, bond_proxies_simple)
              atom2_list.append(a2)
              if result_123 is None: continue
              final_result.append(result_123)
            if atom2_list is None:
              if (sum_vdwr_min1 + eps < d_12 < sum_vdwr_min2 - eps):
                result_123 = find_atom_3(d_12, sum_vdwr, hierarchy, a1, a2, bond_proxies_simple)
                if result_123 is None: continue
                final_result.append(result_123)
    return final_result

# in the twenty one Amino Acids,[arg,his,lys] are positive Amino Acids
#  with electrically charged side chains,the positive charged atom is N
# there are some H atoms around this N atom
# [asp,glu] are negative,the negative charged atom is O
# the N atom and the O atom will make up the iron bond,the o atom and
# one of the H atom make up the H bond
def find_salt_bridge(model,eps = 0.3):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
        sites_cart=model.get_sites_cart())
  hierarchy      = model.get_hierarchy()
  vdwr           = model.get_vdw_radii()
  result         = []
  atom2_list     = []
  positive_acide = ["ARG" , "HIS" , "LYS"]
  negative_acids = ["ASP" , "GLU"]
  H_acides       = ["LYS" , "ARG" , "HIS" , "SER" , "THR" , "TYP"]
  #first find ions bonds from the first following line to ---
  for a1 in hierarchy.atoms():
     if (a1.parent().resname in positive_acide):
        n1 = filter(str.isalpha,a1.name.upper() )
        if n1[0] == "N" :
          for a2 in hierarchy.atoms():
            if (a2.parent().resname in negative_acids):
              if (is_bonded(a1, a2, bond_proxies_simple)): continue
              if (not a1.is_in_same_conformer_as(a2)): continue
              n2 = filter(str.isalpha, a2.name.upper())
              if n1 not in vdwr.keys(): continue
              if n2 not in vdwr.keys(): continue
              d_12      = a1.distance(a2)
              sum_vdwr  = vdwr[n1] + vdwr[n2]
              if ( d_12 < sum_vdwr ):
                 atom2_list.append((a1.id_str(),a2.id_str()))
                 # ---to here,now find all the ions bonds
                 print atom2_list
                 """
                 angle_132 = (a3.angle(a1, a2, deg=True))
                 if(90< angle_132):
                          result.append(group_args(
                              atom_1   = a1,
                              atom_2   = a2,
                              atom_3   = a3,
                              d_i      = d_i,
                              d_h      = d_h,
                              sum_vdwr = sum_vdwr,
                              angle_132= angle_132))
  return result
"""
def define_pi_system(model,eps = 5):
 pi_amino_acids = ["HIS","PRO","PHE","TYR","TYP"]
 geometry = model.get_restraints_manager()
 planarity_proxies_simple = geometry.geometry.planarity_proxies
 hierarchy = model.get_hierarchy()
 vdwr = model.get_vdw_radii()
 for proxy in planarity_proxies_simple :
   print (proxy.i_seqs)






