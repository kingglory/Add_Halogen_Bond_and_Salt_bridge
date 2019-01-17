from __future__ import division
import iotbx.pdb
import iotbx.cif
import scitbx
import numpy as np
import mmtbx.hydrogens.build_hydrogens
from libtbx import group_args
from libtbx import easy_run
import math
import scitbx.matrix

def is_bonded(atom_1, atom_2, bps_dict):
  i12 = [atom_1.i_seq, atom_2.i_seq]
  i12.sort()
  if(not tuple(i12) in bps_dict): return False
  else: return True

def in_plain(atom1,atom2,atom3,atom4,pps_dict):
  i1234 = [atom1.i_seq,atom2.i_seq,atom3.i_seq,atom4.i_seq]
  i1234.sort()
  if (not tuple(i1234) in pps_dict):return False
  else:return True

# first step write codes to find halogen bond in one pdb file
# define a function to try finding the halogen bond pairs

"""
1,when one halogen atoms make halogen bond when other atom,
when the distance and angle limations all fited,
the bond distance more shorter ,more possible;
2, when chooseing the third atoms to make up the angle1,
the angle1 is more near 180,more possible;
3,when chooseing the fourth atoms to make up the angle2,
the angle1 is more near 120,more possible;
4,see Figure 4 in paper "Halogen bonding(X-bonding): 
 A biological perspective"
5,theta_2 angle :Geometry of X-bonds in paper 
"Halogen bond in biological molecules"
"""

def find_halogen_bonds(model, eps = 0.15, emp_scale1 = 0.6,
                       emp_scale2 = 0.75, angle_eps = 40):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                sites_cart = model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ["S", "O", "N","F","CL","BR","I"]
  results = []
  for a1 in hierarchy.atoms():
    e1 = a1.element.strip().upper()
    n1 = a1.name.strip().upper()
    if(e1 in halogens):
      for a2 in hierarchy.atoms():
        e2 = a2.element.strip().upper()
        n2 = a2.name.strip().upper()
        if(not a1.is_in_same_conformer_as(a2)): continue
        if(is_bonded(a1, a2, bps_dict)): continue
        if(a1.parent().parent().resseq ==
             a2.parent().parent().resseq): continue
        if(e2 in halogen_bond_pairs_atom):
          # O2' in 3v04.pdb file will recognized as O2* ,
          #  so replace it
          n1 = n1.replace("'","*")
          n2 = n2.replace("'","*")
          n2 = n2.replace("XT","")
          # 2yj8.pdb  vdwr can't recognize 'OXT'
          d  = a1.distance(a2)
          if(n1 not in vdwr.keys()): continue
          if(n2 not in vdwr.keys()): continue
          sum_vdwr = vdwr[n1] + vdwr[n2]
          sum_vdwr_min2 = sum_vdwr * emp_scale2
          d_x_p = d / sum_vdwr
          if(sum_vdwr_min2-eps < d < sum_vdwr+eps):
            # found HB pairs-candidates
            diff_best = 1.e+9
            result = None
            for a3 in hierarchy.atoms():
              if(not is_bonded(a1, a3, bps_dict)): continue
              # theta_1 angle in paper "Halogen bond in biological molecules"
              angle_312 = (a1.angle(a2, a3, deg=True))
              # Fig.1 in "Halogen bond in biological molecules"
              #for most cases,the angle 1 is between from 130 to 180,
              # we can search in this range
              #if find something ,done!if didn't find something .
              # searching the following range
              #from 90 to 130 degrees
              if(90 < angle_312):
                for a4 in hierarchy.atoms():
                  e4 = a4.element.strip().upper()
                  # Fig.1 in "Halogen bond in biological molecules"
                  if(e4 in ["C", "P", "S"]):
                    if(not is_bonded(a2, a4, bps_dict)): continue
                    # theta_2 angle in "Halogen bond in biological molecules"
                    angle_214 = (a2.angle(a1, a4, deg=True))
                    if(120 - angle_eps < angle_214 < 120 + angle_eps):
                      diff = abs(120 - angle_214)
                      if(diff < diff_best):
                        diff_best = diff
                        result = group_args(
                          atom_1 = a1,
                          atom_2 = a2,
                          atom_3 = a3,
                          atom_4 = a4,
                          d_12 = d,
                          sum_vdwr = sum_vdwr,
                          d_x_p = d_x_p,
                          angle_312 = angle_312,
                          angle_214 = angle_214)
                        if(result is not None):results.append(result)
  return results


def f_halogen_bonds(model, eps = 0.15, emp_scale1 = 0.6,
                       emp_scale2 = 0.75, angle_eps = 40):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
    sites_cart=model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = ["S", "O", "N", "F", "CL", "BR", "I"]
  atom1s = []
  atom2s = []
  atom4s = []
  results = []
  for a in hierarchy.atoms():
    if a.element.strip().upper() in halogens:
      atom1s.append(a)
    if a.element.strip().upper() in halogen_bond_pairs_atom:
      atom2s.append(a)
    if a.element.strip().upper() in ["C","P","S"]:
      atom4s.append(a)
  for a1 in atom1s:
    for a2 in atom2s:
      if (not a1.is_in_same_conformer_as(a2)): continue
      if (is_bonded(a1, a2, bps_dict)): continue
      if (a1.parent().parent().resseq ==
            a2.parent().parent().resseq): continue
      d = a1.distance(a2)
      n1 = a1.name.strip().upper()
      n2 = a2.name.strip().upper()
      if (n1 not in vdwr.keys()): continue
      if (n2 not in vdwr.keys()): continue
      sum_vdwr = vdwr[n1] + vdwr[n2]
      sum_vdwr_min2 = sum_vdwr * emp_scale2
      d_x_p = d / sum_vdwr
      if (sum_vdwr_min2 - eps < d < sum_vdwr + eps):
        diff_best = 1.e+9
        result = None
        for a3 in hierarchy.atoms():
          if (not is_bonded(a1, a3, bps_dict)): continue
          angle_312 = (a1.angle(a2, a3, deg=True))
          if (90 < angle_312):
            for a4 in atom4s:
              if (not is_bonded(a2, a4, bps_dict)): continue
              # theta_2 angle in "Halogen bond in biological molecules"
              angle_214 = (a2.angle(a1, a4, deg=True))
              if (120 - angle_eps < angle_214 < 120 + angle_eps):
                diff = abs(120 - angle_214)
                if (diff < diff_best):
                  diff_best = diff
                  result = group_args(
                    atom_1=a1,
                    atom_2=a2,
                    atom_3=a3,
                    d_12=d,
                    angle_312=angle_312,
                    angle_214=angle_214)
            if (result is not None):results.append(result)
  return results
#Second step,find salt bridge in one pdb filess

Amino_Acids = ["ARG","HIS","LYS","ASP","GLU","SER","THR",
               "ASN","GLU","CYS","SEC","GLY","PRO","ALA",
               "VAL","ILE","LEU","MET","PHE","TYR","TRP"]

def find_hydrogen_bonds(model, eps = 0.15,emp_scale = 0.75):
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                     sites_cart=model.get_sites_cart())
    bps_dict = {}
    [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
    hierarchy = model.get_hierarchy()
    vdwr = model.get_vdw_radii()
    results = []
    for a1 in hierarchy.atoms():
      n1 = a1.name.strip().upper()
      if a1.element_is_hydrogen():
        for a2 in hierarchy.atoms():
          e2 = a2.element.strip().upper()
          n2 = a2.name.strip().upper()
          if (not a1.is_in_same_conformer_as(a2)): continue
          if (is_bonded(a1, a2, bps_dict)): continue
          dict_h_bond_lengh = {"O": 0.98, "N": 1.01, "F": 0.92}
          if e2 in dict_h_bond_lengh.keys():
            if n1 not in vdwr.keys(): continue
            if n2 not in vdwr.keys(): continue
            d_12 = a1.distance(a2)
            sum_vdwr = vdwr[n1] + vdwr[n2]
            d_x_p = d_12 / sum_vdwr
            sum_vdwr_min1 = dict_h_bond_lengh[n2[0]]
            sum_vdwr_min2 = sum_vdwr * emp_scale  # 4 line 31
            if (sum_vdwr_min2 - eps < d_12 < sum_vdwr + eps):
              # found HB pairs-candidates
              for a3 in hierarchy.atoms():
                if (not is_bonded(a1, a3, bps_dict)): continue
                angle_312 = (a1.angle(a2, a3, deg=True))
                if (90 < angle_312):
                  result = group_args(
                    atom_1 = a1,
                    atom_2 = a2,
                    atom_3 = a3,
                    d_12 = d_12,
                    sum_vdwr = sum_vdwr,
                    d_x_p = d_x_p,
                    angle_312 = angle_312)
                  if (result is not None): results.append(result)

    return results



def f_hydrogen_bonds(model, eps = 0.15,emp_scale = 0.75):
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                     sites_cart=model.get_sites_cart())
    bps_dict = {}
    [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
    hierarchy = model.get_hierarchy()
    vdwr = model.get_vdw_radii()
    atom1s = []
    atom2s = []
    results = []
    dict_h_bond_lengh = {"O": 0.98, "N": 1.01, "F": 0.92}
    for a in hierarchy.atoms():
      e = a.element.strip().upper()
      if a.element_is_hydrogen():
        atom1s.append(a)
      if e in dict_h_bond_lengh.keys():
        atom2s.append(a)
    for a1 in atom1s:
      for a2 in atom2s:
        if (not a1.is_in_same_conformer_as(a2)): continue
        if (is_bonded(a1, a2, bps_dict)): continue
        n1 = a1.name.strip().upper()
        n2 = a2.name.strip().upper()
        if n1 not in vdwr.keys(): continue
        if n2 not in vdwr.keys(): continue
        d_12 = a1.distance(a2)
        sum_vdwr = vdwr[n1] + vdwr[n2]
        sum_vdwr_min2 = sum_vdwr * emp_scale
        if (sum_vdwr_min2 - eps < d_12 < sum_vdwr + eps):
          for a3 in hierarchy.atoms():
            if (not is_bonded(a1, a3, bps_dict)): continue
            angle_312 = (a1.angle(a2, a3, deg=True))
            if (90 < angle_312):
              result = group_args(
                atom_1=a1,
                atom_2=a2,
                atom_3=a3,
                d_12=d_12,
                angle_312=angle_312)
              if (result is not None): results.append(result)
    return results



    # in the twenty one Amino Acids,[arg,his,lys] are positive Amino Acids
#  with electrically charged side chains,the positive charged atom is N
# there are some H atoms around this N atom
# [asp,glu] are negative,the negative charged atom is O
# the N atom and the O atom will make up the iron bond,the o atom and
# one of the H atom make up the H bond
def find_salt_bridge(model, eps = 0.15,emp_scale = 0.75):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                  sites_cart=model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  results = []
  ions_bonds_paris_list = []
  positive_acide = ["ARG" , "HIS" , "LYS"]
  negative_acids = ["ASP" , "GLU"]
  H_acides = ["LYS" , "ARG" , "HIS" , "SER" , "THR" , "TYP"]
  #first find ions bonds from the first following line to ---
  for a1 in hierarchy.atoms():
    #if (a1.parent().resname in positive_acide):
    e1 = a1.element.strip().upper()
    n1 = filter(str.isalpha, a1.name.upper())
    if e1 == "N" :
      for a2 in hierarchy.atoms():
        #if (a2.parent().resname in negative_acids):
        if (is_bonded(a1, a2, bps_dict)): continue
        if (not a1.is_in_same_conformer_as(a2)): continue
        n2 = filter(str.isalpha, a2.name.upper())
        if n1 not in vdwr.keys(): continue
        if n2 not in vdwr.keys(): continue
        d_12 = a1.distance(a2)
        sum_vdwr = vdwr[n1] + vdwr[n2]
        if ( d_12 < sum_vdwr ):
          ions_bonds_paris_list.append((a1.id_str(), a2.id_str()))
          print ions_bonds_paris_list
          # ---to here,now find all the ions bonds;
          #the following is to find the hydrogen bonds beside the ions bonds
          for a3 in hierarchy.atoms():
            n3 = a1.name.strip().upper()
            if a3.element_is_hydrogen():
              if a1.distance(a3) > 8:continue
              for a4 in hierarchy.atoms():
                e4 = a4.element.strip().upper()
                n4 = a4.name.strip().upper()
                if (not a3.is_in_same_conformer_as(a4)): continue
                if (is_bonded(a3, a4, bps_dict)): continue
                dict_h_bond_lengh = {"O": 0.98, "N": 1.01, "F": 0.92}
                if e4 in dict_h_bond_lengh.keys():
                  #if the mins distance between ions atoms and hydrogen bonds paris atoms
                  #  is longer than the longest distance between saltbridge :4,
                  # it won't make up saltbridge
                  d13 = a1.distance(a3)
                  d14 = a1.distance(a4)
                  d23 = a2.distance(a3)
                  d24 = a2.distance(a4)
                  if min(d13,d14,d23,d24)< 8:
                    if n3 not in vdwr.keys(): continue
                    if n4 not in vdwr.keys(): continue
                    d_34 = a1.distance(a2)
                    sum_vdwr = vdwr[n3] + vdwr[n4]
                    d_x_p = d_12 / sum_vdwr
                    sum_vdwr_min2 = sum_vdwr * emp_scale  # 4 line 31
                    if (sum_vdwr_min2 - eps < d_12 < sum_vdwr + eps):
                      # found HB pairs-candidates
                      for a5 in hierarchy.atoms():
                        if (not is_bonded(a1, a3, bps_dict)): continue
                        angle_534 = (a3.angle(a5, a4, deg=True))
                        if (90 < angle_534):
                          result = group_args(
                            atom_1 = a1,
                            atom_2 = a2,
                            atom_3 = a3,
                            atom_4 = a4,
                            atom_5 = a5,
                            d_12 = d_12,
                            d_34 = d_34,
                            sum_vdwr = sum_vdwr,
                            d_x_p = d_x_p,
                            angle_534 = angle_534)
                          if (result is not None): results.append(result)
  return results
def f_ions_bonds(model):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
    sites_cart=model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  ions_bonds_paris_list = []
  for a1 in hierarchy.atoms():
    e1 = a1.element.strip().upper()
    n1 = a1.name.strip().upper()
    if e1 == "N":
      for a2 in hierarchy.atoms():
        if a2.element_is_hydrogen():continue
        if (is_bonded(a1, a2, bond_proxies_simple)): continue
        if (not a1.is_in_same_conformer_as(a2)): continue
        n2 = a2.name.strip().upper()
        if n1 not in vdwr.keys(): continue
        if n2 not in vdwr.keys(): continue
        d_12 = a1.distance(a2)
        sum_vdwr = vdwr[n1] + vdwr[n2]
        if (d_12 < sum_vdwr):
          if a1 is None:continue
          if a2 is None:continue
          ions_bonds_paris_list.append((a1, a2))
  return ions_bonds_paris_list
def f_salt_bridge(model,dist_cutoff=4):
  ions_bonds_paris_list = f_ions_bonds(model)
  results = f_hydrogen_bonds(model=model)
  for r in results:
    a3 = list(r.atom_1)
    a4 = list(r.atom_2)
    cm1 = np.array(list[a3, a4]).mean()
    for ions in ions_bonds_paris_list:
      a1 = ions[0]
      a2 = ions[1]
      cm2 = (a1, a2).mean()
      dist = math.sqrt(
        (cm1[0] - cm2[0]) ** 2 +
        (cm1[1] - cm2[1]) ** 2 +
        (cm1[2] - cm2[2]) ** 2)
      if (dist > dist_cutoff): continue
      result = group_args(
        atom_1=a1,
        atom_2=a2,
        atom_3=a3,
        atom_4=a4)
      if (result is not None): results.append(result)
  return results

def define_pi_system(model, dist_cutoff=5):
  geometry = model.get_restraints_manager()
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  results = []
  planes = []
  for proxy in geometry.geometry.planarity_proxies:
    planes.append(atoms.select(proxy.i_seqs))
  for i, pi in enumerate(planes):
    for j, pj in enumerate(planes):
      if(j<=i): continue
      xyzi = pi.extract_xyz()
      xyzj = pj.extract_xyz()
      ni = list(pi.extract_name())
      nj = list(pj.extract_name())
      if(' CA ' in ni and len(ni)==4): continue
      if(' CA ' in nj and len(nj)==4): continue
      cmi = xyzi.mean()
      cmj = xyzj.mean()
      dist = math.sqrt(
        (cmi[0]-cmj[0])**2+
        (cmi[1]-cmj[1])**2+
        (cmi[2]-cmj[2])**2)
      if(dist>dist_cutoff): continue
      ai,bi,ci,di = scitbx.matrix.plane_equation(
        point_1=scitbx.matrix.col(xyzi[0]), 
        point_2=scitbx.matrix.col(xyzi[1]), 
        point_3=scitbx.matrix.col(xyzi[2]))
      aj,bj,cj,dj = scitbx.matrix.plane_equation(
        point_1=scitbx.matrix.col(xyzj[0]), 
        point_2=scitbx.matrix.col(xyzj[1]), 
        point_3=scitbx.matrix.col(xyzj[2]))
      # https://math.tutorvista.com/geometry/angle-between-two-planes.html
      cos_a = abs(ai*aj+bi*bj+ci*cj)/(ai**2+bi**2+ci**2)**0.5/\
        (aj**2+bj**2+cj**2)**0.5
      angle = math.acos(cos_a)*180/math.pi
      if angle < 0:continue
      if angle > 90:continue
      result = group_args( angle = angle,
                           dist = dist,
                           p_i = list(pi.extract_i_seq()),
                           p_j = list(pj.extract_i_seq()),
                           pi_atoms_name = list(pi.extract_name()),
                           pj_atoms_name = list(pj.extract_name())
                           )
      if result in results:continue
      if (result is not None): results.append(result)
  return results





