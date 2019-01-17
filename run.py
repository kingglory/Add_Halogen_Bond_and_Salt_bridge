from __future__ import division
import numpy as np
import mmtbx.hydrogens.build_hydrogens
from libtbx import group_args
from libtbx import easy_run

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
    atom3s = []
    results = []
    dict_h_bond_lengh = {"O": 0.98, "N": 1.01, "F": 0.92}
    for a in hierarchy.atoms():
      e = a.element.strip().upper()
      if a.element_is_hydrogen():
        atom1s.append(a)
      if e in dict_h_bond_lengh.keys():
        atom2s.append(a)
      if not a.element_is_hydrogen():
        atom3s.append(a)
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
          for a3 in atom3s:
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
  print ions_bonds_paris_list
  return ions_bonds_paris_list
def f_salt_bridge(model):
  # return ion bonds atoms paris
  ions_bonds_paris_list = f_ions_bonds(model)
  #return hydrogen bonding atoms paris
  results = f_hydrogen_bonds(model=model)
  #calculate the distance between hydrogen bonds and ions bonds
  #is the distance is shorter enough ,then salt bridge appear
  for r in results:
    a3 = r.atom_1
    a4 = r.atom_2
    for ions in ions_bonds_paris_list:
      a1 = ions[0]
      a2 = ions[1]
      d13 = a1.distance(a3)
      d14 = a1.distance(a4)
      d23 = a2.distance(a3)
      d24 = a2.distance(a4)
      d_ions_H = [d13, d14, d23, d24]
      d_average = sum(d_ions_H)/len(d_ions_H)
      if d_average <=4:
        result = group_args(
          atom_1=a1,
          atom_2=a2,
          atom_3=a3,
          atom_4=a4)
        if (result is not None): results.append(result)
  return results


def distance_atomToArea(atom1_xyz, atom2_xyz, atom3_xyz, atom4_xyz):
  """
  Description: The distance from atom 4 to atom 1, atom 2, atom 3
  Param atom1_xyz: row slices of data boxes, three-dimensional
  Param atom2_xyz:
  Param atom3_xyz:
  Param atom_xyz:
  Return: distance from atom to plane
  """
  #Let's assume that the X of the normal vector (x, y, z) is equal to 1.
  x = 1# X coordinates of normal vectors
  p1 = np.array(atom1_xyz)  #Conversion to Matrix
  p2 = np.array(atom2_xyz)
  p3 = np.array(atom3_xyz)
  p4 = np.array(atom4_xyz)
  v12 = p1 - p2  #Vectors from P1 to P2
  v13 = p1 - p3  #Vectors from P1 to P3
  X = np.vstack((v12, v13))
  #Merge vectors V12 and V13 into 2*3 matrices
  a = -X[:, 0] #The coefficient of a
  yz = np.matrix(X[:, 1:]).I.dot(a)
  #Multiply the inverse matrix of X by a
  y = yz[0, 0]  #Y coordinates of normal vectors
  z = yz[0, 1]  # Z coordinates of normal vectors
  n = np.array([x, y, z])
  #Normal vectors of the plane to which the first three atoms belong
  v41 = p4 - p1  #Vectors from P4 to P1
  return abs(v41.dot(n) / n.dot(n))




def calculate_n(plane1):
  p1_atoms = []
  for atom1 in plane1:
    p1_atoms.append(atom1)
    if len(p1_atoms) == 3:
      a1 = p1_atoms[0]
      a2 = p1_atoms[1]
      a3 = p1_atoms[2]
      x = 1  # X coordinates of normal vectors
      p1 = np.array(a1.xyz)  # Conversion to Matrix
      p2 = np.array(a2.xyz)
      p3 = np.array(a3.xyz)
      v12 = p1 - p2  # Vectors from P1 to P2
      v13 = p1 - p3  # Vectors from P1 to P3
      X = np.vstack((v12, v13))
      # Merge vectors V12 and V13 into 2*3 matrices
      a = -X[:, 0]  # The coefficient of a
      yz = np.matrix(X[:, 1:]).I.dot(a)
      # Multiply the inverse matrix of X by a
      y = yz[0, 0]  # Y coordinates of normal vectors
      z = yz[0, 1]  # Z coordinates of normal vectors
      n = np.array([x, y, z])
      # Normal vectors of the plane to which the first three atoms belong
    #if len(p1_atoms) > 3: continue
  return n

def define_pi_system(model):
  geometry = model.get_restraints_manager()
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  results = []
  planes = []
  d12_list = []
  for proxy in geometry.geometry.planarity_proxies:
    planes.append(atoms.select(proxy.i_seqs))
  for plane in planes:
   plane1 = plane
   n = calculate_n(plane1)#Calculate the normal vector of plane 1
   for atom1 in plane1:
     for plane in planes:
       if plane == plane1: continue
       plane2 = plane
       for atom2 in plane2:
         # Converting atomic coordinates into vectors
         atom1_v = np.array(atom1.xyz)
         atom2_v = np.array(atom2.xyz)
         v_12 = atom1_v - atom2_v
         # Vertors from atom1 to atom2
         d_12 = abs(v_12.dot(n) / n.dot(n))
         d12_list.append(d_12)
         # calculate the average of distance from plane 1 atoms to plane 2 atoms
       d_average = sum(d12_list) / len(d12_list)
       if 2 < d_average < 5:
         result = group_args(plane1 = plane1, plane2 = plane2)
         if result in results:continue
         if (result is not None): results.append(result)
  return results





