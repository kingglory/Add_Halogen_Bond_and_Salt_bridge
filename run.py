from __future__ import division
import iotbx.pdb
import iotbx.cif
import scitbx
import matplotlib as np
from scitbx.array_family import flex
import numpy as np
import mmtbx.hydrogens.build_hydrogens
from libtbx import group_args
from libtbx import easy_run
import math, time
import scitbx.matrix

def is_bonded(atom_1, atom_2, bps_dict):
  i12 = [atom_1.i_seq, atom_2.i_seq]
  i12.sort()
  if(not tuple(i12) in bps_dict): return False
  else: return True

def in_plain(atom1,atom2,atom3,atom4,atom5,pps_dict):
  i12345 = [atom1.i_seq,atom2.i_seq,atom3.i_seq,atom4.i_seq,atom5.iseq]
  i12345.sort()
  if (not tuple(i12345) in pps_dict):return False
  else:return True

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
                       emp_scale2 = 0.75, angle_eps = 30):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                  sites_cart=model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  hierarchy = model.get_hierarchy()
  vdwr = model.get_vdw_radii()
  #main_chain_atoms_plus = ["CA", "N", "O", "C", "CB"]
  halogens = ["CL", "BR", "I", "F"]
  halogen_bond_pairs_atom = [ "O", "S", "N", "F", "CL", "BR", "I"]
  acceptor_pair = ["C","N","P","S"]
  atom1s = []
  atom2s = []
  atom3s = []
  atom4s = []
  results = []
  for a in hierarchy.atoms():
    if a.element.strip().upper() in halogens:
      atom1s.append(a)
    if a.element.strip().upper() in halogen_bond_pairs_atom:
      if a.parent().resname == "HOH":continue
      atom2s.append(a)
    if a.element.strip().upper() in acceptor_pair:
      atom4s.append(a)
    if a.element.strip().upper() == "C":
      atom3s.append(a)
      
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
      if (sum_vdwr_min2 - eps < d < sum_vdwr + eps):
        diff_best = 1.e+9
        result = None
        for a3 in atom3s:
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
                    atom_2=a2)
            if (result is not None):results.append(result)
  return results


def find_hydrogen_bonds(model, min = 1.7, max = 2.2,eps = 0.8):
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                     sites_cart=model.get_sites_cart())
    bps_dict = {}
    [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
    hierarchy = model.get_hierarchy()
    atom1s = []
    atom2s = []
    atom3s = []
    atom4s = []
    results = []
    Accepter_H_pair = ["O","N"]
    for a in hierarchy.atoms():
      e = a.element.strip().upper()
      if a.element_is_hydrogen():
        atom1s.append(a)
      if e == "O":
        if a.parent().resname == "HOH":continue
        atom2s.append(a)
      if e in Accepter_H_pair:
        atom4s.append(a)
      if e == "C":
        atom3s.append(a)

    for a2 in atom2s:
      result = None
      diff_best = 1.e+9
      for a3 in atom4s:
        for a1 in atom1s:
          if (not a1.is_in_same_conformer_as(a2)): continue
          if (not is_bonded(a1, a3, bps_dict)): continue
          if (is_bonded(a1, a2, bps_dict)): continue
          if (a1.parent().parent().resseq ==
              a2.parent().parent().resseq): continue
          d_12 = a1.distance(a2)
          if (min-eps < d_12 < max+eps):
            angle_312 = (a1.angle(a2, a3, deg=True))
            if (100 < angle_312):
              diff = abs( 2 - d_12 )
              if (diff < diff_best):
                  diff_best = diff
                  result = group_args(
                    atom_1=a1,
                    atom_2=a2)
        if (result in results):continue
        if (result is not None): results.append(result)
      if (result is None):
        for a3 in atom3s:
          for a1 in atom1s:
            if (not a1.is_in_same_conformer_as(a2)): continue
            if (not is_bonded(a1, a3, bps_dict)): continue
            if (is_bonded(a1, a2, bps_dict)): continue
            if (a1.parent().parent().resseq ==
                  a2.parent().parent().resseq): continue
            d_12 = a1.distance(a2)
            if (min - eps < d_12 < max + eps):
              angle_312 = (a1.angle(a2, a3, deg=True))
              if (100 < angle_312):
                diff = abs(180 - angle_312)
                if (diff < diff_best):
                  diff_best = diff
                  result = group_args(
                    atom_1=a1,
                    atom_2=a2)
        if (result in results): continue
        if (result is not None): results.append(result)
    # one O atom can make two hydrogen bonds,
    # but one hydrogen atom can just one

    for i,ri in enumerate(results):
      for j,rj in enumerate(results):
        if (j <= i): continue
        ai1 = ri.atom_1
        ai2 = ri.atom_2
        aj1 = rj.atom_1
        aj2 = rj.atom_2
        if ri.atom_1 == rj.atom_1:
          di = ai1.distance(ai2)
          dj = aj1.distance(aj2)
          if di < dj:
            results.remove(rj)
          else:
            results.remove(ri)

    return results


def find_salt_bridge(model, pdb_file_name, min = 1.7, max = 2.2,
                     eps1 = 0.15, eps2 = 0.8, shutoff = 4 ):
  salt_bridge_res = ["resname ARG", "resname HIS",
                     "resname LYS", "resname ASP", "resname GLU"]
  ss = " or ".join(salt_bridge_res)
  m_sel = model.selection(ss)
  new_model = model.select(m_sel)
  hierarchy = new_model.get_hierarchy()
  hierarchy.atoms().reset_i_seq()
  crystal_symmetry = new_model.crystal_symmetry()
  hierarchy.write_pdb_file(file_name=pdb_file_name[0:6] + "new.pdb",
                           crystal_symmetry=crystal_symmetry)
  geometry = new_model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
    sites_cart=new_model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  positive_residues = ["ARG", "HIS", "LYS"]
  negative_residues = ["ASP", "GLU", "HIS"]
  #main_chain_atoms_plus = ["CA","N","O","C","CB"]
  results = []
  positive_atoms = []
  atom1s = []
  atom3s = []
  """ select out the pasitive atoms and hydrogen atoms 
  and negative atoms to lists
  """
  for a in hierarchy.atoms():
    e = a.element.strip().upper()
    if a.element_is_hydrogen():
      #if a.parent().resname == "HOH":continue
      atom1s.append(a)
    if e == "N":
      if not (a.parent().resname in positive_residues):continue
      positive_atoms.append(a)
    if e == "O":
      #if a.parent().resname == "HOH":continue
      if not (a.parent().resname in negative_residues):continue
      atom3s.append(a)

   
      
  """ select out N H pairs in pasitive sites
   if there are more than three hydrohen atoms
   make covalent bond with N in side chains
   then the N atom could be positive atom
  """
  i = 0
  a1_a2_pairs = []
  H_N_pairs = []
  for a2 in positive_atoms:
    for a1 in atom1s:
      if (not is_bonded(a1, a2, bps_dict)): continue
      i = i+1
      a1_a2_pairs.append((a1,a2))
    if i >=2:
      if a1_a2_pairs is None:continue
      H_N_pairs.extend(a1_a2_pairs)
      
  """ select out the negative site that opposite the positive site,
    if the O atom in negative is close enough with N (positive atom),
    and they didn't make covalent bond ,they will have chance to make
    electrostatic interaction.and the negative also will have chance 
    to make hydrogen bonds with th hydrogen atoms nearby the positive
    N atom.then the eletrostatic interaction will make salt bridge
    with hydrogen bonds.
  """
  for a3 in atom3s:
    result = None
    diff_best = 1.e+9
    for r in H_N_pairs:
      a1 = r[0]
      a2 = r[1]
      if (is_bonded(a3, a2, bps_dict)): continue
      if (not a3.is_in_same_conformer_as(a2)): continue
      if (a3.parent().parent().resseq ==
              a2.parent().parent().resseq): continue
      d_32 = a3.distance(a2)
      if(d_32 > shutoff - eps1): continue
      d_13 = a1.distance(a3)
      if (min-eps1 < d_13 < max+eps2):
        angle_312 = (a1.angle(a2, a3, deg=True))
        if (100 < angle_312):
          diff = abs(180 - angle_312)
          if (diff < diff_best):
            diff_best = diff
            '''
            result = group_args(atom_1 = a1.id_str(),
                                atom_2 = a2.id_str(),
                                atom_3 = a3.id_str())
            '''
            result = group_args(atom_1 = a1,
                                atom_2 = a2,
                                atom_3 = a3)

    if (result is not None): results.append(result)

  return results
        

def define_pi_system(model, dist_cutoff_1=6.0,dist_cutoff_2=3.0,
                     dist_h_cutoff=2.5,dist_v_cutoff=1.2,T_angle = 90,
                     P_angle_1 = 0,eps_angle = 30,P_angle_2 =180):

  results = []
  planes = []
  '''
  #Ring_containing_amino_acid(Ring_CAA)
  #Ring_CAA = ["resname HIS","resname PHE",
  #            "resname TYR","resname TRP"]
  #  asc = hierarchy.atom_selection_cache()
  #ss = " or ".join(Ring_CAA)
  #hierarchys = hierarchy.select(asc.selection(ss))
  #m_sel = model.selection(ss + "and not name CB")
  #m_sel = model.selection(ss)
  #new_model = model.select(m_sel)
 
  #hierarchy.atoms().reset_i_seq()
  crystal_symmetry = model.crystal_symmetry()
  #hierarchy.write_pdb_file(file_name=pdb_file_name[0:4]+"new.pdb",
  #                         crystal_symmetry=crystal_symmetry)
  '''
  hierarchy = model.get_hierarchy()
  geometry = model.get_restraints_manager()
  # first limition: the pi is a plane
  atoms = hierarchy.atoms()
  for proxy in geometry.geometry.planarity_proxies:
    planes.append(atoms.select(proxy.i_seqs))
  for i, pi in enumerate(planes):
    for j, pj in enumerate(planes):
      if(j<=i): continue
      xyzi = pi.extract_xyz()
      xyzj = pj.extract_xyz()
      ni = list(pi.extract_name())
      nj = list(pj.extract_name())
      if (' CG ' not in ni): continue
      if (' CG ' not in nj): continue
      if ('CA' in ni): continue
      if ('CA' in nj): continue
      if (len(ni) < 5):continue
      if (len(nj) < 5): continue
      # The first atom name must be CG in a protein ring
      """ second limition: the mean distance between two plane 
      is short than 5
      """
      cmi = xyzi.mean()
      cmj = xyzj.mean()
      dist = math.sqrt(
        (cmi[0]-cmj[0])**2+
        (cmi[1]-cmj[1])**2+
        (cmi[2]-cmj[2])**2)
      if(dist>dist_cutoff_1): continue

      ai,bi,ci,di = scitbx.matrix.plane_equation(
        point_1=scitbx.matrix.col(xyzi[0]), 
        point_2=scitbx.matrix.col(xyzi[1]), 
        point_3=scitbx.matrix.col(xyzi[2]))

      aj,bj,cj,dj = scitbx.matrix.plane_equation(
        point_1=scitbx.matrix.col(xyzj[0]), 
        point_2=scitbx.matrix.col(xyzj[1]), 
        point_3=scitbx.matrix.col(xyzj[2]))
      # calculation the distance from one center point of one plane to another plane
      x, y, z = cmj
      den = math.sqrt(ai ** 2 + bi ** 2 + ci ** 2)
      dist_point_plane = abs(ai * x + bi * y + ci * z + di) / den
      # calculation the dist_h_d between two center of plane in horizontal direction
      # for T type dist_point_plane  also represent one of dist_h_d
      dist_h_d = math.sqrt(dist**2 - dist_point_plane**2)

      """
      third limition : Dihedral angle
      for T(perpendicular)is 180
      for P(In parallel) is 90
      30 degree angular offset
      https://math.tutorvista.com/geometry/angle-between-two-planes.html
      """
      cos_a = abs(ai*aj+bi*bj+ci*cj)/(ai**2+bi**2+ci**2)**0.5/\
        (aj**2+bj**2+cj**2)**0.5
      angle = math.acos(cos_a)*180/math.pi
      # for T type,the angle should be near to 90,but 30 deviation is ok
      if (P_angle_1 + eps_angle) < angle < (T_angle - eps_angle):continue
      if (T_angle - eps_angle) < angle < (T_angle + eps_angle):
        if ((dist_point_plane > dist_v_cutoff )and(dist_h_d > dist_v_cutoff)):continue
      # for P type,the angle shuld be near to 180,but 30 deviaton is ok
      if (P_angle_2 - eps_angle) > angle > (T_angle + eps_angle):continue
      if (eps_angle > angle > P_angle_1)or( P_angle_2 > angle > P_angle_2 - eps_angle):
        if (dist_h_d > dist_h_cutoff):continue
        if (dist_point_plane < dist_cutoff_2): continue
      print angle, dist,dist_point_plane,dist_h_d
      result = group_args( pi  = pi,
                           pj  = pj
                           )
      if result in results:continue
      if (result is not None): results.append(result)
  return results





