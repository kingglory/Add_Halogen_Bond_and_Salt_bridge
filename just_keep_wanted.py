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

def define_pi_system(model, dist_cutoff_1=6.0, dist_cutoff_2=3.0,
                     dist_h_cutoff=2.5, dist_v_cutoff=1.2, T_angle=90,
                     P_angle_1=0, eps_angle=25, P_angle_2=180):
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
      if (j <= i): continue
      xyzi = pi.extract_xyz()
      xyzj = pj.extract_xyz()
      ni = list(pi.extract_name())
      nj = list(pj.extract_name())
      if (' CG ' not in ni): continue
      if (' CG ' not in nj): continue
      if ('CA' in ni): continue
      if ('CA' in nj): continue
      if (len(ni) < 5): continue
      if (len(nj) < 5): continue
      # The first atom name must be CG in a protein ring
      """ second limition: the mean distance between two plane 
      is short than 5
      """
      cmi = xyzi.mean()
      cmj = xyzj.mean()
      dist = math.sqrt(
        (cmi[0] - cmj[0]) ** 2 +
        (cmi[1] - cmj[1]) ** 2 +
        (cmi[2] - cmj[2]) ** 2)
      if (dist > dist_cutoff_1): continue

      ai, bi, ci, di = scitbx.matrix.plane_equation(
        point_1=scitbx.matrix.col(xyzi[0]),
        point_2=scitbx.matrix.col(xyzi[1]),
        point_3=scitbx.matrix.col(xyzi[2]))

      aj, bj, cj, dj = scitbx.matrix.plane_equation(
        point_1=scitbx.matrix.col(xyzj[0]),
        point_2=scitbx.matrix.col(xyzj[1]),
        point_3=scitbx.matrix.col(xyzj[2]))
      # calculation the distance from one center point of one plane to another plane
      x, y, z = cmj
      den = math.sqrt(ai ** 2 + bi ** 2 + ci ** 2)
      dist_point_plane = abs(ai * x + bi * y + ci * z + di) / den
      # calculation the dist_h_d between two center of plane in horizontal direction
      # for T type dist_point_plane  also represent one of dist_h_d
      dist_h_d = math.sqrt(dist ** 2 - dist_point_plane ** 2)

      """
      third limition : Dihedral angle
      for T(perpendicular)is 180
      for P(In parallel) is 90
      30 degree angular offset
      https://math.tutorvista.com/geometry/angle-between-two-planes.html
      """
      cos_a = abs(ai * aj + bi * bj + ci * cj) / (ai ** 2 + bi ** 2 + ci ** 2) ** 0.5 / \
              (aj ** 2 + bj ** 2 + cj ** 2) ** 0.5
      angle = math.acos(cos_a) * 180 / math.pi
      # for T type,the angle should be near to 90,but 30 deviation is ok
      if (P_angle_1 + eps_angle) < angle < (T_angle - eps_angle): continue
      if (T_angle - eps_angle) < angle < (T_angle + eps_angle):
        if ((dist_point_plane > dist_v_cutoff) and (dist_h_d > dist_v_cutoff)): continue
      # for P type,the angle shuld be near to 180,but 30 deviaton is ok
      if (P_angle_2 - eps_angle) > angle > (T_angle + eps_angle): continue
      if (eps_angle > angle > P_angle_1) or (P_angle_2 > angle > P_angle_2 - eps_angle):
        if (dist_h_d > dist_h_cutoff): continue
        if (dist_point_plane < dist_cutoff_2): continue
      # print angle, dist,dist_point_plane,dist_h_d
      result = group_args(pi=pi,
                          pj=pj
                          )
      if result in results: continue
      if (result is not None): results.append(result)
  return results




