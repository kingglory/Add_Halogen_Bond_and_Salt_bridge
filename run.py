from __future__ import division
import scitbx
from libtbx import group_args
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
class get_halogen_bonds(object):
  def __init__(self,model):
    self.model = model
    self.results = self.get_halogen_bonds_pairs()

  def get_halogen_bonds_pairs(self,eps = 0.15, emp_scale1 = 0.6,
                         emp_scale2 = 0.75, angle_eps = 30):
    geometry = self.model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                sites_cart=self.model.get_sites_cart())
    bps_dict = {}
    [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
    hierarchy = self.model.get_hierarchy()
    vdwr = self.model.get_vdw_radii()
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
                      atom_2=a2,
                      atom_3=a3,
                      d=d,
                      angle_312=angle_312)
              if (result is not None):results.append(result)
    return results

  def write_restrains_file(self, pdb_file_name,
                           for_phenix_refine=True,
                           use_defaul_parameters=True):
    str_1 = '''bond{
      atom_selection_1 = %s
      atom_selection_2 = %s
      symmetry_operation = None
      distance_ideal = %f
      sigma = 0.1
      slack = None
      limit = -0.1
      top_out = False
    }
    angle {
      atom_selection_1 = %s
      atom_selection_2 = %s
      atom_selection_3 = %s
      angle_ideal = %f
      sigma = 5
    }
    '''
    str_2 = '''refinement{
  geometry_restraints.edits{
    %s
  }
}
    '''

    i = 1
    sub_fin_str = 'a'
    for r in self.results:
      a1_str = "chain %s and resseq %s and name %s" % (
        r.atom_1.chain().id,
        r.atom_1.parent().parent().resid(),
        r.atom_1.name)
      a2_str = "chain %s and resseq %s and name %s" % (
        r.atom_2.chain().id,
        r.atom_2.parent().parent().resid(),
        r.atom_2.name)
      a3_str = "chain %s and resseq %s and name %s" % (
        r.atom_3.chain().id,
        r.atom_3.parent().parent().resid(),
        r.atom_3.name)
      if (use_defaul_parameters):
        d_ideal = 3.1
        angle_ideal = 165.72
      else:
        d_ideal = r.d
        angle_ideal = r.angle_312
      i = i + 1
      bond_angle_str = str_1 % (a1_str,a2_str,d_ideal,
                                a1_str,a2_str,a3_str,angle_ideal)
      sub_fin_str = sub_fin_str + bond_angle_str
    s_f_str = sub_fin_str[1:]
    str_final = str_2 % (s_f_str)
    file_name = pdb_file_name[:4] + ".eff"
    with open(file_name,'w') as fileobject:
      fileobject.write(str_final)




class get_hydrogen_bonds(object):
  def __init__(self,model):
    self.model = model
    self.results = self.get_hydrogen_bonds_pairs()

  def get_hydrogen_bonds_pairs(self, ideal_angel_YAD = 147.15,
                    angle_AHD_cutoff = 120,eps_angle_AHD = 30,
                    angle_HAY_min = 90,angle_HAY_max = 180,
                    eps_angle_HAY = 10,ideal_dist_A_D = 2.90,
                    sigma_for_angle= 5.0, sigma_for_bond = 0.1,
                                           eps_dist_A_D= 0.5 ):
      # Hydrogen bond  model : Y-A...H-D ;
      geometry = self.model.get_restraints_manager()
      bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                  sites_cart=self.model.get_sites_cart())
      bps_dict = {}
      [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
      hierarchy = self.model.get_hierarchy()
      atom_H = []
      atom_A = []
      atom_D = []
      atom_Y = []
      ress    = []
      results = []
      Accepter_H_pair = ["O","N","S","F","CL"]
      for a in hierarchy.atoms():
        e = a.element.strip().upper()
        if a.element_is_hydrogen():
          atom_H.append(a)
        if e == "O":
          if a.parent().resname == "HOH": continue
          atom_A.append(a)
        if e in Accepter_H_pair:
          if a.parent().resname == "HOH": continue
          atom_D.append(a)
        if e == "C":
          atom_Y.append(a)

      for a_A in atom_A:
        res = None
        diff_best = 1.e+9
        for a_D in atom_D:
          for a_H in atom_H:
            resid_A = a_A.parent().parent().resid()
            resid_D = a_D.parent().parent().resid()
            diff_r_r = abs(int(resid_A) - int(resid_D) )
            if diff_r_r < 2 :continue
            if (not a_H.is_in_same_conformer_as(a_A)): continue
            if (not is_bonded(a_H, a_D, bps_dict)): continue
            if (is_bonded(a_H, a_A, bps_dict)): continue
            if (a_A.parent().parent().resseq ==
                a_D.parent().parent().resseq): continue
            d_A_D = a_D.distance(a_A)
            if (ideal_dist_A_D - eps_dist_A_D  <
                  d_A_D  < ideal_dist_A_D + eps_dist_A_D ):
              angle_AHD = a_H.angle(a_A, a_D, deg=True)
              if (angle_AHD_cutoff - eps_angle_AHD < angle_AHD):
                diff = abs( ideal_dist_A_D - d_A_D )
                if (diff < diff_best):
                    diff_best = diff
                    res = group_args(
                      a_H=a_H,
                      a_A=a_A,
                      a_D=a_D,
                      d_A_D=d_A_D,
                      angle_AHD=angle_AHD,
                      sigma_for_bond=sigma_for_bond,
                      sigma_for_angle=sigma_for_angle,
                      ideal_angel_YAD=ideal_angel_YAD,
                      ideal_dist_A_D=ideal_dist_A_D)
                   # print sigma_for_angle, type(sigma_for_angle)
        if (res in ress):continue
        if (res is not None): ress.append(res)
      for r in ress:
        #print r.sigma_for_angle,type(r.sigma_for_angle)
        a_H = r.a_H
        a_A = r.a_A
        a_D = r.a_D
        d_A_D  = r.d_A_D
        angle_AHD = r.angle_AHD
        sigma_for_angle = r.sigma_for_angle,
        sigma_for_bond = r.sigma_for_bond,
        ideal_angel_YAD = r.ideal_angel_YAD
        ideal_dist_A_D = r.ideal_dist_A_D
        #print r.sigma_for_angle, type(r.sigma_for_angle)
        result = None
        for a_Y in atom_Y :
          if (not is_bonded(a_A, a_Y, bps_dict)): continue
          angle_YAD = a_A.angle(a_Y,a_D,deg=True)
          angle_HAY = a_A.angle(a_H,a_Y,deg=True)
          if (angle_HAY_min - eps_angle_HAY < angle_HAY <
                            angle_HAY_max + eps_angle_HAY):
            #print r.sigma_for_angle, type(r.sigma_for_angle)
            result = group_args(
              a_H=a_H,
              a_A=a_A,
              a_D=a_D,
              a_Y=a_Y,
              d_A_D=d_A_D,
              angle_HAY=angle_HAY,
              angle_AHD=angle_AHD,
              angle_YAD=angle_YAD,
              sigma_for_bond=sigma_for_bond,
              sigma_for_angle=sigma_for_angle,
              ideal_angel_YAD=ideal_angel_YAD,
              ideal_dist_A_D=ideal_dist_A_D)
          if (result in results):continue
          if (result is not None): results.append(result)

      # just keep the more possiable situation for N atom

      for i, ri in enumerate(results):
        for j, rj in enumerate(results):
          if (j <= i): continue
          a_A_i = ri.a_A
          a_D_i = ri.a_D
          a_A_j = rj.a_A
          a_D_j = rj.a_D
          if ri == rj:
            results.remove(rj)
          if a_D_i == a_D_j:
            di = a_A_i.distance(a_D_i)
            dj = a_A_j.distance(a_D_j)
            if di < dj:
              if rj in results:
                results.remove(rj)
            else:
              if ri in results:
                results.remove(ri)
      return results


  def write_restrains_file(self, pdb_file_name,
                           for_phenix_refine=True,
                           use_defaul_parameters=True):
    str_1 = '''bond{
      atom_selection_1 = %s
      atom_selection_2 = %s
      symmetry_operation = None
      distance_ideal = %f
      sigma = %f
      slack = None
      limit = -0.1
      top_out = False
    }
    angle {
      atom_selection_1 = %s
      atom_selection_2 = %s
      atom_selection_3 = %s
      angle_ideal = %f
      sigma = %f
    }
    '''

    str_2 = '''refinement{
  geometry_restraints.edits{
    %s
  }
}
    '''
    i = 1
    sub_fin_str = 'a'
    for r in self.results:
      a1_str = "chain %s and resseq %s and name %s" % (
        r.a_A.chain().id,
        r.a_A.parent().parent().resid(),
        r.a_A.name)
      a2_str = "chain %s and resseq %s and name %s" % (
        r.a_D.chain().id,
        r.a_D.parent().parent().resid(),
        r.a_D.name)
      a3_str = "chain %s and resseq %s and name %s" % (
        r.a_Y.chain().id,
        r.a_Y.parent().parent().resid(),
        r.a_Y.name)
      if (use_defaul_parameters):
        d_ideal = r.ideal_dist_A_D
      else:
        d_ideal = r.d_A_D
      i = i + 1
      if (use_defaul_parameters):
        angle_ideal = r.ideal_angel_YAD
      else:angle_ideal = r.angle_YAD
      sigma_angle = r.sigma_for_angle[0]
      sigma_bond = r.sigma_for_bond[0]
      bond_angle_str = str_1 % (a1_str, a2_str, d_ideal,
                                sigma_bond,a1_str, a2_str,
                                a3_str, angle_ideal,
                                sigma_angle)
      sub_fin_str = sub_fin_str + bond_angle_str
    s_f_str = sub_fin_str[1:]
    str_final = str_2 % (s_f_str)
    file_name = pdb_file_name
    with open(file_name,'w') as fileobject:
      fileobject.write(str_final)


class get_salt_bridge(object):
  def __init__(self, model):
    self.model = model
    self.results = self.get_salt_bridge_pairs()

  def get_salt_bridge_pairs(self, min = 1.7, max = 2.2,
                     eps1 = 0.15, eps2 = 0.8, shutoff = 4 ):
    geometry = self.model.get_restraints_manager()
    hierarchy = self.model.get_hierarchy()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                               sites_cart=self.model.get_sites_cart())
    bps_dict = {}
    [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
    positive_residues = ["ARG", "HIS", "LYS"]
    negative_residues = ["ASP", "GLU", "HIS"]
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
        if a.parent().resname == "HOH":continue
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
              result = group_args(atom_1=a1,
                                  atom_2=a2,
                                  atom_3=a3,
                                  d=d_32,
                                  angle=angle_312)

      if (result is not None): results.append(result)

    return results

  def write_restrains_file(self, pdb_file_name,
                           for_phenix_refine=True,
                           use_defaul_parameters=True):
    str_1 = '''bond{
        atom_selection_1 = %s
        atom_selection_2 = %s
        symmetry_operation = None
        distance_ideal = %f
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      angle {
        atom_selection_1 = %s
        atom_selection_2 = %s
        atom_selection_3 = %s
        angle_ideal = %f
        sigma = 5
      }
      '''
    str_2 = '''refinement{
    geometry_restraints.edits{
      %s
    }
  }
      '''
    sub_fin_str = 'a'
    for r in self.results:
      a1_str = "chain %s and resseq %s and name %s" % (
        r.atom_1.chain().id,
        r.atom_1.parent().parent().resid(),
        r.atom_1.name)
      a2_str = "chain %s and resseq %s and name %s" % (
        r.atom_2.chain().id,
        r.atom_2.parent().parent().resid(),
        r.atom_2.name)
      a3_str = "chain %s and resseq %s and name %s" % (
        r.atom_3.chain().id,
        r.atom_3.parent().parent().resid(),
        r.atom_3.name)
      if (use_defaul_parameters):
        d_ideal = 2.98
        angle_ideal = 165.07
      else:
        d_ideal = r.d
        angle_ideal = r.angle
      bond_angle_str = str_1 % (a1_str, a2_str, d_ideal, a1_str,
                                     a2_str, a3_str, angle_ideal)
      sub_fin_str = sub_fin_str + bond_angle_str
    s_f_str = sub_fin_str[1:]
    str_final = str_2 % (s_f_str)
    file_name = pdb_file_name[:4] + ".eff"
    with open(file_name, 'w') as fileobject:
      fileobject.write(str_final)


def get_stacking_system(model, dist_cutoff_1=6.0,dist_cutoff_2=3.0,
                     dist_h_cutoff=2.2,dist_v_cutoff=1.1,T_angle = 90,
                     P_angle_1 = 0,eps_angle = 26,P_angle_2 =180):

  results = []
  planes = []
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
      # calculation the distance from one center point of
      # one plane to another plane
      x, y, z = cmj
      den = math.sqrt(ai ** 2 + bi ** 2 + ci ** 2)
      dist_point_plane = abs(ai * x + bi * y + ci * z + di) / den
      # calculation the dist_h_d between two center of plane
      # in horizontal direction
      # for T type dist_point_plane
      # also represent one of dist_h_d
      dist_h_d = math.sqrt(dist**2 - dist_point_plane**2)

      """
      third limition : Dihedral angle
      for T(perpendicular)is 180
      for P(In parallel) is 90
      25 degree angular offset
      https://math.tutorvista.com/geometry/angle-between-two-planes.html
      """
      cos_a = abs(ai*aj+bi*bj+ci*cj)/(ai**2+bi**2+ci**2)**0.5/\
        (aj**2+bj**2+cj**2)**0.5
      angle = math.acos(cos_a)*180/math.pi
      # for T type,the angle should be near to 90,but 30 deviation is ok
      if (P_angle_1 + eps_angle) < angle < (T_angle - eps_angle):continue
      if (T_angle - eps_angle) < angle < (T_angle + eps_angle):
        if ((dist_point_plane > dist_v_cutoff )
            and(dist_h_d > dist_v_cutoff)):continue
      # for P type,the angle shuld be near to 180,but 30 deviaton is ok
      if (P_angle_2 - eps_angle) > angle > (T_angle + eps_angle):continue
      if (eps_angle > angle > P_angle_1)\
        or( P_angle_2 > angle > P_angle_2 - eps_angle):
        if (dist_h_d > dist_h_cutoff):continue
        if (dist_point_plane < dist_cutoff_2): continue
      print angle, dist,dist_point_plane,dist_h_d
      result = group_args( pi  = pi,
                           pj  = pj
                           )
      if result in results:continue
      if (result is not None): results.append(result)
  return results





