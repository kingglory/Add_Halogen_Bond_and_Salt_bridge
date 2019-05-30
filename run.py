from __future__ import division
import scitbx
from libtbx import group_args
import math, time
import scitbx.matrix
import iotbx.pdb
import mmtbx.model
from scitbx.array_family import flex
from cctbx import uctbx
from cctbx import crystal
from libtbx.utils import null_out
import mmtbx
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from mmtbx import monomer_library
from libtbx import group_args
import math
import iotbx.pdb.utils


import boost.python
ext = boost.python.import_ext("cctbx_geometry_restraints_ext")

def get_pair_generator(crystal_symmetry, buffer_thickness, sites_cart):
  sst = crystal_symmetry.special_position_settings().site_symmetry_table(
    sites_cart = sites_cart)
  from cctbx import crystal
  conn_asu_mappings = crystal_symmetry.special_position_settings().\
    asu_mappings(buffer_thickness = buffer_thickness)
  conn_asu_mappings.process_sites_cart(
    original_sites      = sites_cart,
    site_symmetry_table = sst)
  conn_pair_asu_table = crystal.pair_asu_table(asu_mappings = conn_asu_mappings)
  conn_pair_asu_table.add_all_pairs(distance_cutoff = buffer_thickness)
  pair_generator = crystal.neighbors_fast_pair_generator(
    conn_asu_mappings, distance_cutoff = buffer_thickness)
  return group_args(
    pair_generator    = pair_generator,
    conn_asu_mappings = conn_asu_mappings)

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



class get_hydrogen_bonds(object):
  def __init__(self,model):
    self.model = model
    self.results = self.get_hydrogen_bonds_pairs()

  def get_hydrogen_bonds_pairs(self, ideal_angle_ADCA = 111.03,
                    ideal_angle_ADC=123.72,ideal_angle_YAD = 147.15,
                    angle_AHD_cutoff = 120,eps_angle_AHD = 30,
                    angle_YAH_min = 90,angle_YAH_max = 180,
                    eps_angle_YAH = 10,ideal_dist_A_D = 2.90,
                    sigma_for_angle = 5.0, sigma_for_bond = 0.1,
                    ideal_angle_AHD = 153.30,eps_dist_A_D= 0.5,
                    ideal_angle_YAH = 120, max_cutoff   = 4.0,
                    min_cutoff   = 1.5,protein_only = True):
      # Hydrogen bond  model : Y-A...H-D-C/CA ;
      # The define of angle and bond is all from left to right in atoms order
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
      atom_C = []
      atom_CA = []
      ress    = []
      resus   = []
      results = []
      hd = ["H", "D"]
      acceptors = ["O", "N", "S", "F", "CL"]
      r = {}
      atoms = list(self.model.get_hierarchy().atoms())
      sites_cart = self.model.get_sites_cart()
      crystal_symmetry = self.model.crystal_symmetry()
      pg = get_pair_generator(
        crystal_symmetry=crystal_symmetry,
        buffer_thickness=max_cutoff,
        sites_cart=sites_cart)
      get_class = iotbx.pdb.common_residue_names_get_class
      for p in pg.pair_generator:
        i, j = p.i_seq, p.j_seq
        ei, ej = atoms[i].element, atoms[j].element
        altloc_i = atoms[i].parent().altloc
        altloc_j = atoms[j].parent().altloc
        resseq_i = atoms[i].parent().parent().resseq
        resseq_j = atoms[j].parent().parent().resseq
        # pre-screen candidates begin
        one_is_hd = ei in hd or ej in hd
        other_is_acceptor = ei in acceptors or ej in acceptors
        dist = math.sqrt(p.dist_sq)
        assert dist <= max_cutoff
        is_candidate = one_is_hd and other_is_acceptor and dist >= min_cutoff and \
                       altloc_i == altloc_j and resseq_i != resseq_j
        if (protein_only):
          for it in [i, j]:
            resname = atoms[it].parent().resname
            is_candidate &= get_class(name=resname) == "common_amino_acid"
        if (not is_candidate): continue
        rt_mx_i = pg.conn_asu_mappings.get_rt_mx_i(p)
        rt_mx_j = pg.conn_asu_mappings.get_rt_mx_j(p)
        rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
        if str(rt_mx_ji) == "x,y,z": continue
        r.setdefault(p.j_seq, []).append(rt_mx_ji)
      for k, v in zip(r.keys(), r.values()):  # remove duplicates!
        r[k] = list(set(v))

      fm = crystal_symmetry.unit_cell().fractionalization_matrix()
      om = crystal_symmetry.unit_cell().orthogonalization_matrix()
      selection = r.keys()
      for rg in hierarchy.residue_groups():
        keep = False
        for i in rg.atoms().extract_i_seq():
          if (i in selection):
            keep = True
            break
        if (keep):
          ops = r[i]
          for op in ops:
            rg_ = rg.detached_copy()
            xyz = rg_.atoms().extract_xyz()
            new_xyz = flex.vec3_double()
            for xyz_ in xyz:
              t1 = fm * flex.vec3_double([xyz_])
              t2 = op * t1[0]
              t3 = om * flex.vec3_double([t2])
              new_xyz.append(t3[0])
            rg_.atoms().set_xyz(new_xyz)
            rg_.link_to_previous = True




















      for a in hierarchy.atoms():
        e = a.element.strip().upper()
        n = a.name.strip().upper()
        if a.element_is_hydrogen():
          atom_H.append(a)
        if e == "O":
          if a.parent().resname == "HOH": continue
          atom_A.append(a)
        if e in acceptors:
          if a.parent().resname == "HOH": continue
          atom_D.append(a)
        if e == "C":
          atom_Y.append(a)
          if n == "C":
            atom_C.append(a)
          if n == "CA":
            atom_CA.append(a)


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
            d_A_H = a_A.distance(a_H)
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
                      d_A_H=d_A_H,
                      angle_AHD=angle_AHD
                    )
        if (res in ress):continue
        if (res is not None): ress.append(res)
      for r in ress:
        a_H = r.a_H
        a_A = r.a_A
        a_D = r.a_D
        d_A_D = r.d_A_D
        d_A_H = r.d_A_H
        angle_AHD = r.angle_AHD
        resu = None
        for a_Y in atom_Y :
          if (not is_bonded(a_A, a_Y, bps_dict)): continue
          angle_YAD = a_A.angle(a_Y,a_D,deg=True)
          angle_YAH = a_A.angle(a_H,a_Y,deg=True)
          if (angle_YAH_min - eps_angle_YAH < angle_YAH <
                            angle_YAH_max + eps_angle_YAH):
            resu = group_args(
              a_H=a_H,
              a_A=a_A,
              a_D=a_D,
              a_Y=a_Y,
              d_A_D=d_A_D,
              d_A_H=d_A_H,
              angle_YAH=angle_YAH,
              angle_AHD=angle_AHD,
              angle_YAD=angle_YAD,
            )
          if (resu in resus):continue
          if (resu is not None): resus.append(resu)
      for r in resus:
        a_H = r.a_H
        a_A = r.a_A
        a_D = r.a_D
        a_Y = r.a_Y
        d_A_D = r.d_A_D
        d_A_H = r.d_A_H
        angle_AHD = r.angle_AHD
        angle_YAH = r.angle_YAH
        angle_YAD = r.angle_YAD
        result = None
        for a_C in atom_C:
          if (not is_bonded(a_C, a_D, bps_dict)): continue
          for a_CA in atom_CA :
            if (not is_bonded(a_CA, a_D, bps_dict)): continue
            angle_ADC = a_D.angle(a_C,a_A,deg=True)
            angle_ADCA = a_D.angle(a_CA,a_A,deg=True)
            result = group_args(
              a_H=a_H,
              a_A=a_A,
              a_D=a_D,
              a_Y=a_Y,
              a_C=a_C,
              a_CA=a_CA,
              d_A_D=d_A_D,
              d_A_H=d_A_H,
              angle_YAH=angle_YAH,
              angle_AHD=angle_AHD,
              angle_YAD=angle_YAD,
              angle_ADC=angle_ADC,
              angle_ADCA=angle_ADCA,
              ideal_dist_A_D=ideal_dist_A_D,
              sigma_for_bond=sigma_for_bond,
              sigma_for_angle=sigma_for_angle,
              ideal_angle_YAH=ideal_angle_YAH,
              ideal_angle_YAD=ideal_angle_YAD,
              ideal_angle_ADC=ideal_angle_ADC,
              ideal_angle_AHD=ideal_angle_AHD,
              ideal_angle_ADCA=ideal_angle_ADCA
            )
            if (result in results): continue
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
    str_bond = '''bond{
      atom_selection_1 = %s
      atom_selection_2 = %s
      symmetry_operation = None
      distance_ideal = %f
      sigma = %f
      slack = None
      limit = -0.1
      top_out = False
    }
    '''
    str_angle='''angle {
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
      a_A_str = "chain %s and resseq %s and name %s" % (
        r.a_A.chain().id,
        r.a_A.parent().parent().resid(),
        r.a_A.name)
      a_D_str = "chain %s and resseq %s and name %s" % (
        r.a_D.chain().id,
        r.a_D.parent().parent().resid(),
        r.a_D.name)
      a_Y_str = "chain %s and resseq %s and name %s" % (
        r.a_Y.chain().id,
        r.a_Y.parent().parent().resid(),
        r.a_Y.name)
      a_H_str = "chain %s and resseq %s and name %s" % (
        r.a_H.chain().id,
        r.a_H.parent().parent().resid(),
        r.a_H.name)
      a_C_str = "chain %s and resseq %s and name %s" % (
        r.a_C.chain().id,
        r.a_C.parent().parent().resid(),
        r.a_C.name)
      a_CA_str = "chain %s and resseq %s and name %s" % (
        r.a_CA.chain().id,
        r.a_CA.parent().parent().resid(),
        r.a_CA.name)
      if (use_defaul_parameters):
        d_ideal = r.ideal_dist_A_D
      else:
        d_ideal = r.d_A_D
      i = i + 1
      if (use_defaul_parameters):
        angle_YAD_ideal = r.ideal_angle_YAD
        angle_AHD_ideal = r.ideal_angle_AHD
        angle_YAH_ideal = r.ideal_angle_YAH
        angle_ADC_ideal = r.ideal_angle_ADC
        angle_ADCA_ideal = r.ideal_angle_ADCA
      else:
        angle_YAD_ideal = r.angle_YAD
        angle_AHD_ideal = r.angle_AHD
        angle_YAH_ideal = r.angle_YAH
        angle_ADC_ideal = r.angle_ADC
        angle_ADCA_ideal = r.angle_ADCA
      sigma_angle = r.sigma_for_angle
      sigma_bond = r.sigma_for_bond
      bond_str = str_bond % (a_A_str, a_D_str, d_ideal,sigma_bond)
      angle_YAD = str_angle % (a_Y_str, a_A_str, a_D_str,
                               angle_YAD_ideal, sigma_angle)
      angle_AHD = str_angle % (a_A_str, a_H_str, a_D_str,
                               angle_AHD_ideal, sigma_angle)
      angle_YAH = str_angle % (a_Y_str, a_A_str, a_H_str,
                               angle_YAH_ideal, sigma_angle)
      angle_ADC = str_angle % (a_A_str, a_D_str, a_C_str,
                               angle_ADC_ideal, sigma_angle)
      angle_ADCA = str_angle % (a_A_str, a_D_str, a_CA_str,
                               angle_ADCA_ideal, sigma_angle)

      sub_fin_str = ( sub_fin_str + bond_str + angle_YAD +
        angle_AHD + angle_YAH + angle_ADC + angle_ADCA)
    s_f_str = sub_fin_str[1:]
    str_final = str_2 % (s_f_str)
    file_name = pdb_file_name
    with open(file_name,'w') as fileobject:
      fileobject.write(str_final)



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





