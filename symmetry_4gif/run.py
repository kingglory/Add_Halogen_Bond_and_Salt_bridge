from __future__ import division
import iotbx.pdb
import mmtbx.model
import math
from libtbx import group_args

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

def run(model,
        max_cutoff   = 4.0, 
        min_cutoff   = 1.5,
        hd           = ["H", "D"],
        acceptors    = ["O","N","S","F","CL"],
        protein_only = True):
  atoms = list(model.get_hierarchy().atoms())
  sites_cart = model.get_sites_cart()
  crystal_symmetry = model.crystal_symmetry()
  pg = get_pair_generator(
    crystal_symmetry = crystal_symmetry, 
    buffer_thickness = max_cutoff, 
    sites_cart       = sites_cart)
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
    if(protein_only):
      for it in [i,j]:
        resname = atoms[it].parent().resname
        is_candidate &= get_class(name=resname) == "common_amino_acid"

    # pre-screen candidates end
    rt_mx_i = pg.conn_asu_mappings.get_rt_mx_i(p)
    rt_mx_j = pg.conn_asu_mappings.get_rt_mx_j(p)
    rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
    print "%5.3f"%math.sqrt(p.dist_sq), \
      "<", atoms[i].parent().resname, resseq_i, ei, atoms[i].name, ">", \
      "<", atoms[j].parent().resname, resseq_j, ej, atoms[j].name, ">", \
      rt_mx_ji, i, j

if(__name__ == "__main__"):
  pdb_inp = iotbx.pdb.input(file_name="4gif_part.pdb")
  model = mmtbx.model.manager(model_input = pdb_inp, build_grm = True)
  run(model = model)

