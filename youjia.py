from __future__ import division
import iotbx.pdb
import iotbx.cif
import scitbx
from scitbx.array_family import flex
import numpy as np
import mmtbx.hydrogens.build_hydrogens
from libtbx import group_args
from libtbx import easy_run
import math, time
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import iotbx.cif
import time

def is_bonded(atom_1, atom_2, bps_dict):
  i12 = [atom_1.i_seq, atom_2.i_seq]
  i12.sort()
  if(not tuple(i12) in bps_dict): return False
  else: return True


def get_model(pdb_file_name, cif_file_name):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
  restraint_objects = None
  if(cif_file_name is not None):
    cif_object = iotbx.cif.reader(cif_file_name).model()
    restraint_objects = [(cif_file_name, cif_object)]
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    restraint_objects = restraint_objects,
    process_input = True,
    log = null_out())
  return model

def get_CNCO_bond_angle(model):
  geometry = model.get_restraints_manager()
  bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
    sites_cart=model.get_sites_cart())
  bps_dict = {}
  [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
  hierarchy = model.get_hierarchy()
  atom1s = []
  atom2s = []
  atom3s = []
  results = []
  for a in hierarchy.atoms():
    if (not a.parent().resname == "MTN"):continue
    if a.element.strip().upper() == "C":
      atom1s.append(a)
    if a.element.strip().upper() == "N":
      atom2s.append(a)
    if a.element.strip().upper() == "O":
      atom3s.append(a)
  for i,a1 in enumerate(atom1s):
    for a2 in atom2s:
      if (not a1.is_in_same_conformer_as(a2)): continue
      if (not is_bonded(a1, a2, bps_dict)): continue
      for j,a4 in enumerate(atom1s):
        if (j <= i): continue
        if (not a4.is_in_same_conformer_as(a2)): continue
        if (not is_bonded(a4, a2, bps_dict)): continue
        if (not a1.is_in_same_conformer_as(a4)): continue
        if (is_bonded(a1, a4, bps_dict)): continue
        for a3 in atom3s:
          if (not a2.is_in_same_conformer_as(a3)): continue
          if (not is_bonded(a3, a2, bps_dict)): continue
          d12 = a1.distance(a2)
          d23 = a2.distance(a3)
          d24 = a2.distance(a4)
          angle123 = a2.angle(a1,a3,deg = True)
          angle124 = a2.angle(a1,a4,deg = True)
          angle324 = a2.angle(a3,a4,deg = True)
          result = group_args(atom_1 = a1,
                              atom_2 = a2,
                              atom_3 = a3,
                              atom_4 = a4,
                              d12 = d12,
                              d23 = d23,
                              d24 = d24,
                              a123 = angle123,
                              a124 = angle124,
                              a324 = angle324)
          if (result is not None): results.append(result)
  return results

def exercise():
  pdb_file_name = raw_input("pdb file name here:")
  easy_run.call("phenix.fetch_pdb {0}".format(pdb_file_name[0:4]))
  easy_run.call("phenix.ready_set {0}".format(pdb_file_name))
  cif_file_name = pdb_file_name[0:4] + ".ligands.cif"
  model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
  results = get_CNCO_bond_angle(model)
  for r in results:
    print ("%s" %r.atom_1.id_str(), r.atom_2.id_str(),r.atom_3.id_str(), r.atom_4.id_str(),
             r.d12, r.d23, r.d24, r.a123, r.a124, r.a324)



if __name__ == '__main__':
    start = time.time()
    exercise()
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost