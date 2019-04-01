from __future__ import division
import iotbx.pdb
import iotbx.cif
import re
import scitbx
from scitbx.array_family import flex
import numpy as np
import pandas as pd
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
                              d_C_L_N = d12,
                              d_O_N = d23,
                              d_C_R_N = d24,
                              a_C_L_N_O = angle123,
                              a_C_N_C = angle124,
                              a_C_R_N_O = angle324)
          if (result is not None): results.append(result)
  return results

def exercise():
  pdb_file_names = raw_input("pdb file name here:")
  p_f = pdb_file_names.split(' ')
  for pdb_file_name in p_f:
    easy_run.call("phenix.fetch_pdb {0}".format(pdb_file_name[0:4]))
    easy_run.call("phenix.ready_set {0}".format(pdb_file_name))
    print (pdb_file_name, "-" * 50)
    cif_file_name = pdb_file_name[0:4] + ".ligands.cif"
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    results = get_CNCO_bond_angle(model)
    for r in results:
      '''print dir(r.atom_1.chain())
      MTN = {"id":["CSN values"],"O-N":[1.272],"C(L)-N":[1.482],
            "C(R)R-N":[1.482],"C-N-C":[115.232],"C(L)-N-O":[122.384],
            "C(R)-N-O":[122.384]}
      df = pd.DataFrame(MTN,columns = ["id","O-N","C(L)-N",
            "C(R)R-N","C-N-C","C(L)-N-O","C(R)-N-O"],index=["one"])
      Series = pd.Series({"O-N":r.d23,"C(L)-N":r.d12,
            "C(R)R-N":r.d24,"C-N-C":r.a124,"C(L)-N-O":r.a123,
            "C(R)-N-O":r.a324})
      data = pd.DataFrame()
      Series = [[r.d12, r.d23, r.d24, r.a123, r.a124, r.a324]]
      df = data.append(Series,ignore_index=True)

        print df
       '''

      print ("%s" %r.atom_1.id_str(), r.atom_2.id_str(),r.atom_3.id_str(),
         r.atom_4.id_str(),r.d_O_N, r.d_C_L_N, r.d_C_R_N, r.a_C_N_C, r.a_C_L_N_O, r.a_C_R_N_O)



if __name__ == '__main__':
  start = time.time()
  exercise()
  end = time.time()
  time_cost = (end - start)
  print "it cost % seconds" % time_cost