from __future__ import division
import sys
sys.path.append("..")
import run
from run import get_stacking_system
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import iotbx.cif
import os
import time
from libtbx import easy_run

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

def exercise():
  files = [  ["1p5j.pdb", None],
             ["1p5j_1.pdb", None],
             ["1p5j_2.pdb", None],
             ["2ito.pdb", "2ito.ligands.cif"],
             ["3v04.pdb", "3v04.ligands.cif"],
             ["5yyf.pdb", "5yyf.ligands.cif"],
             ["5yy2.pdb", "5yyf.ligands.cif"],
             ["5yy3.pdb", "5yyf.ligands.cif"],
          ]



  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    result = get_stacking_system(model = model)
    if result is not None:
     for r in result:
       p_i = list(r.pi.extract_i_seq()),
       p_j = list(r.pj.extract_i_seq()),
       pi_atoms_name = list(r.pi.extract_name()),
       pj_atoms_name = list(r.pj.extract_name())
       print (p_i, p_j)
      #pi_stacking_site =  pi_sites[pdb_file_name]
      #assert (r.p_i, r.p_j) in pi_stacking_site

if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)
