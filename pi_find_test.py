from __future__ import division
from run import define_pi_system
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
  files = [  #["1p5j.pdb", None],
             ["1p51.pdb", None],
             #["1p52.pdb", None],
             ["1p53.pdb", None],
             ["1p54.pdb", None],
             #["2ovp.pdb", None],
             #["5yy1.pdb", "5yyf.ligands.cif"],
             ["5yy2.pdb", "5yyf.ligands.cif"],
             ["5yy3.pdb", "5yyf.ligands.cif"],
             #["6mip.pdb","6mip.ligands.cif"]
          ]
  pi_sites = [

      ]


  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    result = define_pi_system(model = model,pdb_file_name = pdb_file_name)
    if result is not None:
     for r in result:
      print (r.p_i, r.p_j)
     #assert (r.p_i, r.p_j) in pi_sites

if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)
