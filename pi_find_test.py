from __future__ import division
from run import define_pi_system
from run import find_pi
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import iotbx.cif
import os
import time
from libtbx import easy_run
# prepare the cif file if the pdb file needs
def prepare_cif_for_pdb_file():
  pi_file = ["1c14.pdb","3az9.pdb"]
  for pdb_file in pi_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 pi_file = ["1c14.pdb","3az9.pdb"]
 i = 0
 for pdb_file in pi_file:
  pdb_file = pdb_file[0:4] + ".updated.pdb"
  pdb_cif = pdb_file[0:4] + ".ligands.cif"
  if os.path.exists(pdb_cif):
    pi_file[i] = [pdb_file, pdb_cif]
  else:
    pi_file[i] = [pdb_file,None]
  i = i + 1
 return  pi_file

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
  files = [
           ["2vb3.pdb",None],
           ["4bfq.pdb","4bfq.ligands.cif"],
           #["5yyf,pdb","5yyf.ligands.cif"],
           ["6mil.pdb","6mil.ligands.cif"],
           ["6mim.pdb","6mim.ligands.cif"],
           ["6mio.pdb",None]
          ]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    result = define_pi_system(model = model)
    if result is not None:
     for r in result:
      print ("%4.2f"%r.d_14, r.atom_1.id_str(), r.atom_4.id_str())

if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)