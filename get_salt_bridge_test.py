from __future__ import division
from run import find_salt_bridge
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import iotbx.cif
import time

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
  files = [["1jvuh.pdb","1jvuh.ligands.cif"],
           ["2x7gh.pdb","2x7gh.ligands.cif"],
           ["3alph.pdb","3alph.ligands.cif"],
           ["4x7nh.pdb","4x7nh.ligands.cif"]
           ]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    result1 = find_salt_bridge(model=model)
    for r in result1:
      print ("%s"% r.atom_1.id_str(), r.atom_2.id_str(),r.atom_3.id_str())
      # assert  XXX


if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)
