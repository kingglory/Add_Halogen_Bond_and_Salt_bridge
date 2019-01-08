from __future__ import division
from run import find_hydrogen_bonds
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
    model_input       = pdb_inp,
    restraint_objects = restraint_objects,
    process_input     = True,
    log               = null_out())
  return model

def exercise():
  files = [["1f23h.pdb",None],
           ["1jmhh.pdb",None],
           ["1l04h.pdb",None],
           ["1wqmh.pdb",None],
           ["1wqoh.pdb",None],
           ["1wqph.pdb",None]]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    results = find_hydrogen_bonds(model = model)
    for r in results:
      print ("%4.2f"%r.d_12, r.angle_312,r.atom_1.id_str(), r.atom_2.id_str(),
        r.atom_3.id_str())

if __name__ == '__main__':
    start = time.time()
    exercise()
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost