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
  files = [["5v7dh.pdb",None],
           ["2h79h.pdb","2h79h.ligands.cif"],
           ["2itoh.pdb","2itoh.ligands.cif"],
           ["2oxyh.pdb","2oxyh.ligands.cif"],
           ["2vagh.pdb","2vagh.ligands.cif"],
           ["3v04h.pdb","3v04h.ligands.cif"],
           ["4e7rh.pdb","4e7rh.ligands.cif"]]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    final_result = find_hydrogen_bonds(model = model)
    for r in final_result:
      print ("%4.2f"%r.d_12, r.angle_312,r.angle_214,r.atom_1.id_str(), r.atom_2.id_str(),
        r.atom_3.id_str(), r.atom_4.id_str())

if __name__ == '__main__':
    start = time.time()
    exercise()
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost