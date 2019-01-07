from __future__ import division
from run import find_halogen_bonds
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
  files = [["5v7d.pdb",None],
           ["2h79.pdb","2h79.ligands.cif"],
           ["2ito.pdb","2ito.ligands.cif"],
           ["2oxy.pdb","2oxy.ligands.cif"],
           ["2vag.pdb","2vag.ligands.cif"],
           #["2yj8.pdb","2yj8.ligands.cif"],
           ["3v04.pdb","3v04.ligands.cif"],
           ["4e7r.pdb","4e7r.ligands.cif"]
           ]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    results = find_halogen_bonds(model = model)
    for r in results:
      print("%4.2f"%r.d_12, r.angle_312, r.angle_214,r.atom_1.id_str(), 
        r.atom_2.id_str())

if __name__ == '__main__':
  t0 = time.time()
  exercise()
  print "Time: %6.2f"%(time.time()-t0)
