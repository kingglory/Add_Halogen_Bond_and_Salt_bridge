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
    model_input = pdb_inp,
    restraint_objects = restraint_objects,
    process_input = True,
    log = null_out())
  return model

def exercise():
  files = [
           ["1aieh.pdb",None],
           ["1kych.pdb",None],
           ["3q8jh.pdb",None],
           ["4gifh.pdb",None],
           ["5c11h.pdb",None],
           ["2ona_complete_occH0_qm_refined.pdb",None]
           ]
  Hydrogen_atom_pairs = [
    ('pdb=" HD3BARG A 342 "', 'pdb=" OE2 GLU A 346 "'),
    ('pdb=" HG3 LYS A 351 "', 'pdb=" OXT GLY A 356 "'),
    ('pdb="HH22 ARG A  13 "', 'pdb=" OE1 GLU A   9 "'),
    ('pdb="HH21 ARG A  15 "', 'pdb=" O3  SO4 A 101 "'),
    ('pdb=" H   ARG A   4 "', 'pdb=" O1  SIN A   0 "'),
    ('pdb=" H   PHE A   5 "', 'pdb=" OE1 GLU A   8 "'),
    ('pdb=" H   PHE A  14 "', 'pdb=" OD1 ASN A  11 "'),
    ('pdb="HD21 ASN A  11 "', 'pdb=" OD1 ASP A  31 "'),
    ('pdb=" H   SER A 703 "', 'pdb=" OE2 GLU A 706 "'),
    ('pdb=" HD2 ARG A 712 "', 'pdb=" OH  TYR A 708 "'),
    ('pdb="HE21 GLN B   5 "', 'pdb=" OD1 ASP A  17 "'),
    ('pdb=" H   TYR A  46 "', 'pdb=" OD2 ASP A  22 "'),
    ('pdb=" HB3 SER A  37 "', 'pdb=" OE2 GLU A  39 "'),
    ('pdb=" H2  HOH B   7 "', 'pdb=" OXT VAL B   6 "')   
  ]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    results = find_hydrogen_bonds(model=model)
    for r in results:
      #assert (r.atom_1.id_str(), r.atom_2.id_str()) in Hydrogen_atom_pairs
      print ("%s"% r.atom_1.id_str(), r.atom_2.id_str())

if __name__ == '__main__':
    start = time.time()
    exercise()
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost
