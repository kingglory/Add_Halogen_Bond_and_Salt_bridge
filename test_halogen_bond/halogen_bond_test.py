from __future__ import division
import sys
sys.path.append("..")
import run
from run import get_halogen_bonds
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
  files = [["5v7d.pdb",None],
           ["2h79.pdb","2h79.ligands.cif"],
           ["2ito.pdb","2ito.ligands.cif"],
           ["2oxy.pdb","2oxy.ligands.cif"],
           ["3v04.pdb","3v04.ligands.cif"]
           ]
  Halogen_atom_pairs = [
    ('pdb="BR  ABYR A  18 "',
     'pdb=" N   GLY A  28 "'),
    ('pdb="BR  ABYR A  18 "',
     'pdb=" O   GLY A  28 "'),
    ('pdb=" I1   T3 A   1 "',
     'pdb=" O   PHE A 218 "'),
    ('pdb=" I2   T3 A   1 "',
     'pdb=" O   GLY A 290 "'),
    ('pdb=" FAB IRE A2020 "',
     'pdb=" OE2 GLU A 762 "'),
    ('pdb="CL   IRE A2020 "',
     'pdb=" O   LEU A 788 "'),
    ('pdb="BR1  K17 A1001 "',
     'pdb=" O   VAL A 116 "'),
    ('pdb="BR1  K17 B1002 "',
     'pdb=" O   VAL B 116 "'),
    ('pdb=" I1  V04 A 501 "',
     'pdb=" O   VAL A 127 "')
  ]
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name,
                      cif_file_name=cif_file_name)
    get_x_bonds = get_halogen_bonds(model=model)
    results = get_x_bonds.get_halogen_bonds_pairs()
    get_x_bonds.write_restrains_file(pdb_file_name=pdb_file_name)
    for r in results:
      assert (r.atom_1.id_str(), r.atom_2.id_str()) in Halogen_atom_pairs
      print ("%s" % r.atom_1.id_str(), r.atom_2.id_str())


if __name__ == '__main__':
  start = time.time()
  exercise()
  end = time.time()
  time_cost = (end - start)
  print "it cost % seconds" % time_cost
