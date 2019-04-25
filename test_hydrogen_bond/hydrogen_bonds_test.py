from __future__ import division
import sys
sys.path.append("..")
import run
from run import get_hydrogen_bonds_pairs
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
           ["1kych1.pdb",None],
           ["1kych2.pdb",None],
           ["1kych3.pdb",None],
           ["3q8jh1.pdb",None],
           ["3q8jh2.pdb",None],
           ["3q8jh3.pdb",None]
           ]
  Ideal_Hydrogen_Bonds_files = {"1kych1.pdb":[('pdb=" HA  ARG A   5 "', 'pdb=" OE1 GLU A   8 "'),
                                              ('pdb=" H   ARG A   4 "', 'pdb=" O   GLU A   1 "'),
                                              ('pdb=" H   GLU A   8 "', 'pdb=" O   ARG A   4 "'),
                                              ('pdb=" HB2 ARG A   5 "', 'pdb=" OE2 GLU A   1 "'),
                                              ('pdb=" H   ARG A   5 "', 'pdb=" O   GLU A   1 "'),
                                              ('pdb=" HG2 ARG A   4 "', 'pdb=" OE2 GLU A   8 "')],
                                "1kych2.pdb":[('pdb=" HA  GLU A   8 "', 'pdb=" OE1 GLU A  11 "'),
                                              ('pdb=" H   GLU A  11 "', 'pdb=" O   GLU A   8 "'),
                                              ('pdb=" H   ARG A  12 "', 'pdb=" O   GLU A   8 "'),
                                              ('pdb=" HD2 ARG A  12 "', 'pdb=" OE1 GLU A   9 "')],
                                "1kych3.pdb":[('pdb="HH21 ARG A  15 "', 'pdb=" O3  SO4 A 101 "'),
                                              ('pdb=" HE  ARG A  15 "', 'pdb=" O3  SO4 A 101 "')],
                                "3q8jh1.pdb":[('pdb=" H   PHE A   5 "', 'pdb=" OE1 GLU A   8 "'),
                                              ('pdb=" H   GLU A   8 "', 'pdb=" O   PHE A   5 "'),
                                              ('pdb=" HG2 GLU A   8 "', 'pdb=" O   SER A   9 "')],
                                "3q8jh2.pdb":[('pdb=" H   ASN A  11 "', 'pdb=" O   TYR A  15 "'),
                                              ('pdb=" H   TYR A  15 "', 'pdb=" O   ASN A  11 "'),
                                              ('pdb=" H   GLN A  13 "', 'pdb=" OD1 ASN A  11 "'),
                                              ('pdb=" HA  PHE A  14 "', 'pdb=" O   GLN A  13 "'),
                                              ('pdb=" H   PHE A  14 "', 'pdb=" OD1 ASN A  11 "')],
                                "3q8jh3.pdb":[('pdb="HD13 LEU A  23 "', 'pdb=" O   CYS A  18 "'),
                                              ('pdb=" HG  LEU A  21 "', 'pdb=" O   PRO A  19 "'),
                                              ('pdb="HD23 LEU A  21 "', 'pdb=" O   GLY A  20 "')]}
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    results = get_hydrogen_bonds_pairs(model=model)
    for r in results:
      d_12 = r.atom_1.distance(r.atom_2)
      Hydrogen_bond_pairs = Ideal_Hydrogen_Bonds_files[pdb_file_name]
      assert (r.atom_1.id_str(), r.atom_2.id_str()) in Hydrogen_bond_pairs
      print ("%s"% r.atom_1.id_str(), r.atom_2.id_str())

if __name__ == '__main__':
    start = time.time()
    exercise()
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost
