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
           ["2x7gh.db","2x7gh.ligands.cif"],
           ["3alph.pdb","3alph.ligands.cif"],
           ["4x7nh.pdb","4x7nh.ligands.cif"]
           ]
  salt_bridge_sites = [
    ('pdb="HH21 ARG A  10 "', 'pdb=" NH2 ARG A  10 "', 'pdb=" OE1 GLU A   2 "'),
    ('pdb=" HE  ARG A  10 "', 'pdb=" NE  ARG A  10 "', 'pdb=" OE2 GLU A   2 "'),
    ('pdb=" HE  ARG A  85 "', 'pdb=" NE  ARG A  85 "', 'pdb=" OD1 ASP A  83 "'),
    ('pdb="HH21 ARG B  10 "', 'pdb=" NH2 ARG B  10 "', 'pdb=" OE1 GLU B   2 "'),
    ('pdb=" HE  ARG B  10 "', 'pdb=" NE  ARG B  10 "', 'pdb=" OE2 GLU B   2 "'),
    ('pdb=" HE  ARG A  91 "', 'pdb=" NE  ARG A  91 "', 'pdb=" OD1 ASP A 111 "'),
    ('pdb="HH21 ARG A  91 "', 'pdb=" NH2 ARG A  91 "', 'pdb=" OD2 ASP A 111 "'),
    ('pdb=" HZ3 LYS A 121 "', 'pdb=" NZ  LYS A 121 "', 'pdb=" OE1 GLU A 136 "'),
    ('pdb=" HZ3 LYS A 218 "', 'pdb=" NZ  LYS A 218 "', 'pdb=" OD2 ASP A 147 "'),
    ('pdb=" HZ2 LYS A 218 "', 'pdb=" NZ  LYS A 218 "', 'pdb=" OD2 ASP A 150 "'),
    ('pdb="HH22 ARG A 685 "', 'pdb=" NH2 ARG A 685 "', 'pdb=" OE1 GLU A 565 "'),
    ('pdb="HH12 ARG A 685 "', 'pdb=" NH1 ARG A 685 "', 'pdb=" OE2 GLU A 565 "'),
    ('pdb="HH12 ARG A 559 "', 'pdb=" NH1 ARG A 559 "', 'pdb=" OD1 ASP A 606 "'),
    ('pdb="HH22 ARG A  96 "', 'pdb=" NH2 ARG A  96 "', 'pdb=" OD2 ASP A 118 "'),
    ('pdb=" HZ2 LYS A 236 "', 'pdb=" NZ  LYS A 236 "', 'pdb=" OE1 GLU A 186 "'),
    ('pdb=" HZ2 LYS A 166 "', 'pdb=" NZ  LYS A 166 "', 'pdb=" OE2 GLU A 192 "'),
    ('pdb="HH12 ARG A 212 "', 'pdb=" NH1 ARG A 212 "', 'pdb=" OE2 GLU A 194 "'),
    ('pdb="HH11 ARG B 101 "', 'pdb=" NH1 ARG B 101 "', 'pdb=" OD2 ASP B 106 "'),
    ('pdb="HH22 ARG B  96 "', 'pdb=" NH2 ARG B  96 "', 'pdb=" OD1 ASP B 118 "'),
    ('pdb="HH21 ARG B  96 "', 'pdb=" NH2 ARG B  96 "', 'pdb=" OD2 ASP B 118 "'),
    ('pdb=" HZ1 LYS B 178 "', 'pdb=" NZ  LYS B 178 "', 'pdb=" OE2 GLU B 119 "'),
    ('pdb=" HZ2 LYS B 236 "', 'pdb=" NZ  LYS B 236 "', 'pdb=" OE2 GLU B 186 "'),
    ('pdb=" HZ3 LYS B 166 "', 'pdb=" NZ  LYS B 166 "', 'pdb=" OE2 GLU B 192 "'),
    ('pdb="HH21 ARG B 217 "', 'pdb=" NH2 ARG B 217 "', 'pdb=" OD1 ASP B 272 "'),
    ('pdb=" HZ1 LYS A 893 "', 'pdb=" NZ  LYS A 893 "', 'pdb=" OE1 GLU A 609 "'),
    ('pdb=" HE  ARG A 638 "', 'pdb=" NE  ARG A 638 "', 'pdb=" OE1 GLU A 634 "'),
    ('pdb="HH21 ARG A 638 "', 'pdb=" NH2 ARG A 638 "', 'pdb=" OE2 GLU A 634 "'),
    ('pdb="HH12 ARG A1065 "', 'pdb=" NH1 ARG A1065 "', 'pdb=" OE1 GLU A 994 "'),
    ('pdb="HH22 ARG A1065 "', 'pdb=" NH2 ARG A1065 "', 'pdb=" OE2 GLU A 994 "'),
    ('pdb=" HE  ARG A1027 "', 'pdb=" NE  ARG A1027 "', 'pdb=" OE1 GLU A1016 "'),
    ('pdb=" HE  ARG A 911 "', 'pdb=" NE  ARG A 911 "', 'pdb=" OE2 GLU A1050 "')
  ]
                       
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    result1 = find_salt_bridge(model=model)
    for r in result1:
      print ("%s"% r.atom_1.id_str(), r.atom_2.id_str(),r.atom_3.id_str())
      assert  (r.atom_1.id_str(), r.atom_2.id_str(),
               r.atom_3.id_str()) in  salt_bridge_sites
      
       

if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)
