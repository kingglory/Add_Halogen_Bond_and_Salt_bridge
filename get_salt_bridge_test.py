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
  files = [["1jvuh1.pdb", "1jvuh.ligands.cif"],
           ["2x7gh2.pdb", "2x7gh.ligands.cif"],
           ["2x7gh3.pdb", "2x7gh.ligands.cif"],
           ["3alph1.pdb", "3alph.ligands.cif"],
           ["3alph2.pdb", "3alph.ligands.cif"],
           ["3alph3.pdb", "3alph.ligands.cif"],
          ]
  Ideal_Salt_Bridge_files = {"1jvuh1.pdb":[('pdb="HH21 ARG A  10 "', 'pdb=" NH2 ARG A  10 "', 'pdb=" OE1 GLU A   2 "'),
                                           ('pdb=" HE  ARG A  10 "', 'pdb=" NE  ARG A  10 "', 'pdb=" OE2 GLU A   2 "')],
                             "2x7gh2.pdb":[('pdb=" HZ3 LYS A 218 "', 'pdb=" NZ  LYS A 218 "', 'pdb=" OD2 ASP A 147 "'),
                                           ('pdb=" HZ2 LYS A 218 "', 'pdb=" NZ  LYS A 218 "', 'pdb=" OD2 ASP A 150 "')],
                             "2x7gh3.pdb":[('pdb="HH22 ARG A 685 "', 'pdb=" NH2 ARG A 685 "', 'pdb=" OE1 GLU A 565 "'),
                                           ('pdb="HH12 ARG A 685 "', 'pdb=" NH1 ARG A 685 "', 'pdb=" OE2 GLU A 565 "')],
                             "3alph1.pdb":[('pdb="HH12 ARG A  96 "', 'pdb=" NH1 ARG A  96 "', 'pdb=" OD1 ASP A 118 "'),
                                           ('pdb="HH22 ARG A  96 "', 'pdb=" NH2 ARG A  96 "', 'pdb=" OD2 ASP A 118 "')],
                             "3alph2.pdb":[('pdb=" HZ2 LYS A 166 "', 'pdb=" NZ  LYS A 166 "', 'pdb=" OE2 GLU A 192 "')],
                             "3alph3.pdb":[('pdb=" HZ2 LYS A 236 "', 'pdb=" NZ  LYS A 236 "', 'pdb=" OE1 GLU A 186 "')]}
                       
  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    results= find_salt_bridge(model=model,pdb_file_name=pdb_file_name)
    for r in results:
      Salt_Bridge_sites = Ideal_Salt_Bridge_files[pdb_file_name]
      print ("%s"% r.atom_1, r.atom_2,r.atom_3)
      assert  (r.atom_1, r.atom_2,r.atom_3) in  Salt_Bridge_sites
      #print ("%s"% r.atom_1.id_str(), r.atom_2.id_str(),r.atom_3.id_str())
      
       

if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)
