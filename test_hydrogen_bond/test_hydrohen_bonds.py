from __future__ import division
import sys
sys.path.append("..")
import run
from run import get_hydrogen_bonds
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
           #["m-helix.updated.pdb",None],
            ["4gif_part.pdb",None]
           #["1kych.pdb", None]
          ]



  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    #pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)

    model = get_model(pdb_file_name=pdb_file_name,
                      cif_file_name=cif_file_name,
                      )


    get_h_bonds = get_hydrogen_bonds(model=model)

    get_h_bonds.write_restrains_file(pdb_file_name=pdb_file_name[:-4]+'.eff',
                                     use_defaul_parameters=True)

    results = get_h_bonds.get_hydrogen_bonds_pairs()
    #get_h_bonds.get_hydrogen_bonds_from_pair_atoms()

    for r in results:
      print r.a_H.id_str(),r.a_A.id_str()


if __name__ == '__main__':
    start = time.time()
    exercise()
    end = time.time()
    time_cost = (end - start)
    print "it cost % seconds" % time_cost