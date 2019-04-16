from __future__ import division
import iotbx.pdb
import iotbx.cif
import os
import mmtbx.hydrogens.build_hydrogens
from libtbx import easy_run
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
  pdb_file_names = raw_input("pdb file name here,space interval,such as : 2ona.pdb 2whb.pdb")
  p_f = pdb_file_names.split(' ')
  for pdb_file_name in p_f:
    easy_run.call("phenix.fetch_pdb {0}".format(pdb_file_name[0:4]))
    easy_run.call("phenix.ready_set {0}".format(pdb_file_name))
    print (pdb_file_name, "-" * 50)
    cif_file_name = pdb_file_name[0:4] + ".ligands.cif"
    if os.path.exists(cif_file_name):
      cif_file_names = cif_file_name
    else :
      cif_file_names = None
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_names)
    resnames = raw_input('list of resnames ,comma interval,such as: resname PHE,resname HIS')
    resn = resnames.split(',')
    ss = " or ".join(resn)
    m_sel = model.selection(ss)
    new_model = model.select(m_sel)
    hierarchy = new_model.get_hierarchy()
    crystal_symmetry = new_model.crystal_symmetry()
    hierarchy.write_pdb_file(file_name=pdb_file_name[0:4] + "new.pdb",
                             crystal_symmetry=crystal_symmetry)



if __name__ == '__main__':
  start = time.time()
  exercise()
  end = time.time()
  time_cost = (end - start)
  print "it cost % seconds" % time_cost