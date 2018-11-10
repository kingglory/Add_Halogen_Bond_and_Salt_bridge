
from __future__ import division
from libtbx import easy_run
import iotbx.pdb
from iotbx.pdb import hierarchy
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
import os
from libtbx import easy_run

X_bonds_file = ["5v7d.pdb","2h79.pdb", "2ito.pdb", "2oxy.pdb","2vag.pdb",
              "2yj8.pdb", "3v04.pdb", "4e7r.pdb"]
i = 0
for pdb_file in X_bonds_file:
 pdb_cif = pdb_file[0:4] + ".ligands.cif"
 if os.path.exists(pdb_cif):
# easy_run.call(" phenix.ready_set %s " %pdb_file)
  X_bonds_file[i] = [pdb_file,pdb_cif]
 else :
  X_bonds_file[i] = [pdb_file,None]
 i = i + 1
print X_bonds_file


easy_run.call("phenix.geometry_minimization 2yj8.pdb modified.pdb")



