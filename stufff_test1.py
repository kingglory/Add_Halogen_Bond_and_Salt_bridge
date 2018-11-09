
from __future__ import division
import iotbx.pdb
from iotbx.pdb import hierarchy
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
import os


pdb_file = "5v7d.pdb"
pdb_inp = iotbx.pdb.input(file_name=pdb_file)
model = mmtbx.model.manager(
        model_input=pdb_inp,
        process_input=True,
       log=null_out())
hierarchy = model.get_hierarchy()
for atom in hierarchy.atoms():
    e = atom.name.strip().upper()
    if e == "BR":
        print atom.id_str()

try:
    a=2
    c=0
    d = a/c
except BaseException,e :
    if len(e.message)>0:
        print "hello"