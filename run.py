from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex
import math
import mmtbx.model
from libtbx.utils import null_out

def distance(s1, s2):
  return math.sqrt(
    (s1[0]-s2[0])**2 + 
    (s1[1]-s2[1])**2 + 
    (s1[2]-s2[2])**2)

def run(hierarchy, vdwr):
  halogens = ["CL","BR","I","F"] 
  for atom_1 in hierarchy.atoms():
    e = atom_1.element.strip().upper()
    if(e in halogens):
      print e, vdwr[e]
      #for atom_2 in hierarchy.atoms():
      #  d = distance(atom_1.xyz, atom_2.xyz)
      #  print d
        
if __name__ == '__main__':
  pdb_file = "5v7d.pdb"
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  model = mmtbx.model.manager(
    model_input   = pdb_inp, 
    process_input = True,
    log           = null_out())
  vdwr = model.get_vdw_radii()
  run(hierarchy = model.get_hierarchy(),
      vdwr      = vdwr)


