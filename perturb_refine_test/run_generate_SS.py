from __future__ import division
from mmtbx.secondary_structure import build as ssb
from scitbx.array_family import flex

def run():
  h = ssb.secondary_structure_from_sequence(ssb.alpha_helix_str,
    "GGGGGGGGGGGGGGGGGGGG")
  h.write_pdb_file(file_name="m-helix.pdb")

if (__name__ == "__main__"):
  run()
