from run import distance_atomToArea
import numpy as np
import pandas as pd


if __name__ == '__main__':
  atom1_xyz = (0, 0, 1)
  atom2_xyz = (0, 1, 0)
  atom3_xyz = (0, 0, 0)
  atom4_xyz = (1, 0, 0)
  d = distance_atomToArea(atom1_xyz, atom2_xyz,atom3_xyz, atom4_xyz)
  print("The distance from one atom to a plane: " + str(d))
