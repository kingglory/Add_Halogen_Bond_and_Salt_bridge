from __future__ import division
from run import define_pi_system
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import iotbx.cif
import os
import time
from libtbx import easy_run
# prepare the cif file if the pdb file needs
def prepare_cif_for_pdb_file():
  pi_file = ["1c14.pdb","3az9.pdb"]
  for pdb_file in pi_file:
    easy_run.call(" phenix.ready_set %s " %pdb_file)


# prepare the cif file if the pdb file needs
def list_cif_and_pdb_file():
 pi_file = ["1c14.pdb","3az9.pdb"]
 i = 0
 for pdb_file in pi_file:
  pdb_file = pdb_file[0:4] + ".updated.pdb"
  pdb_cif = pdb_file[0:4] + ".ligands.cif"
  if os.path.exists(pdb_cif):
    pi_file[i] = [pdb_file, pdb_cif]
  else:
    pi_file[i] = [pdb_file,None]
  i = i + 1
 return  pi_file

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
  files = [["1p5j.pdb",None],
           ["2ovp.pdb",None],
           ["3cuk.pdb","3cuk.ligands.cif"]
          ]
  pi_sites = [
       ([215, 216, 217, 218], [235, 236, 237, 238]),
       ([865, 866, 867, 868, 869], [926, 927, 928, 929]),
       ([1002, 1003, 1004, 1005], [1038, 1039, 1040, 1041]),
       ([1889, 1890, 1891, 1892], [4012, 4013, 4014, 4015]),
       ([2028, 2029, 2030, 2031, 2032, 2033], [2188, 2189, 2190, 2191]),
       ([2432, 2433, 2434, 2435], [2551, 2552, 2553, 2554]),
       ([388, 389, 390, 391, 392, 393, 394, 395, 396, 397], [521, 522, 523, 524, 525, 526, 527]),
       ([418, 419, 420, 421, 422, 423, 424, 425], [2451, 2452, 2453, 2454, 2455, 2456, 2457, 2458]),
       ([1508, 1509, 1510, 1511], [2262, 2263, 2264, 2265, 2266]),
       ([1615, 1616, 1617, 1618], [2164, 2165, 2166, 2167, 2168]),
       ([1985, 1986, 1987, 1988], [2010, 2011, 2012, 2013]),
       ([2781, 2782, 2783, 2784], [3856, 3857, 3858, 3859, 3860]),
       ([3026, 3027, 3028, 3029, 3030, 3031, 3032, 3033, 3034, 3035], [3159, 3160, 3161, 3162, 3163, 3164, 3165]),
       ([4146, 4147, 4148, 4149], [4900, 4901, 4902, 4903, 4904]),
       ([4656, 4657, 4658, 4659, 4660, 4661], [4832, 4833, 4834, 4835]),
       ([6419, 6420, 6421, 6422, 6423, 6424, 6425, 6426, 6427, 6428], [6459, 6460, 6461, 6462, 6463]),
       ([6765, 6766, 6767, 6768], [7611, 7612, 7613, 7614]),
       ([7294, 7295, 7296, 7297, 7298, 7299], [7470, 7471, 7472, 7473]),
       ([7928, 7929, 7930, 7931, 7932], [8143, 8144, 8145, 8146]),
       ([8302, 8303, 8304, 8305, 8306, 8307, 8308, 8309, 8310, 8311], [8435, 8436, 8437, 8438, 8439, 8440, 8441]),
       ([9422, 9423, 9424, 9425], [10176, 10177, 10178, 10179, 10180]),
       ([9826, 9827, 9828, 9829], [9924, 9925, 9926, 9927]),
       ([9932, 9933, 9934, 9935, 9936, 9937], [10108, 10109, 10110, 10111])
      ]


  for (pdb_file_name, cif_file_name) in files:
    print (pdb_file_name, "-"*50)
    model = get_model(pdb_file_name=pdb_file_name, cif_file_name=cif_file_name)
    result = define_pi_system(model = model)
    if result is not None:
     for r in result:
      print (r.p_i, r.p_j)
      assert (r.p_i, r.p_j) in pi_sites

if __name__ == '__main__':
    t0 = time.time()
    exercise()
    print "Time: %6.2f" % (time.time() - t0)
