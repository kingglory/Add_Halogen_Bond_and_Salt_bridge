# Find_Halogen_Bond_and_Salt_bridge
The project is to write a CCTBX based Python program that would take a PDB file with a protein structure and find all halogen bonds and salt bridges in the structure (meaning finding all pairs of amino-acid residues participating in said interactions). 


FILE_DECRIBE:


run.py                   is the main scripe to implement the functions of Finding halogebonds and salt bridges
get_salt_bridge_test.py  is the test file for get salt bridges function
halogen_bond_test.py     is the test file for get halogen donds function
test.sbatch              is the bash file to submit tasks to cluster
THe others               are the test related files


HALOGEN BONGS FILE TEST:


["5v7d.pdb","2h79.pdb", "2ito.pdb", "2oxy.pdb","2vag.pdb", "2yj8.pdb", "3v04.pdb", "4e7r.pdb"]



          pdb:5v7d
CE2    A  18  BYR;
BR     A  18  BYR;
O      A  28  GLY;
C      A  28  GLY;
X_bond_distance =  2.88
sum_vdwr        =  3.28
percentages     =  0.8780
angle1          =  150.78 degrees
angle2          =  115.49 degrees


          pdb:2h79
C5    A   1  T3;
I1    A   1  T3;
O     A 218  PHE;
C     A 218  PHE;
X_bond_distance =  3.05
sum_vdwr        =  3.38
percentages     =  0.9034
angle1          =  173.30 degrees
angle2          =  121.20 degrees

          pdb:2ito
 
 CAX  A2020   IRE;
 CL   A2020   IRE;
 O    A 788   LEU;
 C    A 788   LEU;
X_bond_distance =  3.29
sum_vdwr        =  3.15
percentages     =  1.04
angle1          =  163.19 degrees
angle2          =  108.71 degrees




          pdb:2oxy
C1    A1001  K17;
BR1   A1001  K17;
O     A 116  VAL;
C     A 116  VAL;
X_bond_distance =  2.95
sum_vdwr        =  3.28
percentages     =  0.9063
angle1          =  176.33 degrees
angle2          =  131.67 degrees


C1    B100２  K17;
BR1   B100２  K17;
O     B 116  VAL;
C     B 116  VAL;
X_bond_distance =  2.9２
sum_vdwr        =  3.28
percentages     =  0.８９０２
angle1          =  17３.８６ degrees
angle2          =  131.８３ degrees


          pdb:2vag
CAP    A1482   V25;
CL2    A1482   V25;
O      A 242   GLU;
C      A 242   GLU;
X_bond_distance =  2.94
sum_vdwr        =  3.15
percentages     =  0.9334
angle1          =  171.18 degrees
angle2          =  158.11 degrees







          pdb:2yj8
C27    A1221  YJ8 ; 
I30    A1221  YJ8 ; 
O      A  61  GLY ; 
C      A  61  GLY ;    
X_bond_distance =  3.12
sum_vdwr        =  3.38
percentages     =  0.9233
angle1          =  174.59 degrees
angle2          =  144.96 degrees

          pdb:3v04  (O A128, less distance ,why not???)
C11     A 501   V04;
I1      A 501   V04;
O       A 127   VAL;
C       A 127   VAL;
X_bond_distance =  3.31
sum_vdwr        =  3.38
percentages     =  0.9792
angle1          =  176.26 degrees
angle2          =  126.81 degrees



          pdb:4e7r
C37　　H 302　ONW
CL2　　H 302　ONW
O　　　H 98 　ASN
C　　　H 98   ASN
X_bond_distance =  3.21
sum_vdwr        =  3.15
percentages     =  1.02
angle1          =  157.84 degrees
angle2          =  83.80 degrees

 
 SALT BRIDGE FILE TEST：
 [“2X7G.pdb”，“1Jvu.pdb”，“4X7N.pdb","3ALP.pdb"]


                2X7G:
  Salt bridge #1: chain A, ARG 632 <> chain A, So4 1700
  Salt bridge #2: chain A, ARG 642 <> chain A, So4 1700
  Salt bridge #3: chain A, ARG 637 <> chain A, So4 1700
  
                1JVU:
  Salt bridge #1: chain A, THR 45 <> chain A, C2P 130
  
                4X7N:
  Salt bridge #1: chain A, ARG 903 <> chain A, TLA 1102
  
                3ALP:
  Salt bridge #1: chain A, ARG 101 <> chain A, CIT 1
  Salt bridge #2: chain B, ARG 101 <> chain B, CIT 2
   
  PI FILE TEST :
  ["1p5j.pdb","2ovp.pdb","3cuk.pdb"]

                1p5j:
  pi stacking #1: chain A,TRP 1159 <> chain A,HIS 2359
  
                20vp:
  pi stacking #1: chain A,TRY 55 <> chain A,TRY 314
  pi stacking #2: chain B,TRY 55 <> chain B,TRY 314
  pi stacking #3: chain C,TRY 55 <> chain C,TRY 314
  pi stacking #4: chain D,TRY 55 <> chain D,TRY 314
  
                3cuk:
  pi stacking #1: chain A,PHE 40 <> chain A,PLP 941
 
   
