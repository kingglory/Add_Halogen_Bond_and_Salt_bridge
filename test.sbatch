#!/bin/bash -l
#SBATCH -o run.%j.%N.out
#SBATCH -p cpu
#SBATCH -N 1
##SBATCH --nodelist=cu[01-12]
#SBATCH --cpus-per-task=4         
#SBATCH -J test_song 

/home/wensong/phenix-1.14rc1-3161/build/bin/python run.py > test_o.log
