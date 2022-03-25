#!/bin/bash -l

#SBATCH --qos batch
#SBATCH --job-name=mm_eval                                                                       
#SBATCH --partition=bigmem                                                                        
#SBATCH --time=05:00:00                                                                       
#SBATCH --ntasks=5
# -- Update to your account number
#SBATCH --account=rcm1

# -- Update to the location of your conda environment
source /scratch2/BMC/rcm1/rhs/anaconda3/bin/activate py36_monet_default

# -- Update the location and name of your run script
cd /scratch2/BMC/rcm1/rhs/MONET/MELODIES-MONET/melodies_monet/examples/submit_jobs/
python run_melodies_monet.py
