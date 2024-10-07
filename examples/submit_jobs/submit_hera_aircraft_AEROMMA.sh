#!/bin/bash -l

#SBATCH --qos batch
#SBATCH --job-name=mm_eval                                                                       
#SBATCH --partition=bigmem                                                                        
#SBATCH --time=05:00:00                                                                       
#SBATCH --ntasks=5
# -- Update to your account number
#SBATCH --account=rcm1

# -- Update to the location of your conda environment
source /scratch2/BMC/rcm1/rhs/miniconda3/bin/activate melodies-monet

# -- Update the location and name of your run script
cd /scratch1/BMC/rcm2/rhs/monet_example/AEROMMA/submit_jobs
python run_aircraft_pairing_loop_AEROMMA.py
