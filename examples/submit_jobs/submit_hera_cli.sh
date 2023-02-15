#!/bin/bash -l

#SBATCH --qos batch
#SBATCH --job-name=mm_eval
#SBATCH --partition=bigmem
#SBATCH --time=05:00:00
#SBATCH --ntasks=5
# -- Update to your account number
#SBATCH --account=rcm2

# -- Update to the location of your conda environment
source /scratch2/BMC/rcm1/rhs/miniconda3/bin/activate mm_dev_cli

# -- Run MELODIES MONET through a command line call instead
cd /scratch2/BMC/rcm1/rhs/MONET/main/wrf-chem/qindan_test/
melodies-monet run control_wrfchem_cli.yaml
