#!/bin/bash -l

#SBATCH --qos batch
#SBATCH --job-name=mm_eval
#SBATCH --partition=bigmem
#SBATCH --time=05:00:00
#SBATCH --ntasks=5
# -- Update to your account number
#SBATCH --account=rcm2

# -- Update to the location of your conda environment
source /scratch4/BMC/rcm1/rhs/miniconda3/bin/activate melodies-monet-test-examples

# -- Run MELODIES MONET through a command line call instead
cd /scratch3/BMC/rcm2/rhs/MM_AEROMMA/test_examples/code2/MELODIES-MONET/examples/yaml/
melodies-monet run control_ufsaqm_airnow_surface-submit.yaml
