#!/bin/sh
#
#SBATCH -J athena_turb_wcool
#SBATCH -o athena_turb_wcool."%j".out
#SBATCH -e athena_turb_wcool."%j".err
#SBATCH --partition=p.test
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:29:00

module purge
module load anaconda 
module load ffmpeg
module list

cd /u/hitesh/athena_plusplus_townsend/plot_scripts/

# ./turb_compile_freya
# cp ../bin/athena ./

srun python load_hdf5.py
# ./turb_run athinput.turb 64 
