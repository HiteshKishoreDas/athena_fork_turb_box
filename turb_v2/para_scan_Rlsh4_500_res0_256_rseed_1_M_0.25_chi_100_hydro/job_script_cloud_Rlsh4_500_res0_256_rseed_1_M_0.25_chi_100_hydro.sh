#!/bin/sh
#
#SBATCH -J h40M01
#SBATCH -o h40M01."%j".out
#SBATCH -e h40M01."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --partition=n0064
#SBATCH --mail-type=ALL
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:59:00

set -e
SECONDS=0

module purge
module load intel/19.1.2
module load impi/2019.8
module load fftw-mpi/3.3.8
module load hdf5-mpi/1.8.21
module list

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_HOME/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FFTW_HOME/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_MPI_ROOT/intel64/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_MPI_ROOT/intel64/lib/release/

cd /ptmp/hitesh/athena_fork_turb_box/turb_v2/para_scan_Rlsh4_500_res0_256_rseed_1_M_0.25_chi_100_hydro/

cp ../para_scan_Rlsh4_500_res0_256_rseed_1_M_0.25_hydro/Turb.hst ./
mv ../para_scan_Rlsh4_500_res0_256_rseed_1_M_0.25_hydro ../Turbulence/

srun ./athena_turb -r ../Turbulence/para_scan_Rlsh4_500_res0_256_rseed_1_M_0.25_hydro/Turb.final.rst -i athinput_cloud_Rlsh4_500_res0_256_rseed_1_M_0.25_chi_100_hydro.turb -t 23:55:00

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
