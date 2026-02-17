#!/usr/bin/env bash
#SBATCH --account="r250142"
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name "Graph500_k1_k2"
#SBATCH --comment "Benchmark nos algos de graph500"
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

romeo_load_x64cpu_env

spack load mpich tbb intel-tbb
TBB_DIR=$(spack location -i tbb 2>/dev/null || spack location -i intel-tbb 2>/dev/null)
TBB_LIB_DIR="$TBB_DIR/lib64"
ls -la "$TBB_LIB_DIR"

export CPPFLAGS="-I$TBB_DIR/include"
export CXXFLAGS="-I$TBB_DIR/include"
export LDFLAGS="-L$TBB_LIB_DIR -Wl,-rpath,$TBB_LIB_DIR"
export LDLIBS="-ltbb -ltbbmalloc -ltbbmalloc_proxy"

make -B -j$(nproc) CPPFLAGS="$CPPFLAGS" CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS"

# Run `./main` with different OMP thread counts and save per-run logs
THREADS_LIST="1 2 4 8 16 32 64"
for t in $THREADS_LIST; do
        echo "=== Running with $t threads ===" | tee -a job.$SLURM_JOBID.out
        OMP_NUM_THREADS=$t srun -n1 -c $t ./main
done