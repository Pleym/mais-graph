#!/usr/bin/env bash
#SBATCH --account="r250142"
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name="Graph500_strong_scaling"
#SBATCH --comment="Strong scaling K1/K2/K3"
#SBATCH --error=strong_%x_%j.err
#SBATCH --output=strong_%x_%j.out

set -euo pipefail

romeo_load_x64cpu_env

spack load mpich tbb intel-tbb
TBB_DIR=$(spack location -i tbb 2>/dev/null || spack location -i intel-tbb 2>/dev/null)
TBB_LIB_DIR="$TBB_DIR/lib64"
ls -la "$TBB_LIB_DIR"

export CPPFLAGS="-I$TBB_DIR/include"
export CXXFLAGS="-I$TBB_DIR/include"
export LDFLAGS="-L$TBB_LIB_DIR -Wl,-rpath,$TBB_LIB_DIR"
export LDLIBS="-ltbb -ltbbmalloc -ltbbmalloc_proxy"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

make -B -j"$(nproc)" CPPFLAGS="$CPPFLAGS" CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS"

RESULTS_DIR="results_strong_${SLURM_JOBID}"
mkdir -p "$RESULTS_DIR"

# Strong scaling: problème fixe, on varie uniquement le nombre de threads
SCALE=12
EDGE_FACTOR=16
ROOTS=16

THREADS_LIST="1 2 4 8 16 32 64"
for t in $THREADS_LIST; do
        echo "=== Strong scaling: threads=$t scale=$SCALE edge_factor=$EDGE_FACTOR roots=$ROOTS ===" | tee -a "$RESULTS_DIR/run_status.log"
        OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
        srun -n1 -c "$t" --cpu-bind=cores ./main "$SCALE" "$EDGE_FACTOR" "$ROOTS" \
                > "$RESULTS_DIR/strong_t${t}.log" 2>&1
done

echo "Strong scaling terminé. Logs: $RESULTS_DIR"