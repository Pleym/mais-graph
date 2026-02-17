#!/usr/bin/env bash
#SBATCH --account="r250142"
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --job-name="Graph500_weak_scaling"
#SBATCH --comment="Weak scaling K1/K2/K3"
#SBATCH --error=weak_%x_%j.err
#SBATCH --output=weak_%x_%j.out

set -euo pipefail

romeo_load_x64cpu_env

spack load mpich tbb intel-tbb
TBB_DIR=$(spack location -i tbb 2>/dev/null || spack location -i intel-tbb 2>/dev/null)
TBB_LIB_DIR="$TBB_DIR/lib64"

export CPPFLAGS="-I$TBB_DIR/include"
export CXXFLAGS="-I$TBB_DIR/include"
export LDFLAGS="-L$TBB_LIB_DIR -Wl,-rpath,$TBB_LIB_DIR"
export LDLIBS="-ltbb -ltbbmalloc -ltbbmalloc_proxy"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

make -j"$(nproc)" CPPFLAGS="$CPPFLAGS" CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS"

RESULTS_DIR="results_weak_${SLURM_JOBID}"
mkdir -p "$RESULTS_DIR"

# Weak scaling: charge par thread approximativement constante
# Ici, SCALE augmente avec log2(threads)
BASE_SCALE=10
EDGE_FACTOR=16
ROOTS=16

# Inclut profils monodomaine (<=96 threads socket) et bi-socket (192 threads)
THREADS_LIST="1 2 4 8 16 32 64 96 128 192"
for t in $THREADS_LIST; do
    exp=0
    tmp=$t
    while [ "$tmp" -gt 1 ]; do
        tmp=$((tmp / 2))
        exp=$((exp + 1))
    done

    SCALE=$((BASE_SCALE + exp))

    echo "=== Weak scaling: threads=$t scale=$SCALE edge_factor=$EDGE_FACTOR roots=$ROOTS ===" | tee -a "$RESULTS_DIR/run_status.log"

    # Affinité: si t<=96 on reste sur un socket; sinon interleave mem + spread cores
    if [ "$t" -le 96 ]; then
        BIND_OPTS=(--cpunodebind=0 --membind=0)
        OMP_BIND=spread
    else
        BIND_OPTS=(--interleave=all)
        OMP_BIND=spread
    fi

    OMP_NUM_THREADS=$t OMP_PROC_BIND=$OMP_BIND OMP_PLACES=cores \
    srun -n1 -c "$t" "${BIND_OPTS[@]}" ./main "$SCALE" "$EDGE_FACTOR" "$ROOTS" \
        > "$RESULTS_DIR/weak_t${t}_scale${SCALE}.log" 2>&1
done

echo "Weak scaling terminé. Logs: $RESULTS_DIR"
