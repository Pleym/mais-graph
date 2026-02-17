#!/usr/bin/env bash
#SBATCH --account="r250142"
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name="Graph500_perf_compare"
#SBATCH --comment="Speedup reference vs optimise"
#SBATCH --error=perf_%x_%j.err
#SBATCH --output=perf_%x_%j.out

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

# Binaries à comparer:
# - REF_BIN: implémentation de référence
# - OPT_BIN: implémentation optimisée
# Tu peux les surcharger à la soumission:
# sbatch --export=ALL,REF_BIN=./main_ref,OPT_BIN=./main_opt kernel_perf_compare.sh
REF_BIN="${REF_BIN:-./main_ref}"
OPT_BIN="${OPT_BIN:-./main_opt}"

if [ ! -x "$REF_BIN" ] || [ ! -x "$OPT_BIN" ]; then
    echo "ERREUR: binaire introuvable ou non exécutable"
    echo "REF_BIN=$REF_BIN"
    echo "OPT_BIN=$OPT_BIN"
    echo "Compile d'abord 2 exécutables (réf/optimisé) ou passe REF_BIN/OPT_BIN via --export"
    exit 1
fi

RESULTS_DIR="results_perf_${SLURM_JOBID}"
mkdir -p "$RESULTS_DIR"
CSV="$RESULTS_DIR/speedup.csv"

echo "threads,k2_ref_ms,k2_opt_ms,k2_speedup,k3_ref_ms,k3_opt_ms,k3_speedup" > "$CSV"

extract_k2_ms() {
    awk -F'avg_time=' '/K2 \(bfs .*\): avg_time=/{split($2,a," "); print a[1]; exit}' "$1"
}

extract_k3_ms() {
    awk -F'avg_time=' '/K3 \(sssp .*\): avg_time=/{split($2,a," "); print a[1]; exit}' "$1"
}

SCALE=12
EDGE_FACTOR=16
ROOTS=16
THREADS_LIST="1 2 4 8 16 32 64"

for t in $THREADS_LIST; do
    echo "=== Perf compare: threads=$t ===" | tee -a "$RESULTS_DIR/run_status.log"

    REF_LOG="$RESULTS_DIR/ref_t${t}.log"
    OPT_LOG="$RESULTS_DIR/opt_t${t}.log"

    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
    srun -n1 -c "$t" --cpu-bind=cores "$REF_BIN" "$SCALE" "$EDGE_FACTOR" "$ROOTS" > "$REF_LOG" 2>&1

    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
    srun -n1 -c "$t" --cpu-bind=cores "$OPT_BIN" "$SCALE" "$EDGE_FACTOR" "$ROOTS" > "$OPT_LOG" 2>&1

    k2_ref=$(extract_k2_ms "$REF_LOG")
    k2_opt=$(extract_k2_ms "$OPT_LOG")
    k3_ref=$(extract_k3_ms "$REF_LOG")
    k3_opt=$(extract_k3_ms "$OPT_LOG")

    k2_speedup=$(awk -v a="$k2_ref" -v b="$k2_opt" 'BEGIN{if (b>0) printf "%.6f", a/b; else print "nan"}')
    k3_speedup=$(awk -v a="$k3_ref" -v b="$k3_opt" 'BEGIN{if (b>0) printf "%.6f", a/b; else print "nan"}')

    echo "$t,$k2_ref,$k2_opt,$k2_speedup,$k3_ref,$k3_opt,$k3_speedup" >> "$CSV"
done

echo "Comparaison terminée. Résultats: $CSV"
cat "$CSV"
