#!/usr/bin/env bash
#SBATCH --account="r250142"
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name="Graph500_optim_pipeline"
#SBATCH --comment="Step-by-step optimization benchmark"
#SBATCH --error=pipeline_%x_%j.err
#SBATCH --output=pipeline_%x_%j.out

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

make -B -j"$(nproc)" CPPFLAGS="$CPPFLAGS" CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS"

RESULTS_DIR="results_pipeline_${SLURM_JOBID}"
mkdir -p "$RESULTS_DIR"
CSV="$RESULTS_DIR/optimization_pipeline.csv"

echo "threads,mode,k2_ms,k3_ms,k2_speedup_vs_ref,k3_speedup_vs_ref" > "$CSV"

extract_k2_ms() {
    awk -F'avg_time=' '/K2 \(bfs .*\): avg_time=/{split($2,a," "); print a[1]; exit}' "$1"
}

extract_k3_ms() {
    awk -F'avg_time=' '/K3 \(sssp .*\): avg_time=/{split($2,a," "); print a[1]; exit}' "$1"
}

SCALE="${SCALE:-12}"
EDGE_FACTOR="${EDGE_FACTOR:-16}"
ROOTS="${ROOTS:-16}"
THREADS_LIST="${THREADS_LIST:-1 2 4 8 16 32 64}"
MODES="bfs_opt sssp_opt all_opt hybrid"

for t in $THREADS_LIST; do
    echo "=== Pipeline: threads=$t (reference) ===" | tee -a "$RESULTS_DIR/run_status.log"
    REF_LOG="$RESULTS_DIR/ref_t${t}.log"

    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
    srun -n1 -c "$t" --cpu-bind=cores ./main "$SCALE" "$EDGE_FACTOR" "$ROOTS" ref > "$REF_LOG" 2>&1

    k2_ref=$(extract_k2_ms "$REF_LOG")
    k3_ref=$(extract_k3_ms "$REF_LOG")
    echo "$t,ref,$k2_ref,$k3_ref,1.000000,1.000000" >> "$CSV"

    for mode in $MODES; do
        echo "=== Pipeline: threads=$t mode=$mode ===" | tee -a "$RESULTS_DIR/run_status.log"
        MODE_LOG="$RESULTS_DIR/${mode}_t${t}.log"

        OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
        srun -n1 -c "$t" --cpu-bind=cores ./main "$SCALE" "$EDGE_FACTOR" "$ROOTS" "$mode" > "$MODE_LOG" 2>&1

        k2_mode=$(extract_k2_ms "$MODE_LOG")
        k3_mode=$(extract_k3_ms "$MODE_LOG")

        k2_speedup=$(awk -v a="$k2_ref" -v b="$k2_mode" 'BEGIN{if (b>0) printf "%.6f", a/b; else print "nan"}')
        k3_speedup=$(awk -v a="$k3_ref" -v b="$k3_mode" 'BEGIN{if (b>0) printf "%.6f", a/b; else print "nan"}')

        echo "$t,$mode,$k2_mode,$k3_mode,$k2_speedup,$k3_speedup" >> "$CSV"
    done
done

echo "Pipeline terminé. CSV: $CSV"
cat "$CSV"
