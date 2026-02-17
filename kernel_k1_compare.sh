#!/usr/bin/env bash
#SBATCH --account="r250142"
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name="Graph500_k1_compare"
#SBATCH --comment="Compare K1 builders vs K2/K3"
#SBATCH --error=k1compare_%x_%j.err
#SBATCH --output=k1compare_%x_%j.out

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

# Build k1_compare (not in Makefile)
mpicxx -Drestrict=__restrict__ -DREUSE_CSR_FOR_VALIDATION -march=znver4 -Ofast -flto -funroll-loops -fomit-frame-pointer -fno-math-errno -DNDEBUG -fopenmp \
    $CPPFLAGS -o k1_compare k1_compare.cpp generator/*.cpp kernels/*.cpp $LDFLAGS $LDLIBS

RESULTS_DIR="results_k1compare_${SLURM_JOBID}"
mkdir -p "$RESULTS_DIR"
CSV="$RESULTS_DIR/k1_compare.csv"

echo "threads,variant,k1_ms,k2_ms,k3_ms,k3_teps" > "$CSV"

SCALE="${SCALE:-12}"
EDGE_FACTOR="${EDGE_FACTOR:-16}"
ROOTS="${ROOTS:-4}"
THREADS_LIST="${THREADS_LIST:-1 2 4 8 16 32 64}"
MODE="${MODE:-all_opt}"

append_variant() {
    local variant="$1" log="$2" threads="$3"
    awk -v var="$variant" -v t="$threads" '
        $0 ~ var {
            k1="nan"; k2="nan"; k3="nan"; teps="nan";
            if (match($0, /K1_ms=([0-9.e+-]+)/, a)) k1=a[1];
            if (match($0, /K2_ms=([0-9.e+-]+)/, a)) k2=a[1];
            if (match($0, /K3_ms=([0-9.e+-]+)/, a)) k3=a[1];
            if (match($0, /K3_teps=([0-9.e+-]+)/, a)) teps=a[1];
            printf("%s,%s,%s,%s,%s,%s\n", t, var, k1, k2, k3, teps);
        }
    ' "$log" >> "$CSV"
}

for t in $THREADS_LIST; do
    LOG="$RESULTS_DIR/k1_t${t}.log"
    echo "=== K1 compare: threads=$t scale=$SCALE edge_factor=$EDGE_FACTOR roots=$ROOTS mode=$MODE ===" | tee -a "$RESULTS_DIR/run_status.log"

    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
    srun -n1 -c "$t" --cpu-bind=cores ./k1_compare "$SCALE" "$EDGE_FACTOR" "$ROOTS" "$MODE" > "$LOG" 2>&1

    append_variant "v2_sort" "$LOG" "$t"
    append_variant "v3_parallel" "$LOG" "$t"

done

echo "Terminé. Logs: $RESULTS_DIR" \
     "CSV: $CSV"
cat "$CSV"
