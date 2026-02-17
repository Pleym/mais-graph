# mais-graph

Version minimaliste centrée sur les implémentations de référence pour les 3 kernels Graph500:

- K1: génération + construction du graphe CSR (`from_edge_list`)
- K2: BFS de référence (`bfs_formal`)
- K3: SSSP de référence (`sssp_bf`)

## Build

```bash
make
```

## Run

```bash
./main [scale] [edge_factor] [roots] [mode]
```

Exemple:

```bash
./main 12 16 16 ref
```

Où:

- `scale` = log2(nombre de sommets)
- `edge_factor` = facteur d'arêtes
- `roots` = nombre de sources testées pour K2/K3
- `mode` = `ref` | `bfs_opt` | `sssp_opt` | `all_opt` | `hybrid`

## Scripts SLURM ROMEO

- [kernel_1_2_experiments.sh](kernel_1_2_experiments.sh): strong scaling (taille fixe, threads variables)
- [kernel_weak_scaling.sh](kernel_weak_scaling.sh): weak scaling (taille augmente avec le nombre de threads)
- [kernel_perf_compare.sh](kernel_perf_compare.sh): comparaison référence vs optimisé avec speedup K2/K3
- [kernel_optimization_pipeline.sh](kernel_optimization_pipeline.sh): benchmark optimisation par optimisation (CSV speedup vs `ref`)

Soumission:

```bash
sbatch kernel_1_2_experiments.sh
sbatch kernel_weak_scaling.sh
sbatch --export=ALL,REF_BIN=./main_ref,OPT_BIN=./main_opt kernel_perf_compare.sh
sbatch kernel_optimization_pipeline.sh
```

Nommage des logs SLURM:

- strong: `strong_Graph500_strong_scaling_<jobid>.out/.err`
- weak: `weak_Graph500_weak_scaling_<jobid>.out/.err`
- perf: `perf_Graph500_perf_compare_<jobid>.out/.err`

Chaque script écrit aussi un `run_status.log` dans son dossier `results_*_<jobid>/`.