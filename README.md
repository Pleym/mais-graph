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
./main [scale] [edge_factor] [roots]
```

Exemple:

```bash
./main 12 16 16
```

Où:

- `scale` = log2(nombre de sommets)
- `edge_factor` = facteur d'arêtes
- `roots` = nombre de sources testées pour K2/K3