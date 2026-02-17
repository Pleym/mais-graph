#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

#include "gen_graph.hpp"
#include <cstdint>

typedef struct {
	int64_t* parent_array;
	float* distance_array;
	double teps;
	double time_ms;

} shortest_path;

shortest_path sssp_bf(graph source_graph, int64_t root);
shortest_path sssp_bf_frontier_omp(graph source_graph, int64_t root);
void shortest_path_destroy(shortest_path& sp);

#endif