#ifndef BREADTH_SEARCH_H
#define BREADTH_SEARCH_H

#include "gen_graph.hpp"
#include <cstdint>

typedef struct bfs_result {
	int64_t* parent_array;
	double teps;
	double time_ms;

	~bfs_result() {
		free(parent_array);
	}
} bfs_result;

bfs_result bfs(graph& g, int64_t source);
bfs_result bfs_formal(graph& g, int64_t source);
bfs_result bfs_top_down_omp_cas(graph& g, int64_t source);
bfs_result bfs_hybrid(graph& g, int64_t source);

#endif