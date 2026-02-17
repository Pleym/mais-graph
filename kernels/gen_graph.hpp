#ifndef GEN_GRAPH_H
#define GEN_GRAPH_H

#include "../generator/generator.hpp"
#include <cstdint>

typedef struct {
	int64_t* slicing_idx;
	int64_t* neighbors;
	float* weights;
	double time_ms;
	int64_t length;
	int64_t nb_nodes;
} graph;

graph from_edge_list(edge_list input_list);

void graph_destroy(graph& g);

// Utilities
template <typename F> static inline void for_each_neighbor(const graph& g, int64_t node, F f) {
	int64_t start = g.slicing_idx[node];
	int64_t end = g.slicing_idx[node + 1];
	for (int64_t i = start; i < end; i++) {
		f(g.neighbors[i], g.weights[i]);
	}
}

template <typename F> static inline void parallel_for_each_neighbor(const graph& g, int64_t node, F f) {
	int64_t start = g.slicing_idx[node];
	int64_t end = g.slicing_idx[node + 1];
#pragma omp parallel for
	for (int64_t i = start; i < end; i++) {
		f(g.neighbors[i], g.weights[i]);
	}
}

static uint64_t degree_of_node(const graph& g, int64_t node) {
	return g.slicing_idx[node + 1] - g.slicing_idx[node];
}

// Compare two graphs for structural equality. Returns true if identical.
bool compare_graphs(const graph& a, const graph& b, double weight_eps = 1e-6);

#endif