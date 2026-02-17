#include "generator.hpp"

#include "make_graph.hpp"
#include <stddef.h>

void print_edge_list(edge_list* list) {
	printf("Number of edges: %zu\n", list->length);
	for (size_t i = 0; i < list->length; i++) {
		printf("(%zu, %zu): %lf\n", list->edges[i].v0, list->edges[i].v1, list->weights[i]);
	}
}

edge_list generate_graph(uint8_t scale, uint8_t edge_factor) {
	size_t nb_edges = (size_t)1 << (scale + edge_factor - 1);
	edge_list list = {
		.length = nb_edges,
		.edges = (packed_edge*)malloc(sizeof(packed_edge) * nb_edges),
		.weights = (float*)malloc(sizeof(float) * nb_edges),
	};

	// Generate the graph edges
	int64_t got_nb_edges;
	int64_t desired_nb_edges = (int64_t)1 << (scale + edge_factor - 1);
	make_graph(scale, desired_nb_edges, 0, 0, &got_nb_edges, &list.edges);

	list.length = (size_t)got_nb_edges;

	// Generate random weights (floats) between 0 and 1 for each edge
	double* weights = (double*)malloc(sizeof(double) * list.length);
	make_random_numbers(got_nb_edges, 0, 0, 0, weights);
	for (size_t i = 0; i < list.length; i++) {
		list.weights[i] = (float)weights[i];
	}
	free(weights);

	return list;
}

edge_list generate_graph_toy() {
	return generate_graph(26, 16);
}

edge_list generate_graph_mini() {
	return generate_graph(29, 16);
}

edge_list generate_graph_small() {
	return generate_graph(32, 16);
}

void edge_list_destroy(edge_list list) {
	free(list.edges);
	free(list.weights);
}