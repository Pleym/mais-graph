/**
 * @file generator/generator.h
 * @brief Header file for edge list generation and manipulation functions.
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include "graph_generator.hpp"

/**
 * @brief The generated edge list.
 */
typedef struct edge_list {
	size_t length;		///< Number of edges in the list
	packed_edge* edges; ///< Array of edges (v0, v1)
	float* weights;		///< Array of edge weights
} edge_list;

/**
 * @brief Print the edge list in stdout
 */
void print_edge_list(edge_list* list);

/**
 * @brief Generate a graph with a given scale (number of nodes = 2^scale) and edge factor (number of edges per node on average)
 */
edge_list generate_graph(uint8_t scale, uint8_t edge_factor);

/**
 * @brief Generate the toy graph from the graph500 specification (scale = 26, edge_factor = 16)
 */
edge_list generate_graph_toy();

/**
 * @brief Generate the mini graph from the graph500 specification (scale = 29, edge_factor = 16)
 */
edge_list generate_graph_mini();

/**
 * @brief Generate the small graph from the graph500 specification (scale = 32, edge_factor = 16)
 */
edge_list generate_graph_small();

/**
 * @brief Free the memory allocated for the edge list
 */
void edge_list_destroy(edge_list list);

#endif // GENERATOR_H