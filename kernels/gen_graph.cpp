#include "gen_graph.hpp"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <vector>

void graph_destroy(graph& g) {
	free(g.neighbors);
	free(g.weights);
	free(g.slicing_idx);
}

void print_slicing_idx(const graph& g) {
	std::cout << "slicing_idx: ";
	for (int64_t i = 0; i <= g.nb_nodes; ++i) {
		std::cout << g.slicing_idx[i] << " ";
	}
	std::cout << std::endl;
}

void print_neighbors(const graph& g) {
	for (int64_t node = 0; node < g.nb_nodes; node++) {
		std::cout << "Node " << node << ": ";
		for (int64_t i = g.slicing_idx[node]; i < g.slicing_idx[node + 1]; ++i) {
			std::cout << "(" << g.neighbors[i] << ", " << g.weights[i] << ") ";
		}
		std::cout << std::endl;
	}
}

void print_neighbors_array(const graph& g) {
	std::cout << "neighbors array: ";
	for (int64_t i = 0; i < g.length; ++i) {
		std::cout << g.neighbors[i] << " ";
	}
	std::cout << std::endl;
}

static void print_neighbors_of_node(const graph& graph_obj, int64_t node) {
	std::vector<int64_t> neighbors;
	for (int64_t i = graph_obj.slicing_idx[node]; i < graph_obj.slicing_idx[node + 1]; ++i) {
		neighbors.push_back(graph_obj.neighbors[i]);
	}
	std::sort(neighbors.begin(), neighbors.end());
	std::cout << "Neighbors of node " << node << ": ";
	for (auto n : neighbors) {
		std::cout << n << " ";
	}
	std::cout << '\n';
}

graph from_edge_list_v1(edge_list input_list) {
	auto start = std::chrono::high_resolution_clock::now();
	int64_t max_node = 0;
	for (int64_t i = 0; i < input_list.length; i++) {
		max_node = std::max(input_list.edges[i].v0, max_node);
		max_node = std::max(input_list.edges[i].v1, max_node);
	}
	graph g;
	g.neighbors = (int64_t*)malloc(2 * input_list.length * sizeof(int64_t));
	g.weights = (float*)malloc(2 * input_list.length * sizeof(float));
	g.slicing_idx = (int64_t*)malloc((max_node + 2) * sizeof(int64_t));
	int64_t cpt = 0;
	int64_t cpt_idx = 1;
	g.slicing_idx[0] = cpt;
	for (int64_t node = 0; node <= max_node; node++) {
		for (int64_t i = 0; i < input_list.length; i++) {
			if (input_list.edges[i].v0 != input_list.edges[i].v1) {

				int64_t neighbor = -1;

				if (input_list.edges[i].v0 == node) {
					neighbor = input_list.edges[i].v1;
				}
				else if (input_list.edges[i].v1 == node) {
					neighbor = input_list.edges[i].v0;
				}

				if (neighbor != -1) {
					bool found = false;
					for (int64_t j = g.slicing_idx[cpt_idx - 1]; j < cpt; j++) {
						if (g.neighbors[j] == neighbor) {
							g.weights[j] = input_list.weights[i];
							found = true;
							break;
						}
					}

					if (!found) {
						g.neighbors[cpt] = neighbor;
						g.weights[cpt] = input_list.weights[i];
						cpt++;
					}
				}
			}
		}
		g.slicing_idx[cpt_idx] = cpt;
		cpt_idx++;
	}
	g.length = cpt;
	g.nb_nodes = max_node + 1;
	g.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
	return g;
}

#include <cmath>

typedef struct {
	int64_t u, v;
	float weight;
} edge;

bool are_same_edge(const edge& e1, const edge& e2) {
	return e1.u == e2.u && e1.v == e2.v;
}

bool is_loop(const edge& e) {
	return e.u == e.v;
}

bool is_loop(const packed_edge& e) {
	return e.v0 == e.v1;
}

graph from_edge_list_v2(edge_list input_list) {
	auto start = std::chrono::high_resolution_clock::now();
	graph g;

	edge* edges = (edge*)malloc(sizeof(edge) * 2 * input_list.length);
	auto original_edges_ptr = edges;

	// Duplicate the edges
	auto start_duplication = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < input_list.length; i++) {
		edges[2 * i] = (edge){
			.u = input_list.edges[i].v1,
			.v = input_list.edges[i].v0,
			.weight = input_list.weights[i],
		};
		edges[(2 * i) + 1] = (edge){
			.u = input_list.edges[i].v0,
			.v = input_list.edges[i].v1,
			.weight = input_list.weights[i],
		};
	}
	printf("Took %f seconds to duplicate edges\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_duplication).count());

	// Sort them lexicographically
	auto start_sort = std::chrono::high_resolution_clock::now();
	std::sort(edges, edges + (input_list.length * 2), [&](const edge& e1, const edge& e2) {
		if (e1.u != e2.u) {
			return e1.u < e2.u;
		}
		return e1.v < e2.v;
	});
	printf("Took %f seconds to sort edges\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_sort).count());

	auto start_build = std::chrono::high_resolution_clock::now();
	int64_t nb_nodes = edges[(input_list.length * 2) - 1].u + 1;
	g.slicing_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));
	g.neighbors = (int64_t*)malloc(2 * input_list.length * sizeof(int64_t));
	g.weights = (float*)malloc(2 * input_list.length * sizeof(float));

	// If there are loops at the beginning
	int nb_skipped = 0;
	while (is_loop(edges[0])) {
		edges++;
		nb_skipped++;
	}

	int64_t nb_nodes_so_far = 0;
	auto first_node = edges[0].u;
	g.slicing_idx[first_node] = 0;
	g.neighbors[0] = edges[0].v;
	g.weights[0] = edges[0].weight;
	int64_t nb_neighbors_so_far = 1;
	for (int64_t i = 1; i < (2 * input_list.length) - nb_skipped; i++) {
		if (are_same_edge(edges[i], edges[i - 1])) {
			continue;
		}

		if (is_loop(edges[i])) {
			continue;
		}

		if (edges[i].u != edges[i - 1].u) {
			// Fill slicing_idx for intermediary nodes with no neighbors
			for (int64_t missing = edges[i - 1].u + 1; missing < edges[i].u; ++missing) {
				nb_nodes_so_far++;
				g.slicing_idx[nb_nodes_so_far] = nb_neighbors_so_far;
			}
			// Different node -> mark a slice
			nb_nodes_so_far++;
			g.slicing_idx[nb_nodes_so_far] = nb_neighbors_so_far;
		}

		// Push v as a neighbor of u
		g.neighbors[nb_neighbors_so_far] = edges[i].v;
		g.weights[nb_neighbors_so_far] = edges[i].weight;
		nb_neighbors_so_far++;
	}
	// Close off the slices
	nb_nodes_so_far++;
	g.slicing_idx[nb_nodes_so_far] = nb_neighbors_so_far;
	g.length = nb_neighbors_so_far;
	g.nb_nodes = nb_nodes;
	printf("Took %f seconds to build CSR\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_build).count());

	g.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();

	// Verification avec vieux algorithme
	// graph g2 = from_edge_list_v1(input_list);
	// print_slicing_idx(g);
	// print_slicing_idx(g2);
	// for (int64_t node = 0; node < g.nb_nodes; ++node) {
	// 	print_neighbors_of_node(g, node);
	// 	print_neighbors_of_node(g2, node);
	// }

	free(original_edges_ptr);

	return g;
}

// Fast two-pass CSR builder (counts -> prefix-sum -> fill).
// - emits undirected edges (u->v and v->u)
// - O(n + m) time, no global sort
graph from_edge_list_v3(edge_list input_list) {
	auto start = std::chrono::high_resolution_clock::now();
	graph g;

	// Compute number of nodes
	auto start_count = std::chrono::high_resolution_clock::now();
	int64_t max_node = -1;
	for (int64_t i = 0; i < input_list.length; ++i) {
		max_node = std::max(max_node, input_list.edges[i].v0);
		max_node = std::max(max_node, input_list.edges[i].v1);
	}
	int64_t nb_nodes = max_node + 1;
	printf("Took %f seconds to compute number of nodes\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_count).count());

	// Count degrees
	auto start_degree = std::chrono::high_resolution_clock::now();
	std::vector<int64_t> degrees(nb_nodes);
	for (int64_t i = 0; i < input_list.length; ++i) {
		// Skip loops
		if (is_loop(input_list.edges[i])) {
			continue;
		}
		degrees[input_list.edges[i].v0]++;
		degrees[input_list.edges[i].v1]++;
	}
	printf("Took %f seconds to count degrees\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_degree).count());

	// Do prefix sum of degrees to get slicing_idx
	auto start_prefix_sum = std::chrono::high_resolution_clock::now();
	g.slicing_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));
	int64_t sum = 0;
	for (int64_t i = 0; i < nb_nodes; ++i) {
		g.slicing_idx[i] = sum;
		sum += degrees[i];
	}
	g.slicing_idx[nb_nodes] = sum;
	int64_t nb_edges = sum;
	printf("Took %f seconds to compute prefix sum\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_prefix_sum).count());

	// Allocate neighbor arrays
	g.neighbors = (int64_t*)malloc(nb_edges * sizeof(int64_t));
	g.weights = (float*)malloc(nb_edges * sizeof(float));

	// Fill using per-node cursors
	auto start_fill = std::chrono::high_resolution_clock::now();
	std::vector<int64_t> next_write_idx(nb_nodes);
	for (int64_t i = 0; i < nb_nodes; ++i) {
		next_write_idx[i] = g.slicing_idx[i];
	}
	for (int64_t i = 0; i < input_list.length; ++i) {
		int64_t u = input_list.edges[i].v0;
		int64_t v = input_list.edges[i].v1;

		// Skip loops
		if (is_loop(input_list.edges[i])) {
			continue;
		}

		float weight = input_list.weights[i];

		// Write v as a neighbor of u
		int64_t write_idx_u = next_write_idx[u];
		g.neighbors[write_idx_u] = v;
		g.weights[write_idx_u] = weight;

		// Write u as a neighbor of v
		int64_t write_idx_v = next_write_idx[v];
		g.neighbors[write_idx_v] = u;
		g.weights[write_idx_v] = weight;

		next_write_idx[u]++;
		next_write_idx[v]++;
	}
	printf("Took %f seconds to fill neighbor arrays\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_fill).count());

	g.length = nb_edges;
	g.nb_nodes = nb_nodes;

	// On-the-fly per-node deduplication using a marker table
	auto start_dedup = std::chrono::high_resolution_clock::now();
	int64_t old_nb_edges = nb_edges;
	int64_t* new_neighbors = (int64_t*)malloc(old_nb_edges * sizeof(int64_t));
	float* new_weights = (float*)malloc(old_nb_edges * sizeof(float));
	int64_t* new_slices_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));

	std::vector<int64_t> marker(nb_nodes, -1); // marker[v] = position in new_neighbors or -1
	std::vector<int64_t> touched;			   // list of nodes whose marker was touched, to reset later
	touched.reserve(64);

	int64_t write_pos = 0;
	for (int64_t node = 0; node < nb_nodes; node++) {
		new_slices_idx[node] = write_pos;
		for_each_neighbor(g, node, [&](int64_t neighbor, float weight) {
			if (marker[neighbor] == -1) {
				marker[neighbor] = write_pos;
				new_neighbors[write_pos] = neighbor;
				new_weights[write_pos] = weight;
				touched.push_back(neighbor);
				write_pos++;
			}
			else {
				// Duplicate neighbor, just keep first weight encountered
			}
		});
		// Reset markers for touched neighbors of this node
		for (int64_t node : touched) {
			marker[node] = -1;
		}
		touched.clear();
	}

	new_slices_idx[nb_nodes] = write_pos;

	// replace arrays
	free(g.neighbors);
	free(g.weights);
	free(g.slicing_idx);
	g.neighbors = new_neighbors;
	g.weights = new_weights;
	g.slicing_idx = new_slices_idx;
	g.length = write_pos;

	printf("Took %f seconds to deduplicate neighbors\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_dedup).count());

	g.time_ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();

	// Verification avec vieux algorithme
	// graph g_ref = from_edge_list_v2(input_list);
	// print_slicing_idx(g);
	// print_slicing_idx(g_ref);
	// for (int64_t node = 0; node < g.nb_nodes; ++node) {
	// 	print_neighbors_of_node(g, node);
	// 	print_neighbors_of_node(g_ref, node);
	// }

	return g;
}

graph from_edge_list_v2_parallel(edge_list input_list) {
	auto start = std::chrono::high_resolution_clock::now();
	graph g;

	edge* edges = (edge*)malloc(sizeof(edge) * 2 * input_list.length);
	auto original_edges_ptr = edges;

	auto start_duplication = std::chrono::high_resolution_clock::now();
	// #pragma omp parallel for schedule(static) // <- Turned out to be way worse in parallel, maybe due to false sharing.
	// Duplicate the edges
	for (size_t i = 0; i < input_list.length; i++) {
		edges[2 * i] = (edge){
			.u = input_list.edges[i].v1,
			.v = input_list.edges[i].v0,
			.weight = input_list.weights[i],
		};
		edges[(2 * i) + 1] = (edge){
			.u = input_list.edges[i].v0,
			.v = input_list.edges[i].v1,
			.weight = input_list.weights[i],
		};
	}
	printf("Took %f seconds to duplicate edges\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_duplication).count());

	// Sort them lexicographically
	auto start_sort = std::chrono::high_resolution_clock::now();
							std::sort(edges, edges + (input_list.length * 2), [&](const edge& e1, const edge& e2) {
		if (e1.u != e2.u) {
			return e1.u < e2.u;
		}
		return e1.v < e2.v;
	});
	printf("Took %f seconds to sort edges\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_sort).count());

	auto start_build = std::chrono::high_resolution_clock::now();
	int64_t nb_nodes = edges[(input_list.length * 2) - 1].u + 1;
	g.slicing_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));
	g.neighbors = (int64_t*)malloc(2 * input_list.length * sizeof(int64_t));
	g.weights = (float*)malloc(2 * input_list.length * sizeof(float));

	// If there are loops at the beginning
	int nb_skipped = 0;
	while (is_loop(edges[0])) {
		edges++;
		nb_skipped++;
	}

	int64_t nb_nodes_so_far = 0;
	auto first_node = edges[0].u;
	g.slicing_idx[first_node] = 0;
	g.neighbors[0] = edges[0].v;
	g.weights[0] = edges[0].weight;
	int64_t nb_neighbors_so_far = 1;
	for (int64_t i = 1; i < (2 * input_list.length) - nb_skipped; i++) {
		if (are_same_edge(edges[i], edges[i - 1])) {
			continue;
		}

		if (is_loop(edges[i])) {
			continue;
		}

		if (edges[i].u != edges[i - 1].u) {
			// Fill slicing_idx for intermediary nodes with no neighbors
			for (int64_t missing = edges[i - 1].u + 1; missing < edges[i].u; ++missing) {
				nb_nodes_so_far++;
				g.slicing_idx[nb_nodes_so_far] = nb_neighbors_so_far;
			}
			// Different node -> mark a slice
			nb_nodes_so_far++;
			g.slicing_idx[nb_nodes_so_far] = nb_neighbors_so_far;
		}

		// Push v as a neighbor of u
		g.neighbors[nb_neighbors_so_far] = edges[i].v;
		g.weights[nb_neighbors_so_far] = edges[i].weight;
		nb_neighbors_so_far++;
	}
	// Close off the slices
	nb_nodes_so_far++;
	g.slicing_idx[nb_nodes_so_far] = nb_neighbors_so_far;
	g.length = nb_neighbors_so_far;
	g.nb_nodes = nb_nodes;
	printf("Took %f seconds to build CSR\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_build).count());

	g.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();

	// Verification avec vieux algorithme
	// graph g2 = from_edge_list_v1(input_list);
	// print_slicing_idx(g);
	// print_slicing_idx(g2);
	// for (int64_t node = 0; node < g.nb_nodes; ++node) {
	// 	print_neighbors_of_node(g, node);
	// 	print_neighbors_of_node(g2, node);
	// }

	free(original_edges_ptr);

	return g;
}

graph from_edge_list_v3_parallel(edge_list input_list) {
	auto start = std::chrono::high_resolution_clock::now();
	graph g;

	// Compute number of nodes
	auto start_count = std::chrono::high_resolution_clock::now();
	int64_t max_node = -1;
#pragma omp parallel for reduction(max : max_node)
	for (int64_t i = 0; i < input_list.length; ++i) {
		max_node = std::max(max_node, input_list.edges[i].v0);
		max_node = std::max(max_node, input_list.edges[i].v1);
	}
	int64_t nb_nodes = max_node + 1;
	printf("Took %f seconds to compute number of nodes\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_count).count());

	// Count degrees (parallel)
	// 	auto start_degree = std::chrono::high_resolution_clock::now();
	// 	std::vector<std::atomic<int64_t>> degrees(nb_nodes);
	// #pragma omp parallel for schedule(static)
	// 	for (int64_t i = 0; i < input_list.length; ++i) {
	// 		// Skip loops
	// 		if (is_loop(input_list.edges[i])) {
	// 			continue;
	// 		}
	// 		degrees[input_list.edges[i].v0].fetch_add(1, std::memory_order_relaxed);
	// 		degrees[input_list.edges[i].v1].fetch_add(1, std::memory_order_relaxed);
	// 	}

	// Count degrees (not parallel) -> was too slow because of atomics or something
	// we also need to remove prefix sum in parallel because it counted on the atomics
	auto start_degree = std::chrono::high_resolution_clock::now();
	std::vector<int64_t> degrees(nb_nodes, 0);
	for (int64_t i = 0; i < input_list.length; ++i) {
		// Skip loops
		if (is_loop(input_list.edges[i])) {
			continue;
		}
		degrees[input_list.edges[i].v0]++;
		degrees[input_list.edges[i].v1]++;
	}
	printf("Took %f seconds to count degrees\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_degree).count());

	// Do prefix sum of degrees to get slicing_idx (parallel scan)
	auto prefix_sum_start = std::chrono::high_resolution_clock::now();

	// g.slicing_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));
	// std::vector<int64_t> partial_sums(omp_get_max_threads() + 1, 0);

	// Compute partial sums in parallel
	// #pragma omp parallel
	// 	{
	// 		int tid = omp_get_thread_num();
	// 		int nthreads = omp_get_num_threads();
	// 		int64_t local_sum = 0;
	// 		int64_t chunk = (nb_nodes + nthreads - 1) / nthreads;
	// 		int64_t start = tid * chunk;
	// 		int64_t end = std::min(nb_nodes, (tid + 1) * chunk);
	// 		for (int64_t i = start; i < end; ++i) {
	// 			local_sum += degrees[i].load(std::memory_order_relaxed);
	// 		}
	// 		partial_sums[tid + 1] = local_sum;
	// #pragma omp barrier
	// 		// Prefix sum of partial_sums (single thread)
	// #pragma omp single
	// 		{
	// 			for (int i = 1; i <= nthreads; ++i) {
	// 				partial_sums[i] += partial_sums[i - 1];
	// 			}
	// 		}
	// #pragma omp barrier
	// 		// Fill slicing_idx in parallel
	// 		int64_t offset = partial_sums[tid];
	// 		for (int64_t i = start; i < end; ++i) {
	// 			g.slicing_idx[i] = offset;
	// 			offset += degrees[i].load(std::memory_order_relaxed);
	// 		}
	// 	}

	// 	g.slicing_idx[nb_nodes] = partial_sums.back();
	// 	int64_t nb_edges = partial_sums.back();

	// Prefix sum (sequential)
	g.slicing_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));
	int64_t sum = 0;
	for (int64_t i = 0; i < nb_nodes; ++i) {
		g.slicing_idx[i] = sum;
		sum += degrees[i];
	}
	g.slicing_idx[nb_nodes] = sum;
	int64_t nb_edges = sum;
	printf("Took %f seconds to compute prefix sum\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prefix_sum_start).count());

	// Allocate neighbor arrays
	g.neighbors = (int64_t*)malloc(nb_edges * sizeof(int64_t));
	g.weights = (float*)malloc(nb_edges * sizeof(float));

	// Fill using per-node cursors (parallel): use atomics to allocate per-node slots
	auto start_fill = std::chrono::high_resolution_clock::now();
	std::vector<std::atomic<int64_t>> next_write_idx(nb_nodes);
	for (int64_t i = 0; i < nb_nodes; ++i) {
		next_write_idx[i].store(g.slicing_idx[i], std::memory_order_relaxed);
	}
#pragma omp parallel for schedule(static)
	for (int64_t i = 0; i < input_list.length; ++i) {
		int64_t u = input_list.edges[i].v0;
		int64_t v = input_list.edges[i].v1;

		// Skip loops
		if (is_loop(input_list.edges[i])) {
			continue;
		}

		float weight = input_list.weights[i];

		// Atomically reserve a slot for u and v
		int64_t write_idx_u = next_write_idx[u].fetch_add(1, std::memory_order_relaxed);
		g.neighbors[write_idx_u] = v;
		g.weights[write_idx_u] = weight;

		int64_t write_idx_v = next_write_idx[v].fetch_add(1, std::memory_order_relaxed);
		g.neighbors[write_idx_v] = u;
		g.weights[write_idx_v] = weight;
	}
	printf("Took %f seconds to fill neighbor arrays\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_fill).count());

	g.length = nb_edges;
	g.nb_nodes = nb_nodes;

	// On-the-fly per-node deduplication using a marker table
	auto start_dedup = std::chrono::high_resolution_clock::now();
	int64_t old_nb_edges = nb_edges;
	int64_t* new_neighbors = (int64_t*)malloc(old_nb_edges * sizeof(int64_t));
	float* new_weights = (float*)malloc(old_nb_edges * sizeof(float));
	int64_t* new_slices_idx = (int64_t*)malloc((nb_nodes + 1) * sizeof(int64_t));

	std::vector<int64_t> marker(nb_nodes, -1); // marker[v] = position in new_neighbors or -1
	std::vector<int64_t> touched;			   // list of nodes whose marker was touched, to reset later
	touched.reserve(64);

	int64_t write_pos = 0;
	for (int64_t node = 0; node < nb_nodes; node++) {
		new_slices_idx[node] = write_pos;
		for_each_neighbor(g, node, [&](int64_t neighbor, float weight) {
			if (marker[neighbor] == -1) {
				marker[neighbor] = write_pos;
				new_neighbors[write_pos] = neighbor;
				new_weights[write_pos] = weight;
				touched.push_back(neighbor);
				write_pos++;
			}
			else {
				// Duplicate neighbor, just keep first weight encountered
			}
		});
		// Reset markers for touched neighbors of this node
		for (int64_t node : touched) {
			marker[node] = -1;
		}
		touched.clear();
	}

	new_slices_idx[nb_nodes] = write_pos;

	// replace arrays
	free(g.neighbors);
	free(g.weights);
	free(g.slicing_idx);
	g.neighbors = new_neighbors;
	g.weights = new_weights;
	g.slicing_idx = new_slices_idx;
	g.length = write_pos;

	printf("Took %f seconds to deduplicate neighbors\n",
		   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_dedup).count());

	g.time_ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();

	// Verification avec vieux algorithme
	// graph g_ref = from_edge_list_v2(input_list);
	// print_slicing_idx(g);
	// print_slicing_idx(g_ref);
	// for (int64_t node = 0; node < g.nb_nodes; ++node) {
	// 	print_neighbors_of_node(g, node);
	// 	print_neighbors_of_node(g_ref, node);
	// }

	return g;
}
#include <omp.h>

void from_edge_list_try_all(edge_list input_list) {
	graph g;

	// No need to repeat the sequential experiment
	if (omp_get_max_threads() == 1) {
		g = from_edge_list_v2(input_list);
		printf("From edge list v2: %fms\n", g.time_ms);
		graph_destroy(g);

		g = from_edge_list_v3(input_list);
		printf("From edge list v3: %fms\n", g.time_ms);
		graph_destroy(g);
	}
	else {
		g = from_edge_list_v2_parallel(input_list);
		printf("From edge list v2 parallel: %fms\n", g.time_ms);
		graph_destroy(g);

		g = from_edge_list_v3_parallel(input_list);
		printf("From edge list v3 parallel: %fms\n", g.time_ms);
		graph_destroy(g);
	}
}

graph from_edge_list(edge_list input_list) {
	return from_edge_list_v2(input_list);
}