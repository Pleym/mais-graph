#include "generator/generator.hpp"
#include "kernels/breadth_search.hpp"
#include "kernels/gen_graph.hpp"
#include "kernels/shortest_path.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <random>

int main(int argc, char const* argv[]) {
	uint8_t scale = 12;
	uint8_t edge_factor = 16;
	int64_t roots = 16;

	if (argc > 1) {
		scale = static_cast<uint8_t>(std::atoi(argv[1]));
	}
	if (argc > 2) {
		edge_factor = static_cast<uint8_t>(std::atoi(argv[2]));
	}
	if (argc > 3) {
		roots = std::max<int64_t>(1, std::atoll(argv[3]));
	}

	std::cout << "Running reference Graph500 kernels with scale=" << static_cast<int>(scale)
			  << ", edge_factor=" << static_cast<int>(edge_factor)
			  << ", roots=" << roots << std::endl;

	edge_list list = generate_graph(scale, edge_factor);
	graph g = from_edge_list(list);
	edge_list_destroy(list);

	std::cout << "K1 (graph): nodes=" << g.nb_nodes
			  << ", edges=" << g.length
			  << ", time=" << g.time_ms << " ms" << std::endl;

	double k2_time_ms = 0.0;
	double k3_time_ms = 0.0;
	double k3_teps = 0.0;

	std::mt19937_64 rng(42);
	std::uniform_int_distribution<int64_t> dist(0, g.nb_nodes - 1);

	for (int64_t i = 0; i < roots;) {
		int64_t source = dist(rng);
		if (degree_of_node(g, source) == 0) {
			continue;
		}
		i++;

		bfs_result bfs_ref = bfs_formal(g, source);
		k2_time_ms += bfs_ref.time_ms;

		shortest_path sssp_ref = sssp_bf(g, source);
		k3_time_ms += sssp_ref.time_ms;
		k3_teps += sssp_ref.teps;
		shortest_path_destroy(sssp_ref);
	}

	k2_time_ms /= roots;
	k3_time_ms /= roots;
	k3_teps /= roots;

	std::cout << "K2 (bfs reference): avg_time=" << k2_time_ms << " ms" << std::endl;
	std::cout << "K3 (sssp reference BF): avg_time=" << k3_time_ms << " ms, avg_teps=" << k3_teps << std::endl;

	graph_destroy(g);
	return 0;
}