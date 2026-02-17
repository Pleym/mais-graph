#include "generator/generator.hpp"
#include "kernels/breadth_search.hpp"
#include "kernels/gen_graph.hpp"
#include "kernels/shortest_path.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>

int main(int argc, char const* argv[]) {
	uint8_t scale = 12;
	uint8_t edge_factor = 16;
	int64_t roots = 16;
	std::string mode = "ref";

	if (argc > 1) {
		scale = static_cast<uint8_t>(std::atoi(argv[1]));
	}
	if (argc > 2) {
		edge_factor = static_cast<uint8_t>(std::atoi(argv[2]));
	}
	if (argc > 3) {
		roots = std::max<int64_t>(1, std::atoll(argv[3]));
	}
	if (argc > 4) {
		mode = argv[4];
	}

	std::cout << "Running Graph500 kernels with scale=" << static_cast<int>(scale)
			  << ", edge_factor=" << static_cast<int>(edge_factor)
			  << ", roots=" << roots
			  << ", mode=" << mode << std::endl;

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

		if (mode == "ref") {
			bfs_result r = bfs_formal(g, source);
			k2_time_ms += r.time_ms;

			shortest_path s = sssp_bf(g, source);
			k3_time_ms += s.time_ms;
			k3_teps += s.teps;
			shortest_path_destroy(s);
		}
		else if (mode == "bfs_opt") {
			bfs_result r = bfs_top_down_omp_cas(g, source);
			k2_time_ms += r.time_ms;

			shortest_path s = sssp_bf(g, source);
			k3_time_ms += s.time_ms;
			k3_teps += s.teps;
			shortest_path_destroy(s);
		}
		else if (mode == "sssp_opt") {
			bfs_result r = bfs_formal(g, source);
			k2_time_ms += r.time_ms;

			shortest_path s = sssp_bf_frontier_omp(g, source);
			k3_time_ms += s.time_ms;
			k3_teps += s.teps;
			shortest_path_destroy(s);
		}
		else if (mode == "all_opt") {
			bfs_result r = bfs_top_down_omp_cas(g, source);
			k2_time_ms += r.time_ms;

			shortest_path s = sssp_bf_frontier_omp(g, source);
			k3_time_ms += s.time_ms;
			k3_teps += s.teps;
			shortest_path_destroy(s);
		}
		else if (mode == "hybrid") {
			bfs_result r = bfs_hybrid(g, source);
			k2_time_ms += r.time_ms;

			shortest_path s = sssp_bf_frontier_omp(g, source);
			k3_time_ms += s.time_ms;
			k3_teps += s.teps;
			shortest_path_destroy(s);
		}
		else {
			std::cerr << "Unknown mode: " << mode << " (supported: ref, bfs_opt, sssp_opt, all_opt, hybrid)" << std::endl;
			graph_destroy(g);
			return 2;
		}
	}

	k2_time_ms /= roots;
	k3_time_ms /= roots;
	k3_teps /= roots;

	std::cout << "K2 (bfs " << mode << "): avg_time=" << k2_time_ms << " ms" << std::endl;
	std::cout << "K3 (sssp " << mode << "): avg_time=" << k3_time_ms << " ms, avg_teps=" << k3_teps << std::endl;

	graph_destroy(g);
	return 0;
}