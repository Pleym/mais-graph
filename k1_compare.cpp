#include "generator/generator.hpp"
#include "kernels/gen_graph.hpp"
#include "kernels/breadth_search.hpp"
#include "kernels/shortest_path.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <omp.h>

struct KernelStats {
	double k2_ms;
	double k3_ms;
	double teps;
};

static KernelStats run_k2k3(graph& g, int64_t roots, const std::string& mode) {
	std::mt19937_64 rng(42);
	std::uniform_int_distribution<int64_t> dist(0, g.nb_nodes - 1);
	double k2_time_ms = 0.0;
	double k3_time_ms = 0.0;
	double k3_teps = 0.0;

	for (int64_t i = 0; i < roots;) {
		int64_t source = dist(rng);
		if (degree_of_node(g, source) == 0) {
			continue;
		}
		++i;

		bfs_result r;
		shortest_path s;
		if (mode == "ref") {
			r = bfs_formal(g, source);
			s = sssp_bf(g, source);
		}
		else if (mode == "bfs_opt") {
			r = bfs_top_down_omp_cas(g, source);
			s = sssp_bf(g, source);
		}
		else if (mode == "sssp_opt") {
			r = bfs_formal(g, source);
			s = sssp_bf_frontier_omp(g, source);
		}
		else {
			// all_opt or hybrid fall back to the fastest combo
			r = bfs_top_down_omp_cas(g, source);
			s = sssp_bf_frontier_omp(g, source);
		}

		k2_time_ms += r.time_ms;
		k3_time_ms += s.time_ms;
		k3_teps += s.teps;

		free(r.parent_array);
		shortest_path_destroy(s);
	}

	return {
		.k2_ms = k2_time_ms / roots,
		.k3_ms = k3_time_ms / roots,
		.teps = k3_teps / roots,
	};
}

int main(int argc, char** argv) {
	int scale = (argc > 1) ? std::atoi(argv[1]) : 12;
	int edge_factor = (argc > 2) ? std::atoi(argv[2]) : 16;
	int64_t roots = (argc > 3) ? std::atoll(argv[3]) : 4;
	std::string mode = (argc > 4) ? argv[4] : "all_opt";

	std::cout << "K1/K2/K3 comparison: scale=" << scale << " edge_factor=" << edge_factor << " roots=" << roots
			  << " mode=" << mode << " threads=" << omp_get_max_threads() << "\n";

	auto t0 = std::chrono::steady_clock::now();
	edge_list list = generate_graph(scale, edge_factor);
	auto t1 = std::chrono::steady_clock::now();
	double edge_gen_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

	std::cout << "Edge list generated in " << edge_gen_ms << " ms (" << list.length << " directed edges)\n";

	graph g_v2 = from_edge_list_v2(list);
	KernelStats stats_v2 = run_k2k3(g_v2, roots, mode);
	std::cout << "v2_sort K1_ms=" << g_v2.time_ms << " K2_ms=" << stats_v2.k2_ms << " K3_ms=" << stats_v2.k3_ms
			  << " K3_teps=" << stats_v2.teps << "\n";
	graph_destroy(g_v2);

	graph g_v3 = from_edge_list_v3_parallel(list);
	KernelStats stats_v3 = run_k2k3(g_v3, roots, mode);
	std::cout << "v3_parallel K1_ms=" << g_v3.time_ms << " K2_ms=" << stats_v3.k2_ms << " K3_ms=" << stats_v3.k3_ms
			  << " K3_teps=" << stats_v3.teps << "\n";
	graph_destroy(g_v3);

	edge_list_destroy(list);
	return 0;
}
