#include "breadth_search.hpp"
#include "bitset.hpp"
#include "gen_graph.hpp"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <queue>
#include <unordered_set>
#include <vector>

static bool verify_bfs_result(const graph& g, int64_t source, int64_t* parents) {
	if (parents[source] != source) {
		std::fprintf(stderr, "verify_bfs_result: source %lld has parent %lld\n", (long long)source, (long long)parents[source]);
		return false;
	}

	// Verify that the parent and child are connected and they have valid IDs
	for (int64_t node = 0; node < g.nb_nodes; node++) {
		if (node == source) {
			continue;
		}

		int64_t parent = parents[node];

		// If the parent is -1 -> node wasn't able to be reached
		if (parent == -1) {
			continue;
		}
		// Verify that the parent has a valid ID
		if (parent < 0 || parent >= g.nb_nodes) {
			std::fprintf(stderr, "verify_bfs_result: node %lld has out-of-range parent %lld\n", (long long)node, (long long)parent);
			return false;
		}

		// Verify that the parent and child are connected
		bool found = false;
		for_each_neighbor(g, parent, [&](int64_t neigh, float) {
			if (neigh == node) {
				found = true;
			}
		});
		if (!found) {
			std::fprintf(
				stderr, "verify_bfs_result: parent %lld of node %lld is not connected to it\n", (long long)parent, (long long)node);
			return false;
		}
	}

	// Ensure no cycles in the result. Treat a node whose parent equals itself
	// as a root (terminal), so it does not form a cycle.
	typedef enum : uint8_t {
		UNVISITED,
		VISITING,
		DONE,
	} cycles_verification;
	std::vector<cycles_verification> state(g.nb_nodes, UNVISITED);
	for (int64_t i = 0; i < g.nb_nodes; ++i) {
		if (parents[i] == -1 || state[i] != UNVISITED) {
			continue;
		}
		int64_t current = i;
		while (current != -1 && state[current] == UNVISITED) {
			state[current] = VISITING;
			int64_t next = parents[current];
			if (next == current) { // self-parent -> treat as root
				current = -1;
				break;
			}
			current = next;
		}

		if (current != -1 && state[current] == VISITING) {
			// collect nodes in the cycle
			std::vector<int64_t> cycle;
			int64_t start = current;
			cycle.push_back(start);
			int64_t cur = parents[start];
			while (cur != start && cur != -1) {
				cycle.push_back(cur);
				cur = parents[cur];
			}
			std::fprintf(stderr, "verify_bfs_result: cycle detected (length %zu):", cycle.size());
			for (size_t k = 0; k < cycle.size(); ++k) {
				std::fprintf(stderr, " %lld", (long long)cycle[k]);
			}
			std::fprintf(stderr, "\n");
			return false;
		}

		// mark the chain as done (stop at self-parent too)
		current = i;
		while (current != -1 && state[current] == VISITING) {
			int64_t next = parents[current];
			state[current] = DONE;
			if (next == current)
				break;
			current = next;
		}
	}

	return true;
}

bfs_result bfs_formal(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	int64_t* parents = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	std::fill(parents, parents + g.nb_nodes, -1);
	parents[source] = source;

	bool* visited = (bool*)malloc(g.nb_nodes * sizeof(bool));
	std::fill(visited, visited + g.nb_nodes, false);

	std::queue<int64_t> queue;
	queue.push(source);
	visited[source] = true;

	while (!queue.empty()) {
		int64_t node = queue.front();
		queue.pop();
		for_each_neighbor(g, node, [&](int64_t neighbor, float _weight) {
			if (!visited[neighbor]) {
				visited[neighbor] = true;
				parents[neighbor] = node;
				queue.push(neighbor);
			}
		});
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	free(visited);

	return (bfs_result){
		.parent_array = parents,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

bfs_result bfs_top_down_omp_cas(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	std::vector<int64_t> frontier;
	std::vector<int64_t> next_frontier;
	frontier.push_back(source);

	std::atomic<int64_t>* parents_atomic = new std::atomic<int64_t>[g.nb_nodes];
	for (int64_t i = 0; i < g.nb_nodes; ++i) {
		parents_atomic[i].store(-1, std::memory_order_relaxed);
	}
	parents_atomic[source].store(source, std::memory_order_relaxed);

	while (!frontier.empty()) {
		int max_threads = omp_get_max_threads();
		std::vector<std::vector<int64_t>> local_next(max_threads);

#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			auto& local = local_next[tid];

#pragma omp for schedule(dynamic, 64)
			for (size_t idx = 0; idx < frontier.size(); ++idx) {
				int64_t node = frontier[idx];
				for_each_neighbor(g, node, [&](int64_t neighbor, float _weight) {
					int64_t expected = -1;
					if (parents_atomic[neighbor].compare_exchange_strong(
							expected,
							node,
							std::memory_order_relaxed,
							std::memory_order_relaxed)) {
						local.push_back(neighbor);
					}
				});
			}
		}

		next_frontier.clear();
		size_t total = 0;
		for (const auto& v : local_next) {
			total += v.size();
		}
		next_frontier.reserve(total);
		for (auto& v : local_next) {
			next_frontier.insert(next_frontier.end(), v.begin(), v.end());
		}

		frontier.swap(next_frontier);
	}

	int64_t* parents = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	for (int64_t i = 0; i < g.nb_nodes; ++i) {
		parents[i] = parents_atomic[i].load(std::memory_order_relaxed);
	}
	delete[] parents_atomic;

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return (bfs_result){
		.parent_array = parents,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

///////////////////////////////////////////////////////////////
//////////// UNORDERED SET FRONTIER IMPLEMENTATION ////////////
///////////////////////////////////////////////////////////////

static void top_down_step(graph& g, std::unordered_set<int64_t>& frontier, std::unordered_set<int64_t>& next, int64_t*& parents) {
	for (auto node : frontier) {
		for_each_neighbor(g, node, [&](int64_t neighbor, float _weight) {
			if (parents[neighbor] == -1) {
				parents[neighbor] = node;
				next.insert(neighbor);
			}
		});
	}
}

static void bottom_up_step(graph& g, std::unordered_set<int64_t>& frontier, std::unordered_set<int64_t>& next, int64_t*& parents) {
	for (int64_t node = 0; node < g.nb_nodes; node++) {
		if (parents[node] == -1) {
			int64_t start = g.slicing_idx[node];
			int64_t end = g.slicing_idx[node + 1];
			for (int64_t i = start; i < end; i++) {
				int64_t neighbor = g.neighbors[i];
				if (frontier.find(neighbor) != frontier.end()) {
					parents[node] = neighbor;
					next.insert(node);
					break;
				}
			}
		}
	}
}

bfs_result bfs_full_top_down(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	auto frontier = std::unordered_set<int64_t>{};
	frontier.insert(source);
	auto next = std::unordered_set<int64_t>{};
	int64_t* parents = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	std::fill(parents, parents + g.nb_nodes, -1);
	parents[source] = source;

	while (!frontier.empty()) {
		top_down_step(g, frontier, next, parents);
		frontier = next;
		next.clear();
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return (bfs_result){
		.parent_array = parents,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

bfs_result bfs_full_bottom_up(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	auto frontier = std::unordered_set<int64_t>{};
	frontier.insert(source);
	auto next = std::unordered_set<int64_t>{};
	int64_t* parents = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	std::fill(parents, parents + g.nb_nodes, -1);
	parents[source] = source;

	while (!frontier.empty()) {
		bottom_up_step(g, frontier, next, parents);
		frontier = next;
		next.clear();
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return (bfs_result){
		.parent_array = parents,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

////////////////////////////////////////////////////////
//////////// BITSET FRONTIER IMPLEMENTATION ////////////
////////////////////////////////////////////////////////

static void top_down_step_bitset(graph& g, BitSet& frontier, BitSet& next, int64_t*& parents) {
	frontier.for_each([&](int64_t node) {
		for_each_neighbor(g, node, [&](int64_t neighbor, float _weight) {
			if (parents[neighbor] == -1) {
				parents[neighbor] = node;
				next.insert(neighbor);
			}
		});
	});
}

static void bottom_up_step_bitset(graph& g, BitSet& frontier, BitSet& next, int64_t*& parents) {
	for (int64_t node = 0; node < g.nb_nodes; node++) {
		if (parents[node] == -1) {
			int64_t start = g.slicing_idx[node];
			int64_t end = g.slicing_idx[node + 1];
			for (int64_t i = start; i < end; i++) {
				int64_t neighbor = g.neighbors[i];
				if (frontier.contains(neighbor)) {
					parents[node] = neighbor;
					next.insert(node);
					break;
				}
			}
		}
	}
}

bfs_result bfs_full_top_down_bitset(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	auto frontier = BitSet(g.nb_nodes);
	auto next = BitSet(g.nb_nodes);
	int64_t* parents = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	std::fill(parents, parents + g.nb_nodes, -1);
	parents[source] = source;
	frontier.insert(source);
	while (!frontier.empty()) {
		top_down_step_bitset(g, frontier, next, parents);
		bitset_swap(frontier, next);
		next.clear();
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return (bfs_result){
		.parent_array = parents,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

bfs_result bfs_full_bottom_up_bitset(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	auto frontier = BitSet(g.nb_nodes);
	auto next = BitSet(g.nb_nodes);
	int64_t* parents = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	std::fill(parents, parents + g.nb_nodes, -1);
	parents[source] = source;
	frontier.insert(source);
	while (!frontier.empty()) {
		bottom_up_step_bitset(g, frontier, next, parents);
		bitset_swap(frontier, next);
		next.clear();
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return (bfs_result){
		.parent_array = parents,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

////////////////////////////////////////////////////////
/////// PARALLEL BITSET FRONTIER IMPLEMENTATION ////////
////////////////////////////////////////////////////////

static void top_down_step_parallel_bitset(graph& g, AtomicBitSet& frontier, AtomicBitSet& next, std::atomic<int64_t>*& parents) {
	frontier.parallel_for_each([&](int64_t node) {
		for_each_neighbor(g, node, [&](int64_t neighbor, float _weight) {
			if (parents[neighbor].load(std::memory_order_acquire) == -1) {
				parents[neighbor].store(node, std::memory_order_release);
				next.insert(neighbor);
			}
		});
	});
}

static void bottom_up_step_parallel_bitset(graph& g, AtomicBitSet& frontier, AtomicBitSet& next, std::atomic<int64_t>*& parents) {
#pragma omp parallel for schedule(dynamic, 1024)
	for (int64_t node = 0; node < g.nb_nodes; node++) {
		if (parents[node].load(std::memory_order_acquire) == -1) {
			int64_t start = g.slicing_idx[node];
			int64_t end = g.slicing_idx[node + 1];
			for (int64_t i = start; i < end; i++) {
				int64_t neighbor = g.neighbors[i];
				if (frontier.contains(neighbor)) {
					parents[node].store(neighbor, std::memory_order_release);
					next.insert(node);
					break;
				}
			}
		}
	}
}

bfs_result bfs_full_top_down_parallel_bitset(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	auto frontier = AtomicBitSet(g.nb_nodes);
	auto next = AtomicBitSet(g.nb_nodes);
	std::atomic<int64_t>* parents = new std::atomic<int64_t>[g.nb_nodes];
	for (int64_t i = 0; i < g.nb_nodes; i++) {
		parents[i].store(-1, std::memory_order_release);
	}
	parents[source].store(source, std::memory_order_release);
	frontier.insert(source);
	while (!frontier.empty()) {
		top_down_step_parallel_bitset(g, frontier, next, parents);
		atomic_bitset_swap(frontier, next);
		next.clear();
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	// Turn atomic parents into regular int64_t array for easier handling by caller, not timing it because we have the result already in the
	// atomic array
	int64_t* parents_array = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	for (int64_t i = 0; i < g.nb_nodes; i++) {
		parents_array[i] = parents[i].load(std::memory_order_acquire);
	}
	delete[] parents;

	return (bfs_result){
		.parent_array = parents_array,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

bfs_result bfs_full_bottom_up_parallel_bitset(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	auto frontier = AtomicBitSet(g.nb_nodes);
	auto next = AtomicBitSet(g.nb_nodes);
	std::atomic<int64_t>* parents = new std::atomic<int64_t>[g.nb_nodes];
	for (int64_t i = 0; i < g.nb_nodes; i++) {
		parents[i].store(-1, std::memory_order_release);
	}
	parents[source].store(source, std::memory_order_release);
	frontier.insert(source);
	while (!frontier.empty()) {
		bottom_up_step_parallel_bitset(g, frontier, next, parents);
		atomic_bitset_swap(frontier, next);
		next.clear();
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	// Turn atomic parents into regular int64_t array for same reason as in top-down version
	int64_t* parents_array = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	for (int64_t i = 0; i < g.nb_nodes; i++) {
		parents_array[i] = parents[i].load(std::memory_order_acquire);
	}
	delete[] parents;

	return (bfs_result){
		.parent_array = parents_array,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

typedef struct hybrid_set {
	enum FrontierType : uint8_t {
		UNORDERED_SET,
		BITSET,
	} currently_using;

	std::unordered_set<int64_t> unordered_set;
	AtomicBitSet bitset;

	hybrid_set(int64_t nb_nodes) : currently_using(UNORDERED_SET), unordered_set(), bitset(nb_nodes) {}

	std::unordered_set<int64_t>& as_unordered_set() {
		if (this->currently_using == FrontierType::UNORDERED_SET) {
			return unordered_set;
		}
		else {
			// convert bitset to unordered_set
			unordered_set.clear();
			bitset.for_each([&](int64_t v) { unordered_set.insert(v); });
			this->currently_using = FrontierType::UNORDERED_SET;
			return unordered_set;
		}
	}

	AtomicBitSet& as_bitset() {
		if (this->currently_using == FrontierType::BITSET) {
			return bitset;
		}
		else {
			// convert unordered_set to bitset
			bitset.clear();
			for (auto v : unordered_set) {
				bitset.insert(v);
			}
			this->currently_using = FrontierType::BITSET;
			return bitset;
		}
	}

	bool empty() const {
		if (this->currently_using == FrontierType::UNORDERED_SET) {
			return unordered_set.empty();
		}
		else {
			return bitset.empty();
		}
	}

	int64_t sum_of_degrees(const graph& g) const {
		int64_t sum = 0;
		if (this->currently_using == FrontierType::UNORDERED_SET) {
			for (auto v : unordered_set) {
				sum += (int64_t)(g.slicing_idx[v + 1] - g.slicing_idx[v]);
			}
		}
		else {
			bitset.for_each([&](int64_t v) { sum += (int64_t)(g.slicing_idx[v + 1] - g.slicing_idx[v]); });
		}
		return sum;
	}
} hybrid_set;

// Hybrid direction-optimizing BFS using bitset frontier.
// Switches between top-down and bottom-up using a simple heuristic based
// on the number of edges to check from the frontier (m_f) and from
// unexplored vertices (m_u). Constants C_TB and C_BT tune switching.
bfs_result bfs_hybrid(graph& g, int64_t source) {
	auto start = std::chrono::high_resolution_clock::now();

	// Hybrid with vector frontier + AtomicBitSet and incremental m_u
	std::vector<int64_t> frontier_vec;
	frontier_vec.push_back(source);
	AtomicBitSet frontier_bits(g.nb_nodes);
	frontier_bits.insert(source);

	std::vector<int64_t> next_vec;
	AtomicBitSet next_bits(g.nb_nodes);

	std::atomic<int64_t>* parents = new std::atomic<int64_t>[g.nb_nodes];
	for (int64_t i = 0; i < g.nb_nodes; ++i) {
		parents[i].store(-1, std::memory_order_relaxed);
	}
	parents[source].store(source, std::memory_order_relaxed);

	std::vector<char> visited(g.nb_nodes, 0);
	visited[(size_t)source] = 1;

	std::atomic<int64_t> remaining_edges{g.slicing_idx[g.nb_nodes]};

	bool top_down = true;
	const double C_TB = 10.0; // threshold for switching top->bottom
	const double C_BT = 40.0; // threshold for switching bottom->top

	while (!frontier_vec.empty()) {
		uint64_t m_f = 0;
		for (int64_t v : frontier_vec) {
			m_f += (uint64_t)(g.slicing_idx[v + 1] - g.slicing_idx[v]);
		}
		uint64_t m_u = (uint64_t)remaining_edges.load(std::memory_order_relaxed);

		if (top_down) {
			if (m_f > m_u / C_TB) {
				top_down = false;
			}
		}
		else {
			if (m_f < m_u / C_BT) {
				top_down = true;
			}
		}

		if (top_down) {
			int max_threads = omp_get_max_threads();
			std::vector<std::vector<int64_t>> locals(max_threads);

#pragma omp parallel
			{
				int tid = omp_get_thread_num();
				auto& local = locals[tid];

#pragma omp for schedule(dynamic, 64)
				for (size_t idx = 0; idx < frontier_vec.size(); ++idx) {
					int64_t node = frontier_vec[idx];
					for_each_neighbor(g, node, [&](int64_t neighbor, float) {
						int64_t expected = -1;
						if (parents[neighbor].compare_exchange_strong(expected, node, std::memory_order_relaxed)) {
							if (!visited[(size_t)neighbor]) {
								visited[(size_t)neighbor] = 1;
								remaining_edges.fetch_sub((int64_t)(g.slicing_idx[neighbor + 1] - g.slicing_idx[neighbor]), std::memory_order_relaxed);
							}
							local.push_back(neighbor);
						}
					});
				}
			}

			size_t total = 0;
			for (auto& v : locals) total += v.size();
			next_vec.clear();
			next_vec.reserve(total);
			for (auto& v : locals) {
				next_vec.insert(next_vec.end(), v.begin(), v.end());
			}

			frontier_bits.clear();
			for (int64_t v : frontier_vec) {
				frontier_bits.insert(v);
			}

			next_bits.clear();
			for (int64_t v : next_vec) {
				next_bits.insert(v);
			}
		}
		else {
			int max_threads = omp_get_max_threads();
			std::vector<std::vector<int64_t>> locals(max_threads);

#pragma omp parallel for schedule(dynamic, 1024)
			for (int64_t node = 0; node < g.nb_nodes; ++node) {
				if (parents[node].load(std::memory_order_relaxed) != -1) {
					continue;
				}
				int64_t start = g.slicing_idx[node];
				int64_t end = g.slicing_idx[node + 1];
				for (int64_t i = start; i < end; ++i) {
					int64_t neigh = g.neighbors[i];
					if (frontier_bits.contains(neigh)) {
						int64_t expected = -1;
						if (parents[node].compare_exchange_strong(expected, neigh, std::memory_order_relaxed)) {
							if (!visited[(size_t)node]) {
								visited[(size_t)node] = 1;
								remaining_edges.fetch_sub((int64_t)(end - start), std::memory_order_relaxed);
							}
							int tid = omp_get_thread_num();
							locals[tid].push_back(node);
						}
						break;
					}
				}
			}

			size_t total = 0;
			for (auto& v : locals) total += v.size();
			next_vec.clear();
			next_vec.reserve(total);
			for (auto& v : locals) {
				next_vec.insert(next_vec.end(), v.begin(), v.end());
			}

			next_bits.clear();
			for (int64_t v : next_vec) {
				next_bits.insert(v);
			}
			frontier_bits = next_bits;
		}

		frontier_vec.swap(next_vec);
	}

	int64_t* parents_out = (int64_t*)malloc(g.nb_nodes * sizeof(int64_t));
	for (int64_t i = 0; i < g.nb_nodes; ++i) {
		parents_out[i] = parents[i].load(std::memory_order_relaxed);
	}
	delete[] parents;

	auto end = std::chrono::high_resolution_clock::now();
	auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	return (bfs_result){
		.parent_array = parents_out,
		.teps = 0.0,
		.time_ms = time_ms,
	};
}

bfs_result bfs(graph& g, int64_t source) {
	bfs_result result = bfs_formal(g, source);
	bool is_correct = verify_bfs_result(g, source, result.parent_array);
	if (!is_correct) {
		printf("Error in formal BFS kernel for node %lu\n", source);
	}
	printf("Formal BFS: time = %.2f ms\n", result.time_ms);
	return result;
}