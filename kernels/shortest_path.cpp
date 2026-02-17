#include "shortest_path.hpp"
#include <assert.h>
#include <chrono>
#include <cstring>
#include <limits>
#include <omp.h>
#include <queue>

void shortest_path_destroy(shortest_path& sp) {
	free(sp.distance_array);
	free(sp.parent_array);
}

shortest_path sssp_dj(graph g, int64_t root) {
	auto start = std::chrono::high_resolution_clock::now();

	const float INF = std::numeric_limits<float>::infinity();
	int64_t N = g.nb_nodes;

	shortest_path sp;
	sp.distance_array = (float*)malloc(N * sizeof(float));
	sp.parent_array = (int64_t*)malloc(N * sizeof(int64_t));

	for (int64_t i = 0; i < N; i++) {
		sp.distance_array[i] = INF;
		sp.parent_array[i] = -1;
	}

	sp.distance_array[root] = 0.0f;
	sp.parent_array[root] = root;

	// (distance, node)
	using pq_elem = std::pair<float, int64_t>;
	std::priority_queue<pq_elem, std::vector<pq_elem>, std::greater<pq_elem>> pq;

	pq.push({0.0f, root});

	while (!pq.empty()) {
		auto [dist_u, u] = pq.top();
		pq.pop();

		// skip si pas la meilleure distance
		if (dist_u > sp.distance_array[u])
			continue;

		for (int64_t i = g.slicing_idx[u]; i < g.slicing_idx[u + 1]; i++) {
			int64_t v = g.neighbors[i];
			float w = g.weights[i];

			float alt = dist_u + w;
			if (alt < sp.distance_array[v]) {
				sp.distance_array[v] = alt;
				sp.parent_array[v] = u;
				pq.push({alt, v});
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	sp.time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	sp.teps = (double)g.length / (sp.time_ms * 1e-3);

	return sp;
}

shortest_path sssp_bf(graph g, int64_t root) {
    auto start = std::chrono::high_resolution_clock::now();

    const float INF = std::numeric_limits<float>::infinity();
    int64_t N = g.nb_nodes;

    shortest_path sp;
    sp.distance_array = (float*)malloc(N * sizeof(float));
    sp.parent_array   = (int64_t*)malloc(N * sizeof(int64_t));

    float* dist_old = (float*)malloc(N * sizeof(float));
    float* dist_new = (float*)malloc(N * sizeof(float));

    for (int64_t i = 0; i < N; i++) {
        dist_old[i] = INF;
        dist_new[i] = INF;
        sp.parent_array[i] = -1;
    }

    dist_old[root] = 0.0f;
    sp.parent_array[root] = root;

    bool changed = true;

    for (int64_t iter = 0; iter < N - 1 && changed; iter++) {
        changed = false;

        std::memcpy(dist_new, dist_old, N * sizeof(float));

        for (int64_t u = 0; u < N; u++) {
            float du = dist_old[u];
            if (du == INF) continue;

            for (int64_t i = g.slicing_idx[u]; i < g.slicing_idx[u + 1]; i++) {
                int64_t v = g.neighbors[i];
                float w = g.weights[i];

                float alt = du + w;
                if (alt < dist_new[v]) {
                    dist_new[v] = alt;
                    sp.parent_array[v] = u;
                    changed = true;
                }
            }
        }

        std::swap(dist_old, dist_new);
    }

    std::memcpy(sp.distance_array, dist_old, N * sizeof(float));

    free(dist_old);
    free(dist_new);

    auto end = std::chrono::high_resolution_clock::now();
    sp.time_ms =
        std::chrono::duration<double, std::milli>(end - start).count();

    sp.teps = (double)g.length / (sp.time_ms * 1e-3);

    return sp;
}

shortest_path sssp_bf_frontier_omp(graph g, int64_t root) {
	auto start = std::chrono::high_resolution_clock::now();

	const float INF = std::numeric_limits<float>::infinity();
	int64_t N = g.nb_nodes;

	shortest_path sp;
	sp.distance_array = (float*)malloc(N * sizeof(float));
	sp.parent_array = (int64_t*)malloc(N * sizeof(int64_t));

	std::atomic<float>* dist = new std::atomic<float>[N];
	std::atomic<int64_t>* parent = new std::atomic<int64_t>[N];
	for (int64_t i = 0; i < N; ++i) {
		dist[i].store(INF, std::memory_order_relaxed);
		parent[i].store(-1, std::memory_order_relaxed);
	}
	dist[root].store(0.0f, std::memory_order_relaxed);
	parent[root].store(root, std::memory_order_relaxed);

	std::vector<int64_t> frontier;
	frontier.push_back(root);

	std::atomic<int>* seen_epoch = new std::atomic<int>[N];
	for (int64_t i = 0; i < N; ++i) {
		seen_epoch[i].store(-1, std::memory_order_relaxed);
	}
	int epoch = 0;

	for (int64_t iter = 0; iter < N - 1 && !frontier.empty(); ++iter, ++epoch) {
		std::vector<std::vector<int64_t>> local_next(omp_get_max_threads());

#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			auto& local = local_next[tid];

#pragma omp for schedule(dynamic, 64)
			for (size_t fi = 0; fi < frontier.size(); ++fi) {
				int64_t u = frontier[fi];
				float du = dist[u].load(std::memory_order_relaxed);
				if (du == INF) {
					continue;
				}

				int64_t begin = g.slicing_idx[u];
				int64_t end = g.slicing_idx[u + 1];
				for (int64_t ei = begin; ei < end; ++ei) {
					int64_t v = g.neighbors[ei];
					float alt = du + g.weights[ei];

					float cur = dist[v].load(std::memory_order_relaxed);
					bool updated = false;
					while (alt < cur) {
						if (dist[v].compare_exchange_weak(cur, alt, std::memory_order_relaxed)) {
							updated = true;
							break;
						}
					}

					if (updated) {
						parent[v].store(u, std::memory_order_relaxed);
						int prev = seen_epoch[v].exchange(epoch, std::memory_order_relaxed);
						if (prev != epoch) {
							local.push_back(v);
						}
					}
				}
			}
		}

		size_t total = 0;
		for (const auto& vec : local_next) total += vec.size();
		std::vector<int64_t> next_frontier;
		next_frontier.reserve(total);
		for (auto& vec : local_next) {
			next_frontier.insert(next_frontier.end(), vec.begin(), vec.end());
		}
		frontier.swap(next_frontier);
	}

	for (int64_t i = 0; i < N; ++i) {
		sp.distance_array[i] = dist[i].load(std::memory_order_relaxed);
		sp.parent_array[i] = parent[i].load(std::memory_order_relaxed);
	}
	delete[] dist;
	delete[] parent;
	delete[] seen_epoch;

	auto end = std::chrono::high_resolution_clock::now();
	sp.time_ms = std::chrono::duration<double, std::milli>(end - start).count();
	sp.teps = (double)g.length / (sp.time_ms * 1e-3);

	return sp;
}

shortest_path sssp_parallel(graph g, int64_t root) {
	auto start = std::chrono::high_resolution_clock::now();

	const float INF = std::numeric_limits<float>::infinity();
	int64_t N = g.nb_nodes;

	shortest_path sp;
	sp.distance_array = (float*)malloc(N * sizeof(float));
	sp.parent_array = (int64_t*)malloc(N * sizeof(int64_t));

	float* dist_old = (float*)malloc(N * sizeof(float));
	float* dist_new = (float*)malloc(N * sizeof(float));

	for (int64_t i = 0; i < N; i++) {
		dist_old[i] = INF;
		dist_new[i] = INF;
		sp.parent_array[i] = -1;
	}

	dist_old[root] = 0.0f;
	sp.parent_array[root] = root;

	bool changed = true;

	for (int64_t iter = 0; iter < N - 1 && changed; iter++) {
		changed = false;

		// copie des distances
		std::memcpy(dist_new, dist_old, N * sizeof(float));

#pragma omp parallel for schedule(static)
		for (int64_t u = 0; u < N; u++) {
			float du = dist_old[u];
			if (du == INF)
				continue;

			for (int64_t i = g.slicing_idx[u]; i < g.slicing_idx[u + 1]; i++) {
				int64_t v = g.neighbors[i];
				float w = g.weights[i];
				float alt = du + w;

				if (alt < dist_new[v]) {
#pragma omp critical
					{
						if (alt < dist_new[v]) {
							dist_new[v] = alt;
							sp.parent_array[v] = u;
							changed = true;
						}
					}
				}
			}
		}

		// swap buffers
		std::swap(dist_old, dist_new);
	}

	// résultat final
	std::memcpy(sp.distance_array, dist_old, N * sizeof(float));

	free(dist_old);
	free(dist_new);

	auto end = std::chrono::high_resolution_clock::now();
	sp.time_ms = std::chrono::duration<double, std::milli>(end - start).count();

	sp.teps = (double)g.length / (sp.time_ms * 1e-3);

	return sp;
}