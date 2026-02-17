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