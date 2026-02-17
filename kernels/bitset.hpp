#ifndef KERNELS_BITSET_HPP
#define KERNELS_BITSET_HPP

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <unordered_set>

// Simple bitset backed by contiguous 64-bit words.
// - Call bs_init to allocate/initialize
// - Call bs_free to release
// - bs_set / bs_clear / bs_test operate on single bits
// - bs_reset clears all bits
// - bs_any checks whether any bit is set
// - bs_for_each_set iterates over set bits, using CTZ for speed
struct BitSet {
	uint64_t* data = nullptr;
	size_t words = 0;
	int64_t nb_bits = 0;

	BitSet(int64_t nb_bits) {
		this->nb_bits = nb_bits;
		this->words = (nb_bits + 63) / 64;
		this->data = (uint64_t*)calloc(this->words, sizeof(uint64_t));
	}

	BitSet(const std::unordered_set<int64_t>& indices, int64_t nb_bits) : BitSet(nb_bits) {
		for (int64_t idx : indices) {
			this->insert(idx);
		}
	}

	~BitSet() {
		if (this->data != nullptr) {
			free(this->data);
		}
	}

	void insert(int64_t idx) {
		this->data[(size_t)idx >> 6] |= (1ULL << (idx & 63));
	}

	void remove(int64_t idx) {
		this->data[(size_t)idx >> 6] &= ~(1ULL << (idx & 63));
	}

	bool contains(int64_t idx) const {
		return ((this->data[(size_t)idx >> 6] >> (idx & 63)) & 1ULL) != 0;
	}

	void clear() {
		if (this->data) {
			memset(this->data, 0, this->words * sizeof(uint64_t));
		}
	}

	bool empty() const {
		for (size_t i = 0; i < this->words; ++i) {
			if (this->data[i] != 0U) {
				return false;
			}
		}
		return true;
	}

	template <typename Func> void for_each(Func&& func) const {
		for (size_t w = 0; w < this->words; ++w) {
			uint64_t word = this->data[w];
			while (word) {
				const int tz = __builtin_ctzll(word);
				int64_t idx = static_cast<int64_t>(w) * 64 + tz;
				if (idx >= this->nb_bits) {
					break;
				}
				func(idx);
				word &= word - 1;
			}
		}
	}
};

void bitset_swap(BitSet& a, BitSet& b) {
	std::swap(a.data, b.data);
	std::swap(a.words, b.words);
	std::swap(a.nb_bits, b.nb_bits);
}

// Iterate over set bits in index order. Callback receives (int64_t idx).
// Uses a fast loop relying on __builtin_ctzll to locate least-significant set bit.
// Macro iterator over set bits. Usage:
// bs_each_set(bs, idx) expands to a loop where `idx` is an `int64_t` variable
// holding the current set-bit index. Example:
// bs_each_set(bs, i) { do_work(i); if (i==target) break; }
// #define bs_for_each(bs, idx) \
// 	for (size_t _bs_w = 0; _bs_w < (bs).words; _bs_w++)                                                                                    \
// 		for (uint64_t _bs_word = (bs).data[_bs_w]; _bs_word; _bs_word &= _bs_word - 1)                                                     \
// 			for (int64_t idx = (int64_t)_bs_w * 64 + (int64_t)__builtin_ctzll(_bs_word); (idx) < (bs).nb_bits; (idx) = (int64_t)-1)

// Class-based atomic bitset backed by contiguous 64-bit atomic words.
// - Construct to allocate/initialize
// - insert/remove/contains operate atomically
// - try_insert returns true when this call changed the bit from 0->1
// - clear resets all words to zero (atomic stores)
// - empty checks whether any bit is set (atomic loads)
// - for_each_set iterates over a snapshot of set bits and calls a provided callable
struct AtomicBitSet {
	std::atomic<uint64_t>* data;
	size_t words;
	int64_t nb_bits;

	AtomicBitSet(int64_t nb_bits) : data(nullptr), words(0), nb_bits(nb_bits) {
		assert(nb_bits >= 0);
		words = (nb_bits + 63) / 64;
		data = new std::atomic<uint64_t>[words];
		for (size_t i = 0; i < words; i++) {
			data[i].store(0ULL, std::memory_order_relaxed);
		}
	}

	AtomicBitSet(const std::unordered_set<int64_t>& indices, int64_t nb_bits) : AtomicBitSet(nb_bits) {
		for (int64_t idx : indices) {
			insert(idx);
		}
	}

	~AtomicBitSet() {
		delete[] data;
		data = nullptr;
	}

	void insert(int64_t idx) {
		assert(0 <= idx && idx < nb_bits);
		const size_t w = static_cast<size_t>(idx) >> 6;
		const uint64_t mask = (1ULL << (idx & 63));
		data[w].fetch_or(mask, std::memory_order_acq_rel);
	}

	void remove(int64_t idx) {
		assert(0 <= idx && idx < nb_bits);
		const size_t w = static_cast<size_t>(idx) >> 6;
		const uint64_t mask = ~(1ULL << (idx & 63));
		data[w].fetch_and(mask, std::memory_order_acq_rel);
	}

	bool contains(int64_t idx) const {
		assert(0 <= idx && idx < nb_bits);
		const size_t w = static_cast<size_t>(idx) >> 6;
		const uint64_t mask = (1ULL << (idx & 63));
		return (data[w].load(std::memory_order_acquire) & mask) != 0ULL;
	}

	void clear() {
		for (size_t i = 0; i < words; ++i) {
			data[i].store(0ULL, std::memory_order_release);
		}
	}

	bool empty() const {
		for (size_t i = 0; i < words; ++i) {
			if (data[i].load(std::memory_order_acquire) != 0ULL)
				return false;
		}
		return true;
	}

	// Iterate over a snapshot of set bits and call `func(int64_t idx)`.
	// The callable `func` must accept one int64_t parameter.
	template <typename Func> void for_each(Func&& func) const {
		for (size_t w = 0; w < words; ++w) {
			uint64_t word = data[w].load(std::memory_order_acquire);
			while (word) {
				const int tz = __builtin_ctzll(word);
				int64_t idx = static_cast<int64_t>(w) * 64 + tz;
				if (idx >= nb_bits) {
					break;
				}
				func(idx);
				word &= word - 1;
			}
		}
	}

	template <typename Func> void parallel_for_each(Func&& func) const {
#pragma omp parallel for schedule(guided)
		for (size_t w = 0; w < words; ++w) {
			uint64_t word = data[w].load(std::memory_order_acquire);
			while (word) {
				const int tz = __builtin_ctzll(word);
				int64_t idx = static_cast<int64_t>(w) * 64 + tz;
				if (idx >= nb_bits) {
					break;
				}
				func(idx);
				word &= word - 1;
			}
		}
	}
};

void atomic_bitset_swap(AtomicBitSet& a, AtomicBitSet& b) {
	std::swap(a.data, b.data);
	std::swap(a.words, b.words);
	std::swap(a.nb_bits, b.nb_bits);
}

#endif // KERNELS_BITSET_HPP
