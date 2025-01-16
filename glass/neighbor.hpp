#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <queue>
#include <vector>

#include "glass/memory.hpp"

namespace glass {

namespace searcher {

template <typename Block = uint64_t> struct Bitset {
  constexpr static int block_size = sizeof(Block) * 8;
  int nbytes;
  Block *data;
  explicit Bitset(int n)
      : nbytes((n + block_size - 1) / block_size * sizeof(Block)),
        data((uint64_t *)alloc64B(nbytes)) {
    memset(data, 0, nbytes);
  }
  ~Bitset() { free(data); }
  void set(int i) {
    data[i / block_size] |= (Block(1) << (i & (block_size - 1)));
  }
  bool get(int i) {
    return (data[i / block_size] >> (i & (block_size - 1))) & 1;
  }

  void *block_address(int i) { return data + i / block_size; }
};

template <typename dist_t = float> struct Neighbor {
  int id;
  dist_t distance;

  Neighbor() = default;
  Neighbor(int id, dist_t distance) : id(id), distance(distance) {}

  inline friend bool operator<(const Neighbor &lhs, const Neighbor &rhs) {
    return lhs.distance < rhs.distance ||
           (lhs.distance == rhs.distance && lhs.id < rhs.id);
  }
  inline friend bool operator>(const Neighbor &lhs, const Neighbor &rhs) {
    return !(lhs < rhs);
  }
};

template <typename dist_t> struct MaxHeap {
  explicit MaxHeap(int capacity) : capacity(capacity), pool(capacity) {}
  void push(int u, dist_t dist) {
    if (size < capacity) {
      pool[size] = {u, dist};
      std::push_heap(pool.begin(), pool.begin() + ++size);
    } else if (dist < pool[0].distance) {
      sift_down(0, u, dist);
    }
  }
  int pop() {
    std::pop_heap(pool.begin(), pool.begin() + size--);
    return pool[size].id;
  }
  void sift_down(int i, int u, dist_t dist) {
    pool[0] = {u, dist};
    for (; 2 * i + 1 < size;) {
      int j = i;
      int l = 2 * i + 1, r = 2 * i + 2;
      if (pool[l].distance > dist) {
        j = l;
      }
      if (r < size && pool[r].distance > std::max(pool[l].distance, dist)) {
        j = r;
      }
      if (i == j) {
        break;
      }
      pool[i] = pool[j];
      i = j;
    }
    pool[i] = {u, dist};
  }
  int size = 0, capacity;
  std::vector<Neighbor<dist_t>, align_alloc<Neighbor<dist_t>>> pool;
};

template <typename dist_t> struct MinMaxHeap {
  explicit MinMaxHeap(int capacity) : capacity(capacity), pool(capacity) {}
  bool push(int u, dist_t dist) {
    if (cur == capacity) {
      if (dist >= pool[0].distance) {
        return false;
      }
      if (pool[0].id >= 0) {
        size--;
      }
      std::pop_heap(pool.begin(), pool.begin() + cur--);
    }
    pool[cur] = {u, dist};
    std::push_heap(pool.begin(), pool.begin() + ++cur);
    size++;
    return true;
  }
  dist_t max() { return pool[0].distance; }
  void clear() { size = cur = 0; }

  int pop_min() {
    int i = cur - 1;
    for (; i >= 0 && pool[i].id == -1; --i)
      ;
    if (i == -1) {
      return -1;
    }
    int imin = i;
    dist_t vmin = pool[i].distance;
    for (; --i >= 0;) {
      if (pool[i].id != -1 && pool[i].distance < vmin) {
        vmin = pool[i].distance;
        imin = i;
      }
    }
    int ret = pool[imin].id;
    pool[imin].id = -1;
    --size;
    return ret;
  }

  int size = 0, cur = 0, capacity;
  std::vector<Neighbor<dist_t>, align_alloc<Neighbor<dist_t>>> pool;
};

template <typename dist_t> struct LinearPool {
  LinearPool(int n, int capacity, int = 0)
      : nb(n), capacity_(capacity), data_(capacity_ + 1), vis(n) {}

  int find_bsearch(dist_t dist) {
    int lo = 0, hi = size_;
    while (lo < hi) {
      int mid = (lo + hi) / 2;
      if (data_[mid].distance > dist) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  bool insert(int u, dist_t dist) {
    if (size_ == capacity_ && dist >= data_[size_ - 1].distance) {
      return false;
    }
    int lo = find_bsearch(dist);
    std::memmove(&data_[lo + 1], &data_[lo],
                 (size_ - lo) * sizeof(Neighbor<dist_t>));
    data_[lo] = {u, dist};
    if (size_ < capacity_) {
      size_++;
    }
    if (lo < cur_) {
      cur_ = lo;
    }
    return true;
  }

  int pop() {
    set_checked(data_[cur_].id);
    int pre = cur_;
    while (cur_ < size_ && is_checked(data_[cur_].id)) {
      cur_++;
    }
    return get_id(data_[pre].id);
  }

  bool has_next() const { return cur_ < size_; }
  int id(int i) const { return get_id(data_[i].id); }
  int size() const { return size_; }
  int capacity() const { return capacity_; }

  constexpr static int kMask = 2147483647;
  int get_id(int id) const { return id & kMask; }
  void set_checked(int &id) { id |= 1 << 31; }
  bool is_checked(int id) { return id >> 31 & 1; }

  int nb, size_ = 0, cur_ = 0, capacity_;
  std::vector<Neighbor<dist_t>, align_alloc<Neighbor<dist_t>>> data_;
  Bitset<uint64_t> vis;
};

template <typename dist_t> struct HeapPool {
  HeapPool(int n, int capacity, int topk)
      : nb(n), capacity_(capacity), candidates(capacity), retset(topk), vis(n) {
  }
  bool insert(int u, dist_t dist) {
    retset.push(u, dist);
    return candidates.push(u, dist);
  }
  int pop() { return candidates.pop_min(); }
  bool has_next() const { return candidates.size > 0; }
  int id(int i) const { return retset.pool[i].id; }
  int capacity() const { return capacity_; }
  int nb, size_ = 0, capacity_;
  MinMaxHeap<dist_t> candidates;
  MaxHeap<dist_t> retset;
  Bitset<uint64_t> vis;
};

template <typename dist_t> struct DynamicDiffPool {
  DynamicDiffPool(int n, int init_capacity, int reserved_capacity, float ipdiff, float init_max)
      : nb(n), capacity_(init_capacity), data_(reserved_capacity), 
      ipdiff_(ipdiff), init_max_(-init_max), init_size_(init_capacity), vis(n) {}

  void prepare_search() {
    threshold_ = std::min(init_max_, head_dist()) + ipdiff_;
  }

  int find_bsearch(dist_t dist) {
    int lo = 0, hi = size_;
    while (lo < hi) {
      int mid = (lo + hi) >> 1;
      if (data_[mid].distance > dist) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  int find_bsearch_dist(const dist_t dist, int hi) {
    int lo = 0;
    while (lo < hi) {
      int mid = (lo + hi) >> 1;
      if (dist < data_[mid].distance) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  bool insert(int u, dist_t dist) {
    if (size_ == capacity_ && dist >= data_[size_ - 1].distance) {
      return false;
    }
    int lo = find_bsearch(dist);
    std::memmove(&data_[lo + 1], &data_[lo],
                 (size_ - lo) * sizeof(Neighbor<dist_t>));
    data_[lo] = {u, dist};
    if (lo < cur_) {
      cur_ = lo;
    }

    if (size_ >= capacity_) {
      if (data_[static_cast<int>(capacity_ * 0.75)].distance <= threshold_) [[unlikely]] {  // first half all good
        capacity_ *= 2;
        if (data_.size() <= capacity_) [[unlikely]] {
          data_.resize(capacity_ + 1);
        }
        size_++;
      }
    }
    else {
      size_++;
    }

    if (lo == 0) [[unlikely]] {
      if (dist < init_max_) {
        threshold_ = dist + ipdiff_;
        shrink_flag_ = true;
      }
    }
    return true;
  }

  void clear() {
    size_ = 0;
    cur_ = 0;
  }

  void shrink() {
    if (!shrink_flag_) return;
    shrink_flag_ = false;
    while (capacity_ > _init_size &&
           (capacity_ / 2 >= size_ || data_[capacity_ / 2].distance > threshold_)) {
      capacity_ /= 2;
    }
    data_.resize(capacity_ + 1);
    size_ = std::min(size_, capacity_);
  }

  int pop() {
    set_checked(data_[cur_].id);
    int pre = cur_;
    while (cur_ < size_ && is_checked(data_[cur_].id)) {
      cur_++;
    }
    return get_id(data_[pre].id);
  }

  float head_dist() { return data_[0].distance; }
  float dist(size_t idx) { return data_[idx].distance; }
  int pick_size() {
    return find_bsearch_dist(threshold_, size_);
  }

  bool has_next() const { return cur_ < size_; }
  int id(int i) const { return get_id(data_[i].id); }
  int size() const { return size_; }
  int capacity() const { return capacity_; }

  constexpr static int kMask = 2147483647;
  int get_id(int id) const { return id & kMask; }
  void set_checked(int &id) { id |= 1 << 31; }
  bool is_checked(int id) { return id >> 31 & 1; }

  int nb, size_ = 0, cur_ = 0, capacity_;
  std::vector<Neighbor<dist_t>, align_alloc<Neighbor<dist_t>>> data_;
  Bitset<uint64_t> vis;
  float ipdiff_, init_max_, threshold_ = 0;
  bool shrink_flag_ = false;
  size_t init_size_;
};

} // namespace searcher

struct Neighbor {
  int id;
  float distance;
  bool flag;

  Neighbor() = default;
  Neighbor(int id, float distance, bool f)
      : id(id), distance(distance), flag(f) {}

  inline bool operator<(const Neighbor &other) const {
    return distance < other.distance;
  }
};

struct Node {
  int id;
  float distance;

  Node() = default;
  Node(int id, float distance) : id(id), distance(distance) {}

  inline bool operator<(const Node &other) const {
    return distance < other.distance;
  }
};

inline int insert_into_pool(Neighbor *addr, int K, Neighbor nn) {
  // find the location to insert
  int left = 0, right = K - 1;
  if (addr[left].distance > nn.distance) {
    memmove(&addr[left + 1], &addr[left], K * sizeof(Neighbor));
    addr[left] = nn;
    return left;
  }
  if (addr[right].distance < nn.distance) {
    addr[K] = nn;
    return K;
  }
  while (left < right - 1) {
    int mid = (left + right) / 2;
    if (addr[mid].distance > nn.distance) {
      right = mid;
    } else {
      left = mid;
    }
  }
  // check equal ID

  while (left > 0) {
    if (addr[left].distance < nn.distance) {
      break;
    }
    if (addr[left].id == nn.id) {
      return K + 1;
    }
    left--;
  }
  if (addr[left].id == nn.id || addr[right].id == nn.id) {
    return K + 1;
  }
  memmove(&addr[right + 1], &addr[right], (K - right) * sizeof(Neighbor));
  addr[right] = nn;
  return right;
}

} // namespace glass
