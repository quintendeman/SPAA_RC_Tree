#include <iostream>
#include <parlay/primitives.h>
#include "treeGen.h"
#include <parlay/hash_table.h>
#include <chrono>

using wedge = std::tuple<long,long,double>;

struct hash_tuplic {
  using kType = std::tuple<long, long, double>; // Key type: the 3-element tuple
  using eType = void*; // Value type: pointer to the tuple

  // Define an empty value as nullptr
  eType empty() { return nullptr; }

  // Extract the key (tuple) from the value (pointer)
  kType getKey(eType v) {
    return *static_cast<wedge*>(v); // Dereference pointer to get the tuple
  }

    // Hash function for the tuple
    size_t hash(kType v) {
        // Hash each element of the tuple and combine
        static const long random_prime1 = 1299059;
        static const long random_prime2 = 1149769;
        static const long random_prime3 = 602057;

        size_t h1 = std::hash<long>()(std::get<0>(v)) * random_prime1;
        size_t h2 = std::hash<long>()(std::get<1>(v)) * random_prime2;
        size_t h3 = std::hash<double>()(std::get<2>(v)) * random_prime3;
        return h1 ^ h2; // Combine hashes using XOR and shifts
    }

  // Comparison function for tuples
  int cmp(kType v1, kType v2) {
    // Lexicographically compare the tuples
    if (v1 > v2) return 1;
    if (v1 == v2) return 0;
    return -1;
  }

  // Replace query: always return false (no replacement logic here)
  bool replaceQ(eType, eType) { return false; }

  // Update: return the current value
  eType update(eType v, eType) { return v; }

  // Atomic compare-and-swap (CAS) operation
  bool cas(eType* p, eType o, eType n) {
    return std::atomic_compare_exchange_strong_explicit(
      reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
  }
};


struct TupleFirstHash {
    size_t operator()(const long& t) const {
        return t + 1; // Hash the first element
    }
};

// Custom equality function for std::tuple<long, long, double>
struct TupleFirstEqual {
    bool operator()(const long& t1, const long& t2) const {
        return t1 == t2; // Compare the first element
    }
};


int main()
{
    const long graph_size_bits = 25;
    const long graph_size = 1l << graph_size_bits;

    const double min_weight = 0.0;
    const double max_weight = 100.0f;

    TreeGen<long, double> TG(graph_size, min_weight, max_weight, 0.00001, 2, constant);

    TG.generateInitialEdges();

    auto retedges = TG.getAllEdges();

    auto ret_edges_ = parlay::delayed_tabulate(retedges.size(), [&] (long i){
        auto retedge = retedges[i];
        std::swap(std::get<0>(retedge), std::get<1>(retedge));
        return retedge;
    });

    auto output_groups = group_by_key(ret_edges_, TupleFirstHash(), TupleFirstEqual()); // by default hashes and equalities according to first value and I don't know how to fix
    // slow to swap but still same order of complexity so whatever

    // for(auto& og : output_groups)
    // {
    //     std::cout << og.first << ": ";
    //     for(auto& val : og.second)
    //         std::cout << val << " ";
    //     std::cout << std::endl;
    // }

    using type_for_hash = wedge;

    auto HT = parlay::hashtable<hash_tuplic> (graph_size, hash_tuplic(), 5.0f);
    
    auto HT_entries = new wedge[graph_size];
    parlay::parallel_for(0, ret_edges_.size(), [&] (long i){
        HT_entries[i] = ret_edges_[i];
    });

    std::cout << "Doing " << ret_edges_.size() << " inserts takes... ";

    auto start = std::chrono::high_resolution_clock::now();

    parlay::parallel_for(0, ret_edges_.size(), [&] (long i ){
        HT.insert(&HT_entries[i]);
    });

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Print the duration
    std::cout << duration.count() << " ms" << std::endl;



    
    // for(long i = 0; i < ret_edges_.size(); i++)
    // {
    //     auto ptr_to_ret_edge = HT.find(ret_edges_[i]);
    //     wedge W = *static_cast<wedge*>(ptr_to_ret_edge);
    //     std::cout << std::get<0>(W) << " " << std::get<1>(W) << " " << std::get<2>(W) << std::endl;
    // }




    delete[] HT_entries;


}