#include <iostream>
#include <parlay/primitives.h>
#include "treeGen.h"
#include <parlay/hash_table.h>
#include <chrono>
#include <cstdlib>
#include "ternarizer.h"
#include "RC.h"
#include "RCdynamic.h"
#include "random_trees.h"
#include "path_query.h"

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
        size_t h = static_cast<size_t>(t);
        h ^= (h >> 30);  // First mix
        h *= 0xbf58476d1ce4e5b9ULL; // A simple prime multiplier
        h ^= (h >> 27);  // Second mix
        h *= 0x94d049bb133111ebULL; // Another prime multiplier
        h ^= (h >> 31);  // Final mix
        return h;
    }
};

// Custom equality function for std::tuple<long, long, double>
struct TupleFirstEqual {
    bool operator()(const long& t1, const long& t2) const {
        return t1 == t2; // Compare the first element
    }
};

parlay::sequence<wedge> generate_high_degree_graphs(long graph_size)
{
  auto retedges =  parlay::tabulate(graph_size, [graph_size] (long i) {
    long middle = graph_size/2;

    return std::tuple<long,long,double>(i, middle, (double) ((i * 12312 + i) % 12132));

  });

  return parlay::filter(retedges, [] (wedge we) {
      return std::get<0>(we) != std::get<1>(we);
  });
}

auto generate_delete_edges(long graph_size)
{
  return parlay::tabulate(graph_size/4, [graph_size] (long i) {
    return std::pair<long,long>(i, graph_size/2);
  });
}

static void print_help() {
    std::cout << "Usage: ./program [options]\n"
              << "Options:\n"
              << "  --graph_size <num>    Set graph size (default: 100000)\n"
              << "  --dist <e|g|c|u>      Set distribution: exponential (e), geometric (g), constant (c), uniform (u) (default: e)\n"
              << "  --ln <value>          Set ln parameter (default: 0.1)\n"
              << "  --mean <value>        Set mean value (default: 20.0)\n"
              << "  --randomized          Enable randomized mode\n"
              << "  --help                Show this help message\n";
}

int main(int argc, char* argv[])
{
    
    long graph_size = 100000;
    auto distribution = exponential; // either exponential, geometric, constant or linear
    double ln = 0.1;
    double mean = 20.0;
    bool randomized = false;

    // std::cout << "Argc " << argc << std::endl;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph_size" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "--dist" && i + 1 < argc) {
            std::string dist_arg = argv[++i];
            if (dist_arg == "e") {
                distribution = exponential;
            } else if (dist_arg == "g") {
                distribution = geometric;
            } else if (dist_arg == "c") {
                distribution = constant;
            } else if (dist_arg == "u") {
                distribution = uniform;
            } else {
                std::cerr << "Invalid distribution type\n";
                return 1;
            }
        } else if (arg == "--randomized") {
          randomized = true;
        } else if (arg == "--ln" && i + 1 < argc) {
            ln = std::stod(argv[++i]);
        } else if (arg == "--mean" && i + 1 < argc) {
            mean = std::stod(argv[++i]);
        }
    }
        

    const double min_weight = 0.0;
    const double max_weight = 100.0f;

    // std::cout << "Graph size " << graph_size << std::endl;

    double static_gen_time = 0.0f;
    double dynamic_gen_time = 0.0f;
    double static_tern_time = 0.0f;
    double dynamic_tern_time = 0.0f;
    int static_counter = 0;
    int dynamic_counter = 0;

    const int maxII = 8;

    for(int II = 0; II < maxII; II++)
    {
        // std::cout << "II " << II << std::endl;
        TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true, II);

        TG.generateInitialEdges();

        auto retedges = TG.getAllEdges();

        // auto retedges = generate_high_degree_graphs(graph_size);
        
        auto static_creation_start = std::chrono::high_resolution_clock::now(); //

        ternarizer<long, double> TR(graph_size, 0);

        auto ret_edge_modified = TR.add_edges(retedges);

        auto static_creation_middle = std::chrono::high_resolution_clock::now();
        
        // TR.print_state();
        // TR.verify_simple_tree();

        const long max_degree = 3;
        parlay::sequence<cluster<long, double>> clusters; 
        double defretval = 0.0;

        parlay::sequence<wedge> empty_edges_sequence;

        create_base_clusters(clusters, ret_edge_modified, max_degree, graph_size * extra_tern_node_factor);

        create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);

        auto static_creation_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double>  dur1 = static_creation_middle - static_creation_start;
        std::chrono::duration<double>  dur2 = static_creation_end - static_creation_start;

        
        deleteRCtree(clusters); 

      if(II > maxII/2)
      {
        static_gen_time+=dur2.count();
        static_tern_time+=dur1.count();
        static_counter++;
      }
    }

    // std::cout << "STATIC TIME " << static_gen_time/static_counter << " of which " << static_tern_time/static_counter << " was spent ternerizing" << std::endl;
    // std::cout << "DYNAMIC TIME " << dynamic_gen_time/dynamic_counter << " of which " << dynamic_tern_time/dynamic_counter << " was spent ternerizing" << std::endl;
    std::cout << parlay::internal::init_num_workers() << ",";
    std::cout << graph_size << "," << ln << "," << mean << ",";
    switch(distribution){
      case exponential:
        std::cout << "e";
        break;
      case geometric:
        std::cout << "g";
        break;
      case constant:
        std::cout << "c";
        break;
      case uniform:
        std::cout << "u";
        break;
    }

    
    
    std::cout << "," << static_gen_time/static_counter << "," << static_tern_time/static_counter << "," << dynamic_gen_time/dynamic_counter << "," << dynamic_tern_time/dynamic_counter << ",";
    if(randomized)
      std::cout << "true";
    else
      std::cout << "false";
    std::cout << std::endl;
}