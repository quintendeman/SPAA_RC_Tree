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

        create_base_clusters(clusters, ret_edge_modified.second, max_degree, graph_size * extra_tern_node_factor);

        create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);

        auto static_creation_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double>  dur1 = static_creation_middle - static_creation_start;
        std::chrono::duration<double>  dur2 = static_creation_end - static_creation_start;

        // std::cout << "static elapsed time: " << dur2.count() << " seconds of which " << dur1.count() << " were spent ternerizing \n";

        // auto del_edges = generate_delete_edges(graph_size);

        // auto retpair = TR.delete_edges(del_edges);

        // batchInsertEdge(retpair.second, retpair.first, clusters, (double) 0.0f, [] (double A, double B) {return A+B;});

        testPathQueryValid(clusters, TG.parents, TG.weights, TG.random_perm_map, graph_size);

        dynamic_gen_time = 0.0f;
        dynamic_counter = 0;
        dynamic_tern_time = 0.0f;

        const int maxi = 3;
        long total_edge_size = retedges.size();
        for (unsigned int i = 0; i < maxi; i++) 
        {
          
          double prob = (static_cast<double>(rand() % 10000)) / 10000.0f;
          // std::cout << "i: " << i << std::endl;

          double increment = (static_cast<double>(2 * i))/maxi;
          prob+=increment;
          if(prob > 1.0f)
            prob = 1.0f;

          // Skip measuring this part
          auto del_pairs = TG.generateDeleteEdges(prob); // Don't measure this

          // std::cout << "Deleting " << del_pairs.size() << " edges " << std::endl;

          // Measure TR.delete_edges(del_pairs)
          auto delete_start = std::chrono::high_resolution_clock::now();
          auto ret_seqs = TR.delete_edges(del_pairs);
          auto delete_end = std::chrono::high_resolution_clock::now();

          // Skip measuring this part

          // Measure batchInsertEdge(ret_seqs.second, ret_seqs.first, clusters, (double) 0.0f, ...);
          
          prob = (static_cast<double>(rand() % 10000)) / 10000.0f;
          prob+=increment;
          if(prob > 1.0f)
            prob = 1.0f;

          // Skip measuring this part
          auto add_truples = TG.generateAddEdges(prob); // Don't measure this
          
          // std::cout << "Adding " << add_truples.size() << " edges" << std::endl;

          // std::cout << "Total edges modified -" << del_pairs.size() << "+" << add_truples.size() << ":" << del_pairs.size() + add_truples.size() << std::endl;

          // Measure TR.add_edges(add_truples)
          auto add_start = std::chrono::high_resolution_clock::now();
          add_truples = TR.add_edges(add_truples);
          add_truples.append(ret_seqs.first);
          auto add_end = std::chrono::high_resolution_clock::now();

          // Measure batchInsertEdge(empty_pairs, add_truples, clusters, (double) 0.0f, ...);
          auto final_insert_start = std::chrono::high_resolution_clock::now();
          parlay::sequence<std::pair<long, long>> empty_pairs;
          // batchInsertEdge(ret_seqs.second, ret_seqs.first, clusters, (double) 0.0f, [] (double A, double B) {return A+B;});
          batchInsertEdge(ret_seqs.second, add_truples, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
          auto final_insert_end = std::chrono::high_resolution_clock::now();

          testPathQueryValid(clusters, TG.parents, TG.weights, TG.random_perm_map, graph_size); // Don't measure this
          
          // Calculate and print elapsed times for the relevant operations
          std::chrono::duration<double> delete_elapsed = delete_end - delete_start;
          std::chrono::duration<double> add_elapsed = add_end - add_start;
          std::chrono::duration<double> final_insert_elapsed = final_insert_end - final_insert_start;

          // Print the time taken for each measured section
          // std::cout << "  Time spent ternerizing: " << delete_elapsed.count() + add_elapsed.count() << " seconds" << std::endl;
          // std::cout << "  batchInsertEdge: " << final_insert_elapsed.count() << " seconds" << std::endl;
          // std::cout << "  TR.add_edges: " << add_elapsed.count() << " seconds" << std::endl;
          // std::cout << "  batchInsertEdge (2nd call): " << final_insert_elapsed.count() << " seconds" << std::endl;
          if(i > maxi/2)
          {
            dynamic_gen_time+= delete_elapsed.count() + add_elapsed.count() + final_insert_elapsed.count();
            dynamic_tern_time+=delete_elapsed.count() + add_elapsed.count();
            dynamic_counter++;
          }
      }
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