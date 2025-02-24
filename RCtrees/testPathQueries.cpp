#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <cmath>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../include/parlay/internal/get_time.h"
#include "../examples/samplesort.h"
#include "RCdynamic.h"
#include "utils.h"
#include "random_trees.h"
#include "path_query.h"
#include "treeGen.h"
#include "ternarizer.h"

// **************************************************************
// Driver
// **************************************************************
using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

using datatype = double;

const vertex max_degree = 3;

parlay::sequence<std::pair<vertex,vertex>> generate_random_pairs(const parlay::sequence<vertex>& parents, const parlay::sequence<vertex>& map, long num_gens = 1)
{
    parlay::random_generator gen(0);
    std::uniform_int_distribution<long> dis(0, parents.size()-1);


    return std::move(parlay::filter(parlay::delayed_tabulate(num_gens, [&] (long i) {
        auto r = gen[i];
        long random_value = dis(r);
        long child = random_value;
        long parent = map[random_value];

        return std::pair<long, long>(child, parent);
    }), [] (auto pr){return std::get<0>(pr) != std::get<1>(pr);}));

}

int main(int argc, char* argv[]) {
    auto usage = "Usage: testMST.out [--graph-size <graph-size>] [-n <graph-size>] [--num-queries <number of path queries to perform>]\nPrints compressed tree creation time, MST time, insertion time";

    srand(time(NULL));

    const vertex max_rand_size =10000000l; //100000000l;
    // vertex graph_size = rand() % max_rand_size; 
    vertex graph_size = 1; 
    vertex num_queries = 1;
    for(unsigned short i = 0; i < 3; i++) // random but leaning more towards higher values
    {
        vertex tg = rand() % max_rand_size;
        if (tg > graph_size)
            graph_size = tg;
    }
    
   

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph-size" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "--num-queries" && i + 1 < argc) {
            num_queries = std::stol(argv[++i]);
        } else if (i == 1 && arg.find_first_not_of("0123456789") == std::string::npos) {
            // Handle the case where a single number is provided as graph_size
            graph_size = std::stol(arg);
        } else {
            std::cout << usage << std::endl;
            return -1;
        }
    }

    if (graph_size == 0) {
        std::cout << "graph_size should be an integer greater than 1" << std::endl;
        return -1;
    }

    // std::cout << "Working with a graph of size " << graph_size << std::endl; 

    TreeGen<long, double> TG(graph_size, 0.0f, 1.0f, 0.1, 20.0, exponential, true, 0);

    TG.generateInitialEdges();

    auto retedges = TG.getAllEdges();

    ternarizer<long, double> TR(graph_size, 0.0f);

    auto ternerized_initial_edges = TR.add_edges(retedges);
    
    parlay::sequence<cluster<vertex, datatype>> clusters;

    // // Measure creation time
    // auto start_creation = std::chrono::high_resolution_clock::now();
    create_base_clusters(clusters, ternerized_initial_edges, static_cast<long>(3), graph_size*extra_tern_node_factor);
    create_RC_tree(clusters, graph_size*extra_tern_node_factor, static_cast<double>(0.0f), [] (double A, double B) {return A+B;}, false);

    auto random_pairs = generate_random_pairs(TG.parents, TG.random_perm_map, num_queries);

    // std::cout << "Adding " << random_add_edges.size() << " new edges " << std::endl;
    
    parlay::internal::timer t1;

    t1.start();

    auto results = parlay::tabulate(random_pairs.size(), [&] (long i) {
        auto& rp = random_pairs[i]; 
        const double identity = 0.0f;
        return pathQuery(rp.first, rp.second, clusters, identity, [] (double A, double B) {return A+B;});
    });

    std::cout << t1.next_time() << std::endl;
    
     
    if(graph_size <= 100)
        printTree(clusters);

    // if(graph_size <= 100)
    //     printTree(clusters);

    deleteRCtree(clusters);


    return 0;
}
