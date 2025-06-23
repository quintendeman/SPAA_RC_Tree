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
#include "RCdynamic.h"
#include "../examples/samplesort.h"
#include "incMST.h"
#include "utils.h"
#include "random_trees.h"
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


// This wipes out all existing weights, replacing them with random ones
// The lambda function is MAX
void insert_random_weights(const parlay::sequence<vertex>& parents, parlay::sequence<cluster<vertex, datatype>>& clusters)
{
    std::random_device rd;                         // Seed generator
    std::mt19937 gen(rd());                        // Standard Mersenne Twister engine seeded with rd()
    std::uniform_int_distribution<vertex> dist(0, clusters.size() * 100);   

    parlay::sequence<std::tuple<vertex, vertex, datatype>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, datatype>(i, parents[i], (datatype) ( dist(gen)));
    }); 
    
    weighted_edges = parlay::filter(weighted_edges, [&] (auto edge) {
        return std::get<0>(edge) != std::get<1>(edge);
    });

    batchModifyEdgeWeights(weighted_edges, [] (datatype a, datatype b) {
        return a > b ? a : b; // MAX (not min)
    }, clusters);

    return;
}


parlay::sequence<std::tuple<vertex,vertex,double>> generate_random_edges(const parlay::sequence<vertex>& parents, const parlay::sequence<vertex>& map, long num_gens = 1)
{
    parlay::random_generator gen(0);
    std::uniform_int_distribution<long> dis(0, parents.size()-1);
    std::uniform_real_distribution<double> dis_ur(0.0f, 1.0f);

    return std::move(parlay::filter(parlay::delayed_tabulate(num_gens, [&] (long i) {
        auto r = gen[i];
        long random_value = dis(r);
        long child = map[random_value];
        long parent = map[parents[random_value]];
        double random_weight = dis_ur(r);
        return std::tuple<long, long, double>(child, parent, random_weight);
    }), [] (auto tpl){return std::get<0>(tpl) != std::get<1>(tpl);}));

}

int main(int argc, char* argv[]) {
    auto usage = "Usage: testMST.out [--graph-size <graph-size>] [-n <graph-size>] [--num-additions <extra edges added>]\nPrints compressed tree creation time, MST time, insertion time";

    srand(time(NULL));

    const vertex max_rand_size =10000000l; //100000000l;
    // vertex graph_size = rand() % max_rand_size; 
    vertex graph_size = 1; 
    vertex num_additions = 1;
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
        } else if (arg == "--num-additions" && i + 1 < argc) {
            num_additions = std::stol(argv[++i]);
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
    create_base_clusters(clusters, ternerized_initial_edges.second, static_cast<long>(3), graph_size*extra_tern_node_factor);
    create_RC_tree(clusters, graph_size*extra_tern_node_factor, static_cast<double>(0.0f), [] (double A, double B) {return A+B;}, false);

    auto random_add_edges = generate_random_edges(TG.parents, TG.random_perm_map, num_additions);

    // std::cout << "Adding " << random_add_edges.size() << " new edges " << std::endl;
    
    incrementMST(clusters, random_add_edges, TR);

     
    if(graph_size <= 100)
        printTree(clusters);

    // if(graph_size <= 100)
    //     printTree(clusters);

    deleteRCtree(clusters);


    return 0;
}
