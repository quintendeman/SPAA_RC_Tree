#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../include/parlay/internal/get_time.h"
#include "RC.h"
#include "../examples/samplesort.h"

// **************************************************************
// Driver
// **************************************************************
using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

using datatype = float;

const vertex max_degree = 3;

static const vertex batch_insertion_size = 1000;

int main(int argc, char* argv[]) {
    auto usage = "Usage: RC [--graph-size <graph-size>] [-n <graph-size>] [--num-queries <num-queries>] [--print-creation] [--do-height <true|false>] [--randomized <true|false>]";

    vertex graph_size = 100; // Default value
    vertex num_queries = 100; // Default value
    bool print_creation = false;
    bool randomized = false; // Default value

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph-size" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "--num-queries" && i + 1 < argc) {
            num_queries = std::stol(argv[++i]);
        } else if (arg == "--print-creation") {
            print_creation = true;
        } else if (arg == "--randomized" && i + 1 < argc) {
            std::string value = argv[++i];
            randomized = (value == "true");
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

    if (num_queries > graph_size)
        num_queries = graph_size;

    if (graph_size < 10)
        graph_size = 10;

    std::cout << "Working with a graph of size " << graph_size << std::endl; 

    auto parents = generate_tree_graph(graph_size);
    degree_cap_parents(parents, max_degree);
    graph G;
    G = convert_parents_to_graph(G, parents);

    parlay::sequence<std::tuple<vertex, vertex, datatype>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, datatype>(i, parents[i], (datatype) (i + parents[i] * 10 ));
    });
    weighted_edges = parlay::filter(weighted_edges, [&] (std::tuple<vertex, vertex, datatype> wedge) {
        return std::get<0>(wedge) != std::get<1>(wedge);
    });

    parlay::sequence<cluster<vertex, datatype>> clusters;

    // Measure creation time
    auto start_creation = std::chrono::high_resolution_clock::now();
    create_base_clusters(G, clusters, max_degree);
    create_RC_tree(clusters, graph_size, randomized);
    set_heights(clusters);
    auto end_creation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> creation_time = end_creation - start_creation;

    std::cout << "There are " << weighted_edges.size() << " edges being inserted" << std::endl;

    batchModifyEdgeWeights(weighted_edges, [] (datatype a, datatype b) {
        return a + b;
    }, clusters);

    printTree(clusters);

    // const datatype defretval = 0.0f;

    // PathQuery(&clusters[0], &clusters[G.size()/2], 0.0f, [] (datatype a, datatype b) {
    //     return a + b;
    // });

    // parlay::random_generator gen;
    // std::uniform_int_distribution<vertex> dis(0, graph_size-1);

    // auto random_indices = parlay::tabulate(batch_insertion_size < graph_size ? batch_insertion_size : graph_size, [&] (vertex i) 
    // {
    //     auto r = gen[i];
    //     auto random_index = dis(r);
    //     return random_index;
    // });

    // auto delete_pairs = parlay::tabulate(random_indices.size(), [&] (vertex index)
    // {
    //     auto& i = random_indices[index];
    //     parents[i] = i;    
    //     if (i & 1)
    //         {
    //             return std::pair<vertex, vertex>(i, parents[i]);
    //         }

    //     return std::pair<vertex, vertex>(parents[i], i);        
    // });



    // auto add_edges = parlay::tabulate(random_indices.size(), [&] (vertex index)
    // {   
    //     auto& i = random_indices[index];
    //     vertex& v = i;
    //     auto r = gen[i];
    //     vertex w;
    //     if(v > 0)
    //         w = dis(r) % v; // something to the left of w
    //     else
    //         w = 0;
    //     datatype weight = (double) (dis(r) * dis(r));
        
    //     parents[i] = w;
    //     if (i & 1)
    //     {
    //         return std::make_tuple(v, w, weight);
    //     }
    //     return std::make_tuple(w, v, weight);
    // });
    
    // degree_cap_parents(parents, max_degree);

    // add_edges = parlay::filter(add_edges, [&] (auto edge_tuple) {
    //     const vertex& i = std::get<0>(edge_tuple);
    //     return i != parents[i];
    // });

    // batchInsertEdge(delete_pairs, add_edges, clusters);

    // printTree(clusters);

    deleteRCtree(clusters);

    if (print_creation) {
        std::cout << graph_size << "," << std::setprecision(6) << creation_time.count() << std::endl;
    } else {
        std::cout << "Creation time: " << std::setprecision(6) << creation_time.count() << " seconds" << std::endl;
    }

    return 0;
}
