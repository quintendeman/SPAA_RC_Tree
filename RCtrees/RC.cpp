#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
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

const vertex max_degree = 3;

int main(int argc, char* argv[]) {
    auto usage = "Usage: RC [--graph-size <graph-size>] [-n <graph-size>] [--num-queries <num-queries>] [--print-creation] [--do-height <true|false>] [--randomized <true|false>]";

    vertex graph_size = 100; // Default value
    vertex num_queries = 100; // Default value
    bool print_creation = false;
    bool do_height = true; // Default value
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
        } else if (arg == "--do-height" && i + 1 < argc) {
            std::string value = argv[++i];
            do_height = (value == "true");
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

    if (graph_size < 16)
        graph_size = 16;

    auto parents = generate_tree_graph(graph_size);
    degree_cap_parents(parents, max_degree);
    graph G;
    G = convert_parents_to_graph(G, parents);

    parlay::sequence<std::tuple<vertex, vertex, vertex>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, vertex>(i, parents[i], i + 123210 * i % 12312);
    });
    weighted_edges = parlay::filter(weighted_edges, [&] (std::tuple<vertex, vertex, vertex> wedge) {
        return std::get<0>(wedge) != std::get<1>(wedge);
    });

    parlay::sequence<cluster<vertex>> clusters;

    // Measure creation time
    auto start_creation = std::chrono::high_resolution_clock::now();
    create_base_clusters(G, clusters, max_degree);
    create_RC_tree(clusters, graph_size, do_height, randomized);
    auto end_creation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> creation_time = end_creation - start_creation;

    // printTree(clusters);

    deleteRCtree(clusters);

    if (print_creation) {
        std::cout << graph_size << "," << std::setprecision(6) << creation_time.count() << std::endl;
    } else {
        std::cout << "Creation time: " << std::setprecision(6) << creation_time.count() << " seconds" << std::endl;
    }

    return 0;
}
