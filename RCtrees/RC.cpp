/*
  This code generates a tree and then creates an RC tree out of this tree.
  had help from chatgpt for the cli, i am lazy
*/

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

const vertex max_degree = 8;

int main(int argc, char* argv[]) {
    auto usage = "Usage: RC [--graph-size=<graph-size>] [--num-queries=<num-queries>] [--print-creation] [--print-query]";

    vertex graph_size = 100; // Default value
    vertex num_queries = 100; // Default value
    bool print_creation = false;
    bool print_query = false;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--graph-size=") == 0) {
            graph_size = std::stol(arg.substr(13));
        } else if (arg.find("--num-queries=") == 0) {
            num_queries = std::stol(arg.substr(14));
        } else if (arg == "--print-creation") {
            print_creation = true;
        } else if (arg == "--print-query") {
            print_query = true;
        } else {
            std::cout << usage << std::endl;
            return -1;
        }
    }

    
    if (graph_size == 0) {
        std::cout << "graph_size should be an integer greater than 1" << std::endl;
        return -1;
    }

    // auto old_graph_size = graph_size;
    // graph_size = pow(2, floor(log2(graph_size)));

    if(num_queries > graph_size)
        num_queries = graph_size;


    if (graph_size < 16)
        graph_size = 16;

    // std::cout << "Setting graph_size to closest (lower) power of 2 greater than 8, so "  << old_graph_size << " => " << graph_size <<  std::endl;

    auto parents = generate_tree_graph(graph_size);
    degree_cap_parents(parents, max_degree);
    graph G;
    G = convert_parents_to_graph(G, parents);

    parlay::sequence<std::tuple<vertex, vertex, vertex>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, vertex>(i, parents[i], i+123210*i % 12312);
    });
    weighted_edges = parlay::filter(weighted_edges, [&] (std::tuple<vertex, vertex, vertex> wedge) {
        return std::get<0>(wedge) != std::get<1>(wedge);
    });

    parlay::sequence<cluster<vertex> > clusters;

    // Measure creation time
    auto start_creation = std::chrono::high_resolution_clock::now();
    create_base_clusters(G, clusters, max_degree);
    create_RC_tree(clusters, graph_size);
    adjust_weights(clusters, weighted_edges, [] (vertex a, vertex b) {
        return a < b ? b : a;
    });
    auto end_creation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> creation_time = end_creation - start_creation;

    // Measure query time
    auto start_query = std::chrono::high_resolution_clock::now();
    parlay::parallel_for(0, num_queries, [&] (vertex i) {
        queryPath((vertex)((i + num_queries / 10) % num_queries), i, (vertex)-1, clusters, [] (vertex a, vertex b) {
            return a < b ? b : a;
        });
    });
    auto end_query = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> query_time = end_query - start_query;

    delete_RC_Tree_edges(clusters);

    if (print_creation) {
        // std::cout << "Creation time: " << std::setprecision(6) << creation_time.count() << " seconds" << std::endl;
        std::cout << graph_size << "," << std::setprecision(6) << creation_time.count() << std::endl;
    }
    if (print_query) {
        // std::cout << "Query time: " << std::setprecision(6) << query_time.count() << " seconds" << std::endl;
        std::cout << num_queries << "," << std::setprecision(6) << query_time.count() << std::endl;
    }
    if (!print_creation && !print_query) {
        std::cout << "Creation time: " << std::setprecision(6) << creation_time.count() << " seconds" << std::endl;
        std::cout << "Query time: " << std::setprecision(6) << query_time.count() << " seconds" << std::endl;
    }

    return 0;
}
