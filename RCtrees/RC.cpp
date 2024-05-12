/*
  This code generates a tree and then creates an RC tree out of this tree.
*/


#include <iostream>
#include <string>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/include/parlay/internal/get_time.h"
#include "RC.h"
#include "/scratch/parlaylib/examples/samplesort.h"
#include <chrono>
#include <iomanip>
#include <fstream>

#include <cmath>



// **************************************************************
// Driver
// **************************************************************
using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

const vertex max_degree = 8;


int main(int argc, char* argv[]) {
    auto usage = "Usage: RC <n>";
   
    if (argc != 2) {std::cout << usage << std::endl;
        return -1;
    }

    vertex n = 0;
    try { n = std::stoi(argv[1]); }
    catch (...) {}
    if (n == 0) {
        std::cout << "n should be an integer greater than 1" << std::endl;
        return -1;
    }

    auto old_n = n;

    n = pow(2, floor(log2(n)));

    // std::cout << "Setting n to closest (lower) power of 2, so "  << old_n << " => " << n <<  std::endl;

    auto parents = generate_tree_graph(n);

    degree_cap_parents(parents, max_degree);

    graph G;

    G = convert_parents_to_graph(G, parents);

    parlay::sequence<std::tuple<vertex, vertex, vertex>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, vertex> (i, parents[i], i+123210*i % 12312);
    });
    weighted_edges = parlay::filter(weighted_edges, [&] (std::tuple<vertex, vertex, vertex> wedge) {
        return std::get<0>(wedge) != std::get<1>(wedge);
    });

    parlay::sequence<cluster<vertex> > clusters;

    create_base_clusters(G, clusters, max_degree);

    
    // auto start = std::chrono::high_resolution_clock::now();

    create_RC_tree(clusters, n);

    adjust_weights(clusters, weighted_edges);


    delete_RC_Tree_edges(clusters);

    // auto end = std::chrono::high_resolution_clock::now();

    // std::chrono::duration<double> duration = end - start;

    // std::cout << std::fixed << std::setprecision(9) << n << "," << duration.count() << std::endl;
    
    // // This writing code is borrowed from chatgpt
    // std::ofstream file("data.csv", std::ios_base::app);
    // if (file.is_open()) {
    //     file << std::fixed << std::setprecision(9) << n << "," << duration.count() << std::endl;
    //     file.close();
    // } else {
    //     std::cerr << "Error opening file!" << std::endl;
    // }


    return 0;
}
