/*
  This code is a heavily modified version of the hotspots assignment
*/


#include <iostream>
#include <string>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/include/parlay/internal/get_time.h"
#include "RC.h"
#include "/scratch/parlaylib/examples/samplesort.h"

#include <cmath>

/*
    MAXIMUM DEGREE
    I will figure out how to make this dynamic later
*/


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

    std::cout << "Setting n to closest (lower) power of 2, so "  << old_n << " => " << n <<  std::endl;

    auto parents = generate_tree_graph(n);

    degree_cap_parents(parents, max_degree);

    graph G;

    G = convert_parents_to_graph(G, parents);

    utils::print_graph_stats(G);

    // degree_cap_graph(G, max_degree);
    parlay::sequence<cluster<vertex> > clusters;

    create_base_clusters(G, clusters);

    create_RC_tree(clusters, n);

    // for(uint i = 0; i < clusters.size(); i++)
    // {
    //     std::cout << i<< " "<<clusters[i].data.size() <<  " " << clusters[i].final_colour << " " << " ";
    //     if(clusters[i].state & live)
    //     {
    //         std::cout << "live ";
    //     }
    //     else if (clusters[i].state & nullary_cluster)
    //     {
    //         std::cout << "nullary ";
    //     }
    //     else if (clusters[i].state & binary_cluster)
    //     {
    //         std::cout << "binary ";
    //     }
    //     else if (clusters[i].state & unary_cluster)
    //     {
    //         std::cout << "unary ";
    //     }
    //     for(uint j = 0; j < clusters[i].data.size(); j++)
    //     {
    //         if(clusters[i].data[j] == NULL)
    //             std::cout << "null ";
    //         else
    //             std::cout << clusters[i].data[j]->index << " ";
    //     }
    //     if(clusters[i].is_MIS)
    //     {
    //         std::cout << "\u2713";
    //     }
    //     std::cout << std::endl;
        
    // }

    return 0;
}

// graph G = utils::rmat_symmetric_graph(n, 1*n);

    // degree_cap_graph(G, max_degree);

    // utils::print_graph_stats(G);

    // auto colours = colour_chains_to_logn(G, max_degree);

    // auto vertices = parlay::tabulate(n, [&] (vertex v) {return v;});

    // is_valid_colouring(G, colours); // Takes a while but dw about it

    // auto result = vertices;

    // parlay::sequence<unsigned long> offsets = counting_sort(vertices.begin(), vertices.end(), result.begin(), colours.begin(), 8 * sizeof(vertex));

    // std::cout << "offsets ";
    // for(uint i = 1; i < offsets.size(); i++)
    // {
    //     std::cout << offsets[i] << " ";
    // }
    // std::cout << std::endl;

    // parlay::sequence<bool> MIS = get_MIS(G, result, offsets);

    // std::cout << "MIS ";
    // for(uint i = 1; i < MIS.size(); i++)
    // {
    //     std::cout << MIS[i] << " ";
    // }
    // std::cout << std::endl;

    // auto base_clusters = 

