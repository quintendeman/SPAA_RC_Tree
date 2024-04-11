/*
  This code is a heavily modified version of the hotspots assignment
*/


#include <iostream>
#include <string>

#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/include/parlay/internal/get_time.h"

#include "MIS.h"
#include "/scratch/parlaylib/examples/samplesort.h"

/*
    MAXIMUM DEGREE
    I will figure out how to make this dynamic later
*/

const int max_degree = 3;

// **************************************************************
// Driver
// **************************************************************
using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;


int main(int argc, char* argv[]) {
    auto usage = "Usage: BFS <n>";
   
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

    graph G;


    auto parents = generate_tree_graph(n, true, false);

    
    
    parlay::sequence<vertex> colours = parlay::tabulate(n, [&] (vertex i) {return i;});

    colours = six_colour_rooted_tree(parents, colours);

    
    G = convert_parents_to_graph(G, parents);

    utils::print_graph_stats(G);

    colours = convert_6_to_3_tree(G, parents, colours);

    sample_sort(colours);

    std::cout << "Unique colours are: " << parlay::unique(colours).size() << std::endl;


    return 0;
}