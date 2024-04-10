/*
  This code is a heavily modified version of the hotspots assignment
*/


#include <iostream>
#include <string>

#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/include/parlay/internal/get_time.h"

#include "MIS.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"
#include "/scratch/parlaylib/examples/samplesort.h"

/*
    MAXIMUM DEGREE
    I will figure out how to make this dynamic later
*/

const int max_degree = 3;

// **************************************************************
// Driver
// **************************************************************
using vertex = int;
using nested_seq = parlay::sequence<parlay::sequence<vertex>>;
using graph = nested_seq;
using utils = graph_utils<vertex>;



int main(int argc, char* argv[]) {
    auto usage = "Usage: BFS <n>";
   
    if (argc != 2) {std::cout << usage << std::endl;
        return -1;
    }

    int n = 0;
    graph G;
    try { n = std::stoi(argv[1]); }
    catch (...) {}
    if (n == 0) {
        std::cout << "n should be an integer greater than 1" << std::endl;
        return -1;
    }

    long max_degree_count = 0;
    G = utils::rmat_symmetric_graph(n, 10*n);
    return_degree_capped_graph(G, max_degree);
    G = utils::symmetrize(G);
    utils::print_graph_stats(G);    

    auto colours = colour_graph(G, max_degree);
    int original_colour_size = colours.size();
    // extra colours were:
    std::cout << "unique colours are: " << original_colour_size << std::endl;

    sample_sort(colours); // TODO, replace with bucket sort and get counts of each colour

    // for(uint i = 0; i < colours.size(); i++)
    // {
    //     print_string(colours[i].first);
    // }

    auto mis_vertices = generate_MIS(G, colours);

    std::cout << "Number of nodes in the mis are " << mis_vertices.size() << std::endl;

    return 0;
}