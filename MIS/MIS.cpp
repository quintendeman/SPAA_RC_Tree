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
// using utils = graph_utils<vertex>;



int main(int argc, char* argv[]) {
    auto usage = "Usage: BFS <n>";
   
    if (argc != 2) {std::cout << usage << std::endl;
        return -1;
    }

    long n = 0;
    try { n = std::stoi(argv[1]); }
    catch (...) {}
    if (n == 0) {
        std::cout << "n should be an integer greater than 1" << std::endl;
        return -1;
    }


    auto parents = generate_tree_graph(n, true, false);

    
    
    parlay::sequence<long> colours = parlay::tabulate(n, [&] (long i) {return i;});

    colours = six_colour_rooted_tree(parents, colours);

    sample_sort(colours);

    std::cout << "Unique colours are: " << parlay::unique(colours).size() << std::endl;
    

    return 0;
}