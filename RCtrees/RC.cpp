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

const vertex max_degree = 2;


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

    graph G = utils::rmat_symmetric_graph(n, 1*n);

    degree_cap_graph(G, max_degree);

    utils::print_graph_stats(G);

    auto colours = colour_chains_to_logn(G, max_degree);

    is_valid_colouring(G, colours); // Takes a while but dw about it

    


    return 0;
}



// auto parents = generate_tree_graph(n, true, false);
    // vertex invalid_vertex = -1;
    

    
    
    // parlay::sequence<vertex> colours = parlay::tabulate(n, [&] (vertex i) {return i;});

    // auto is_valid = verify_colouring(parents, colours, &invalid_vertex);

    // std::cout << "Initial colouring is " << (is_valid ? "\033[32mVALID\033[0m" : "\033[31mINVALID\033[0m") << std::endl;

    // colours = six_colour_rooted_tree(parents, colours);

    // is_valid = verify_colouring(parents, colours, &invalid_vertex);  

    // std::cout << "Six colouring is " << (is_valid ? "\033[32mVALID\033[0m" : "\033[31mINVALID\033[0m") << std::endl;  
    
    // G = convert_parents_to_graph(G, parents);

    // is_valid = verify_colouring(G, parents, colours, &invalid_vertex);  

    // std::cout << "Six colouring is still " << (is_valid ? "\033[32mVALID\033[0m" : "\033[31mINVALID\033[0m") << std::endl;  

    // // sample_sort(colours);

    // // std::cout << "Unique colours are " << std::endl;

    // // auto six_unique_colours = parlay::unique(colours);
    // // for(uint i = 0; i < six_unique_colours.size(); i++)
    // //     std::cout << six_unique_colours[i] << " ";
    // // std::cout << std::endl;

    // // // utils::print_graph_stats(G);

    // colours = convert_6_to_3_tree(G, parents, colours, true);

    // is_valid = verify_colouring(parents, colours, &invalid_vertex);  

    // std::cout << "Three colouring is " << (is_valid ? "\033[32mVALID\033[0m" : "\033[31mINVALID\033[0m") << std::endl;  
    // if(!is_valid)
    //     std::cout << "invalid index is: " << invalid_vertex << std::endl;

    
    // // for(uint i = 0; i < G.size(); i++)
    // // {
    // //     std::cout << colours[i] << "\t";
    // // }
    // // std::cout << std::endl;

    // // sample_sort(colours);

    // // std::cout << "Unique colours are: " << parlay::unique(colours).size() << std::endl;