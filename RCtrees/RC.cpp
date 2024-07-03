#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../include/parlay/internal/get_time.h"
#include "RCdynamic.h"
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

const double espilon = 0.001;

static bool isNearlyEqual(double a, double b, double epsilon = espilon) {
    return std::abs(a - b) < epsilon;
}

void test_rc_valid(parlay::sequence<vertex>& parents, parlay::sequence<cluster<vertex, datatype>>& clusters)
{
    parlay::sequence<std::tuple<vertex, vertex, datatype>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, datatype>(i, parents[i], (datatype) ( clusters.size()));
    });

    // std::cout << "Working with " << blue << weighted_edges.size() << reset << " edges" << std::endl;

    // pick a random pair to have a path between
    vertex starting_index, ending_index;
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> disint(0, clusters.size()-1);
    starting_index = disint(gen);

    // chase edge until it hits a dead end then stop
    vertex current_index = starting_index;
    vertex parent_index; // anything other than current index will work

    auto count = 0;

    std::cout << "Test Path: " << bold << bright_yellow << current_index << reset;



    while(true) // we hit a root
    {
        auto edge_tuple = weighted_edges[current_index];
        parent_index = std::get<1>(edge_tuple);
        if(parent_index == current_index)
            break;
        weighted_edges[current_index] = std::tuple<vertex, vertex, datatype>(current_index, parent_index, (datatype) 1.0f);
        std::cout << bold << bright_yellow << " -> " << parent_index << reset;
        current_index = parent_index;
        if(count > (10 * clusters.size()))
            {std::cout << red << "Cycle?" << std::endl; exit(1);}
        count++;
    }


    std::cout << std::endl;



    ending_index = parent_index;

    // std::cout << blue << "Checking from " << ending_index << " to " << starting_index << reset << std::endl;


    weighted_edges = parlay::filter(weighted_edges, [&] (auto edge) {
        return std::get<0>(edge) != std::get<1>(edge);
    });

    batchModifyEdgeWeights(weighted_edges, [] (datatype a, datatype b) {
        return a + b;
    }, clusters);

    auto retval = PathQuery(&clusters[starting_index], &clusters[ending_index], (datatype) -1, [] (datatype a, datatype b) {
        return a + b;
    });

    if(isNearlyEqual(retval, count) || (starting_index == ending_index && retval == -1))
        std::cout << green << "The count should have been " << (starting_index == ending_index ? -1 : count) << " and it was " << retval << reset << std::endl;
    else
        std::cout << red << "The count should have been " << count << " but it was " << retval << reset << std::endl;

    return;
}



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


    parlay::sequence<cluster<vertex, datatype>> clusters;

    // Measure creation time
    auto start_creation = std::chrono::high_resolution_clock::now();
    create_base_clusters(G, clusters, max_degree);
    create_RC_tree(clusters, graph_size, randomized);
    auto end_creation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> creation_time = end_creation - start_creation;


    test_rc_valid(parents, clusters);
    


    if(graph_size <= 100)
        printTree(clusters);

    auto ret_adj = getAdjacencyAtLevel(&clusters[(vertex) graph_size/2], 0, clusters);

    for(auto& i : ret_adj)
        std::cout << i << " ";
    std::cout << std::endl;

    ret_adj = getAdjacencyAtLevel(&clusters[(vertex) graph_size/2], 1, ret_adj, 0, clusters);

    for(auto& i : ret_adj)
        std::cout << i << " ";
    std::cout << std::endl;


    deleteRCtree(clusters);

    if (print_creation) {
        std::cout << graph_size << "," << std::setprecision(6) << creation_time.count() << std::endl;
    } else {
        std::cout << "Creation time: " << std::setprecision(6) << creation_time.count() << " seconds" << std::endl;
    }

    return 0;
}


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