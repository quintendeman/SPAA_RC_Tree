#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <cmath>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../include/parlay/internal/get_time.h"
#include "RCdynamic.h"
#include "../examples/samplesort.h"
#include "incMST.h"


// **************************************************************
// Driver
// **************************************************************
using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

using datatype = double;

const vertex max_degree = 3;

const double espilon = 0.001;

static bool isNearlyEqual(double a, double b, double epsilon = espilon) {
    return std::abs(a - b) < epsilon;
}


// This wipes out all existing weights, replacing them with random ones
// The lambda function is MAX
void insert_random_weights(const parlay::sequence<vertex>& parents, parlay::sequence<cluster<vertex, datatype>>& clusters)
{
    std::random_device rd;                         // Seed generator
    std::mt19937 gen(rd());                        // Standard Mersenne Twister engine seeded with rd()
    std::uniform_int_distribution<vertex> dist(0, clusters.size() * 100);   

    parlay::sequence<std::tuple<vertex, vertex, datatype>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, datatype>(i, parents[i], (datatype) ( dist(gen)));
    }); 
    
    weighted_edges = parlay::filter(weighted_edges, [&] (auto edge) {
        return std::get<0>(edge) != std::get<1>(edge);
    });

    batchModifyEdgeWeights(weighted_edges, [] (datatype a, datatype b) {
        return a > b ? a : b; // MAX (not min)
    }, clusters);

    return;
}

// Warning, tree properties not preserved but degree preserved (hopefully)
parlay::sequence<std::tuple<vertex,vertex,double>> generate_random_edges(const parlay::sequence<vertex>& parents, bool do_just_one = false)
{
    parlay::sequence<std::tuple<vertex,vertex,double>> new_edges;
    std::random_device rd;                         // Seed generator
    std::mt19937 gen(rd());                        // Standard Mersenne Twister engine seeded with rd()
    std::uniform_int_distribution<vertex> dist(0, parents.size() * 100);   

    parlay::sequence<std::atomic<unsigned char>> counts = parlay::sequence<std::atomic<unsigned char>>(parents.size());
    
    parlay::parallel_for(0, counts.size(), [&] (vertex i) {
        counts[i] = 0;
    });

    parlay::parallel_for(0, counts.size(), [&] (vertex i) {
        if(parents[i] != i)
        {
            counts[parents[i]].fetch_add(1);
            counts[i].fetch_add(1);
        }
    });


    int num_samples = 3 + (rand() % 100);
    new_edges = parlay::tabulate(parents.size() / num_samples, [&] (vertex i) {
        auto my_index = i * num_samples;
        auto parent_index = my_index - num_samples/2;
        if(parent_index < 0 || my_index >= parents.size() || parent_index == my_index || counts[my_index] == max_neighbours || counts[parent_index] == max_neighbours)
            return std::tuple<vertex,vertex,double>(0, 0, 0);
        return std::tuple<vertex, vertex, double> (my_index, parent_index, (double)dist(gen));
    });
    new_edges = parlay::filter(new_edges, [&] (auto& edge) {
        return std::get<0>(edge) != std::get<1>(edge);
    });
    parlay::sequence<std::tuple<vertex,vertex,double>> ret_edges;
    if(new_edges.size() && do_just_one)
    {
        ret_edges.push_back(new_edges[0]);
        ret_edges.push_back(new_edges[new_edges.size()-1]);
        return ret_edges;
    }
    return new_edges;
}

int main(int argc, char* argv[]) {
    auto usage = "Usage: RC [--graph-size <graph-size>] [-n <graph-size>] [--num-queries <num-queries>] [--print-creation] [--do-height <true|false>] [--randomized <true|false>] [--do-just-one]";

    srand(time(NULL));

    const vertex max_rand_size = 100000000l;
    // vertex graph_size = rand() % max_rand_size; 
    vertex graph_size = 1; 
    for(unsigned short i = 0; i < 3; i++) // random but leaning more towards higher values
    {
        vertex tg = rand() % max_rand_size;
        if (tg > graph_size)
            graph_size = tg;
    }
    
    if((graph_size % 2) == 1)
        graph_size++;
    vertex num_queries = 100; // Default value
    bool print_creation = false;
    bool randomized = false; // Default value
    bool do_just_one = false;

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
        } else if (arg == "--do-just-one") {
            do_just_one = true;
        } 
        else if (arg == "--randomized" && i + 1 < argc) {
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
    // create_RC_tree(clusters, graph_size, randomized);

    double defretval = 0.0;
    create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;}, randomized);

    auto end_creation = std::chrono::high_resolution_clock::now();

    insert_random_weights(parents, clusters);

    checkEdgePtrValid(clusters);

    auto new_edges_as_truples = generate_random_edges(parents, do_just_one);

    if(graph_size <= 100)
    {
        std::cout << green << "Old edges: " << std::endl;
        for(auto& edge : convertClustersToTruples(clusters))
            std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
        std::cout << reset << std::endl;
    }

    if(graph_size <= 100 || do_just_one)
    {
        std::cout << blue << "New edges: " << std::endl;
        for(auto& edge: new_edges_as_truples)
            std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
        std::cout << reset << std::endl;
    }

    // parlay::sequence<vertex> edge_counts = parlay::sequence<vertex>(clusters.size(), 0); // scratch space to be allocated once
    // parlay::sequence<vertex> vertex_counts = parlay::sequence<vertex>(clusters.size(), 0);

    std::cout << "There are " << new_edges_as_truples.size() << " new edges " << std::endl;

    incrementMST(clusters, new_edges_as_truples);

    std::chrono::duration<double> creation_time = end_creation - start_creation;
    
    
    if(graph_size <= 100)
        printTree(clusters);

    // if(graph_size <= 100)
    //     printTree(clusters);

    deleteRCtree(clusters);

    if (print_creation) {
        std::cout << graph_size << "," << std::setprecision(6) << creation_time.count() << std::endl;
    } else {
        std::cout << "Creation time: " << std::setprecision(6) << creation_time.count() << " seconds" << std::endl;
    }

    return 0;
}
