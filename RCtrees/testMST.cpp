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
#include "utils.h"
#include "random_trees.h"
#include "treeGen.h"
#include "ternarizer.h"

// **************************************************************
// Driver
// **************************************************************
using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

using datatype = double;

const vertex max_degree = 3;

parlay::sequence<std::tuple<vertex,vertex,double>> generate_random_edges(const parlay::sequence<vertex>& parents, const parlay::sequence<vertex>& map, long num_gens = 1)
{
    static int seed = 0;
    parlay::random_generator gen(seed++);
    std::uniform_int_distribution<long> dis(0, parents.size()-1);
    std::uniform_real_distribution<double> dis_ur(0.0f, 1.0f);

    static auto random_map = parlay::random_permutation(parents.size());

    return std::move(parlay::filter(parlay::tabulate(num_gens, [&] (long i) {
        auto r = gen[i];
        long child = i;
        long parent = random_map[i];
        double random_weight = dis_ur(r);
        return std::tuple<long, long, double>(child, parent, random_weight);
    }), [] (auto tpl){return std::get<0>(tpl) != std::get<1>(tpl);}));

}

int main(int argc, char* argv[]) {
    auto usage = "Usage: testMST.out [--graph_size <graph_size>] [-n <graph_size>] [--num-additions <extra edges added>] [--ln <0.0-1.0>] [--mean <1.0f+>] [--dist <u,e,g,c>] [--randomized]\nPrints compressed tree creation time, MST time, insertion time";

    srand(time(NULL));


    // vertex graph_size = rand() % max_rand_size; 
    vertex graph_size = 10000; 
    vertex num_additions = -1;
    auto distribution = exponential; // either exponential, geometric, constant or linear
    double ln = 0.1;
    double mean = 20.0;
    bool randomized = false;
   

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph_size" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
        } else if (arg == "--num-additions" && i + 1 < argc) {
            num_additions = std::stol(argv[++i]);
        } else if (i == 1 && arg.find_first_not_of("0123456789") == std::string::npos) {
            // Handle the case where a single number is provided as graph_size
            graph_size = std::stol(arg);
        } else if (arg == "--dist" && i + 1 < argc) {
            std::string dist_arg = argv[++i];
            if (dist_arg == "e") {
                distribution = exponential;
            } else if (dist_arg == "g") {
                distribution = geometric;
            } else if (dist_arg == "c") {
                distribution = constant;
            } else if (dist_arg == "u") {
                distribution = uniform;
            } else {
                std::cerr << "Invalid distribution type\n";
                return 1;
            }
        } else if (arg == "--randomized") {
          randomized = true;
        } else if (arg == "--ln" && i + 1 < argc) {
            ln = std::stod(argv[++i]);
        } else if (arg == "--mean" && i + 1 < argc) {
            mean = std::stod(argv[++i]);
        } else {
            std::cout << "Invalid argument: " << arg << std::endl;
            std::cout << usage << std::endl;
            return -1;
        }
    }

    if (graph_size < 0) {
        std::cout << "graph_size should be an integer greater than 1" << std::endl;
        return -1;
    }


    // std::cout << "Working with a graph of size " << graph_size << std::endl; 

    TreeGen<long, double> TG(graph_size, 0.0f, 1.0f, 0.1, 1.1, exponential, true, 0);

    TG.generateInitialEdges();

    auto retedges = TG.getAllEdges();

    ternarizer<long, double> TR(graph_size, 0.0f);

    auto ternerized_initial_edges = TR.add_edges(retedges);
    
    parlay::sequence<cluster<vertex, datatype>> clusters;

    // // Measure creation time
    // auto start_creation = std::chrono::high_resolution_clock::now();
    create_base_clusters(clusters, ternerized_initial_edges.second, static_cast<long>(3), graph_size*extra_tern_node_factor);
    create_RC_tree(clusters, graph_size*extra_tern_node_factor, static_cast<double>(0.0f), [] (double A, double B) {return A+B;}, false);

    if(num_additions >= 0)
    {

        std::cout << parlay::internal::init_num_workers() << ",";
        std::cout << graph_size << ",";
        std::cout << ln << "," << mean << ",";
            switch(distribution){
            case exponential:
                std::cout << "e";
                break;
            case geometric:
                std::cout << "g";
                break;
            case constant:
                std::cout << "c";
                break;
            case uniform:
                std::cout << "u";
                break;
            }
        std::cout << ",";
        std::cout << num_additions << ",";

        double cst_time = 0.0;
        double mst_time = 0.0;
        double insert_time = 0.0;
        for(short a = 0; a < 10; a++) {
            auto random_add_edges = generate_random_edges(TG.parents, TG.random_perm_map, num_additions);
            auto rets = incrementMST(clusters, random_add_edges, TR);
            double c,m,i;
            std::tie(c,m,i) = rets;
            cst_time+=c;
            mst_time+=m;
            insert_time+=i;
        }
        std::cout << cst_time/10 << "," << mst_time/10 << "," << insert_time;

    }
    else
    {
        // sweep
        long multiplier = 5;
        num_additions = 1;

        while(num_additions <= graph_size)
        {
            std::cout << parlay::internal::init_num_workers() << ",";
            std::cout << graph_size << ",";
            std::cout << ln << "," << mean << ",";
                switch(distribution){
                case exponential:
                    std::cout << "e";
                    break;
                case geometric:
                    std::cout << "g";
                    break;
                case constant:
                    std::cout << "c";
                    break;
                case uniform:
                    std::cout << "u";
                    break;
                }
            std::cout << ",";
            std::cout << num_additions << ",";
            double cst_time = 0.0;
            double mst_time = 0.0;
            double insert_time = 0.0;
            for(short a = 0; a < 10; a++) {
                auto random_add_edges = generate_random_edges(TG.parents, TG.random_perm_map, num_additions);
                auto rets = incrementMST(clusters, random_add_edges, TR);
                double c,m,i;
                std::tie(c,m,i) = rets;
                cst_time+=c;
                mst_time+=m;
                insert_time+=i;
            }
            std::cout << cst_time/10 << "," << mst_time/10 << "," << insert_time;
            std::cout << std::endl;

            num_additions *= multiplier;
            if(multiplier == 5)
                multiplier = 2;
            else
                multiplier = 5;
        }

    }
        



    // if(graph_size <= 100)
    //     printTree(clusters);

    // if(graph_size <= 100)
    //     printTree(clusters);

    deleteRCtree(clusters);


    return 0;
}
