#include <iostream>
#include <parlay/primitives.h>
#include "treeGen.h"
#include <parlay/hash_table.h>
#include <chrono>
#include <cstdlib>
#include "ternarizer.h"
#include "RC.h"
#include "RCdynamic.h"
#include "random_trees.h"
#include "subtree_query.h"

using wedge = std::tuple<long,long,double>;

long graph_size = 100000;
auto distribution = exponential; // either exponential, geometric, constant or linear
double ln = 0.1;
double mean = 20.0;
bool randomized = false;
long num_queries = -1;

static void print_help() {
    std::cout << "Usage: ./program [options]\n"
              << "Options:\n"
              << "  --graph_size <num>    Set graph size (default: 100000)\n"
              << "  --dist <e|g|c|u>      Set distribution: exponential (e), geometric (g), constant (c), uniform (u) (default: e)\n"
              << "  --ln <value>          Set ln parameter (default: 0.1)\n"
              << "  --mean <value>        Set mean value (default: 20.0)\n"
              << "  --randomized          Enable randomized mode\n"
              << "  --num-queries         Number of subtree queries to perform\n"
              << "  --help                Show this help message\n";
}

void test_simple_subtree_queries(const long& num_queries, parlay::sequence<cluster<long,double>>& clusters,  TreeGen<long, double>& TG,  ternarizer<long, double>& TR)
{
    // assert(clusters.size());
    static int seed = 15213;
    parlay::random_generator gen(clusters.size() + seed++);
    std::uniform_int_distribution<long> dis(0, TG.num_vertices - 1);

    auto assocfunc = [] (double A, double B) {return A + B;};
    double identity = 0.0f;

    auto test_pairs = parlay::filter(parlay::delayed_tabulate(num_queries, [&] (long i) {
        auto r = gen[i];
        long v = dis(r);

        return std::pair<long,long>(TG.random_perm_map[v], TG.random_perm_map[TG.parents[v]]);
    }), [] (auto pr){ return pr.first != pr.second;});

    test_pairs = parlay::random_shuffle(test_pairs);

    // std::cout << "Testing " << test_pairs.size() << " values " << std::endl;

    auto manual_values = parlay::delayed_tabulate(test_pairs.size(), [&] (long i){
        auto& pr = test_pairs[i];
        return manual_subtree_sum(pr.first, pr.second, clusters, TR, identity, assocfunc);
    });

    

    auto subtree_sum_values = parlay::delayed_tabulate(test_pairs.size(), [&] (long i){
        auto& pr = test_pairs[i];
        return subtree_query(pr.first, pr.second, clusters, TR, identity, assocfunc);
    });

    parlay::parallel_for(0, test_pairs.size(), [&] (long i) {
        assert(isNearlyEqual(manual_values[i],subtree_sum_values[i]) && "manual subtree values don't match");
    });
    return;
}

void test_batched_subtree_queries(const long& num_queries, parlay::sequence<cluster<long,double>>& clusters,  TreeGen<long, double>& TG,  ternarizer<long, double>& TR)
{
    // assert(clusters.size());
    static int seed = 15213;
    parlay::random_generator gen(clusters.size() + seed++);
    std::uniform_int_distribution<long> dis(0, TG.num_vertices - 1);

    auto assocfunc = [] (double A, double B) {return A + B;};
    double identity = 0.0f;

    auto test_pairs = parlay::filter(parlay::delayed_tabulate(num_queries, [&] (long i) {
        auto r = gen[i];
        long v = dis(r);
        if(i & 1)
            return std::pair<long,long>(TG.random_perm_map[TG.parents[v]], TG.random_perm_map[v]);

        return std::pair<long,long>(TG.random_perm_map[v], TG.random_perm_map[TG.parents[v]]);
    }), [] (auto pr){ return pr.first != pr.second;});

    test_pairs = parlay::random_shuffle(test_pairs);

    // std::cout << "Testing queries " << test_pairs.size() << " in a batch " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    
    double total_time = 0.0;

    for(int a = 0; a < 10; a++) {
        dis(gen); //resettle the generator 

        auto test_pairs = parlay::filter(parlay::delayed_tabulate(num_queries, [&] (long i) {
            auto r = gen[i];
            long v = dis(r);
            if(i & 1)
                return std::pair<long,long>(TG.random_perm_map[TG.parents[v]], TG.random_perm_map[v]);
    
            return std::pair<long,long>(TG.random_perm_map[v], TG.random_perm_map[TG.parents[v]]);
        }), [] (auto pr){ return pr.first != pr.second;});
    
        test_pairs = parlay::random_shuffle(test_pairs);

        //std::cout << "test pairs " << test_pairs[0].first << " " << test_pairs[0].second << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        auto subtree_sum_values = parlay::tabulate(test_pairs.size(), [&] (long i){
            auto& pr = test_pairs[i];
            return subtree_query(pr.first, pr.second, clusters, TR, identity, assocfunc);
        });
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double, std::milli>(end - start).count();   
    }

    total_time /= 10;
    
    

    std::cout << parlay::internal::init_num_workers() << ",";
    std::cout << clusters.size()/3 << ",";
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
    std::cout << num_queries << ",";
    std::cout << total_time
            << ",";

    total_time = 0.0;
    for(int a = 0; a < 10; a++) {
        dis(gen); //resettle the generator
        auto test_pairs = parlay::filter(parlay::delayed_tabulate(num_queries, [&] (long i) {
            auto r = gen[i];
            long v = dis(r);
            if(i & 1)
                return std::pair<long,long>(TG.random_perm_map[TG.parents[v]], TG.random_perm_map[v]);
    
            return std::pair<long,long>(TG.random_perm_map[v], TG.random_perm_map[TG.parents[v]]);
        }), [] (auto pr){ return pr.first != pr.second;});
    
        test_pairs = parlay::random_shuffle(test_pairs);
        
        auto start_ = std::chrono::high_resolution_clock::now();
        auto subtree_batched_values = batched_subtree_query(test_pairs, clusters, TR, identity, assocfunc);
        auto end_ = std::chrono::high_resolution_clock::now();
        total_time += (std::chrono::duration<double, std::milli>(end_ - start_).count());
    }
    std::cout << total_time/10
            << std::endl;

    // assert(subtree_sum_values.size() == subtree_batched_values.size());

    // for(long i = 0; i < test_pairs.size(); i++)
    // {
    //     auto translated_pair = TR.translate_edge(test_pairs[i].first, test_pairs[i].second);
    //     std::cout << translated_pair.first << "--" << translated_pair.second << " -> " << subtree_sum_values[i] << " " << subtree_batched_values[i] << std::endl;
    // }



    // parlay::parallel_for(0, subtree_sum_values.size(), [&] (long i) {
    //     if(!isNearlyEqual(subtree_sum_values[i], subtree_batched_values[i]))
    //     {
    //         std::cout << subtree_sum_values[i] << " != " << subtree_batched_values[i] << std::endl;
    //         assert("Subtree queries don't match" && false);
    //     }
    // });



}

int main(int argc, char* argv[])
{
    


    // std::cout << "Argc " << argc << std::endl;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph_size" && i + 1 < argc) {
            graph_size = std::stol(argv[++i]);
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
        } else if (arg == "--num-queries" && i + 1 < argc) {
            num_queries = std::stod(argv[++i]);
        } else if (arg == "--mean" && i + 1 < argc) {
            mean = std::stod(argv[++i]);
        } else {
            print_help();
            return -1;
        }
    }

        

    const double min_weight = 0.0;
    const double max_weight = 100.0f;

    // std::cout << "II " << II << std::endl;
    TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true, 0);

    TG.generateInitialEdges();

    auto retedges = TG.getAllEdges();
    
    ternarizer<long, double> TR(graph_size, 0);

    auto ret_edge_modified = TR.add_edges(retedges);

    const long max_degree = 3;
    parlay::sequence<cluster<long, double>> clusters; 
    double defretval = 0.0;
    
    create_base_clusters(clusters, ret_edge_modified.second, max_degree, graph_size * extra_tern_node_factor);

    create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false, false);

    // parlay::parallel_for(0, clusters.size(), [&] (long i) {
    //     clusters[i].partial_sum_complete[0] = 0;
    //     clusters[i].partial_sum_complete[1] = 0;
    //     clusters[i].partial_sum_complete[2] = 0;
    // });

    // test_simple_subtree_queries(clusters.size() / 10, clusters, TG, TR);

    if(num_queries >= 0)
        test_batched_subtree_queries(num_queries, clusters, TG, TR);
    else
    {
        num_queries = 1;
        long multiplier = 5;
        while(num_queries <= graph_size)
        {
            test_batched_subtree_queries(num_queries, clusters, TG, TR);
            num_queries *= multiplier;
            if(multiplier == 5)
                multiplier = 2;
            else
                multiplier = 5;
        }
    }

    deleteRCtree(clusters); 



}