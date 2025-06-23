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
#include "path_query.h"

using wedge = std::tuple<long,long,double>;

static void print_help() {
    std::cout << "Usage: ./program [options]\n"
              << "Options:\n"
              << "  --graph_size <num>    Set graph size (default: 100000)\n"
              << "  --dist <e|g|c|u>      Set distribution: exponential (e), geometric (g), constant (c), uniform (u) (default: e)\n"
              << "  --ln <value>          Set ln parameter (default: 0.1)\n"
              << "  --mean <value>        Set mean value (default: 20.0)\n"
              << "  --randomized          Enable randomized mode\n"
              << "  --num-insertions      Number of (add+delete) edges i.e. batch insertion size -1 for sweep\n"
              << "  --help                Show this help message\n";
}

int main(int argc, char* argv[])
{
    
    long graph_size = 100000;
    auto distribution = constant; // either exponential, geometric, constant or uniform
    double ln = 0.1;
    double mean = 1.01f;
    bool randomized = false;
    long batch_insertion = -1;

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
        } else if (arg == "--num-insertions") {
          batch_insertion = std::stol(argv[++i]);;
        } else if (arg == "--ln" && i + 1 < argc) {
            ln = std::stod(argv[++i]);
        } else if (arg == "--mean" && i + 1 < argc) {
            mean = std::stod(argv[++i]);
        }
    }
        
    const double min_weight = 0.0;
    const double max_weight = 100.0f;

    // std::cout << "Graph size " << graph_size << std::endl;

    double dynamic_gen_time = 0.0f;
    double dynamic_tern_time = 0.0f;
    int dynamic_counter = 0;

    const int maxII = 4;

    for(int II = 0; II < maxII; II++) // settle the fragmentation
    {

        TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true, II);

        TG.generateInitialEdges();

        auto retedges = TG.getAllEdges();

        ternarizer<long, double> TR(graph_size, 0);

        auto ret_edge_modified = TR.add_edges(retedges);

        const long max_degree = 3;
        parlay::sequence<cluster<long, double>> clusters; 
        double defretval = 0.0;

        parlay::sequence<wedge> empty_edges_sequence;

        create_base_clusters(clusters, ret_edge_modified.second, max_degree, graph_size * extra_tern_node_factor);

        create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);

        deleteRCtree(clusters); 
    }

    // now for the dynamic stuff

    TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true, 12345);

    TG.generateInitialEdges();

    auto retedges = TG.getAllEdges();

    ternarizer<long, double> TR(graph_size, 0);

    auto ret_edge_modified = TR.add_edges(retedges);

    const long max_degree = 3;
    parlay::sequence<cluster<long, double>> clusters; 
    double defretval = 0.0;

    parlay::sequence<wedge> empty_edges_sequence;

    create_base_clusters(clusters, ret_edge_modified.second, max_degree, graph_size * extra_tern_node_factor);

    create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);

    // std::cout << "Num dynamic edges " << TG.get_num_dynamic_edges() << std::endl;
    if(batch_insertion != -1)
    {
      double desired_prob = (static_cast<double>(batch_insertion)/TG.get_num_dynamic_edges());
      if(desired_prob > 1.0f)
        desired_prob = 1.0f;

      parlay::sequence<std::pair<long, long>> empty_pairs;
      long batch_insert_size = 0;
      auto del_pairs = TG.generateDeleteEdges(desired_prob);
      auto add_truples = TG.generateAddEdges(desired_prob);
      batch_insert_size = del_pairs.size() + add_truples.size();
      auto delete_start = std::chrono::high_resolution_clock::now();
      auto del_pairs__ = TR.delete_edges(del_pairs);
      auto add_pairs__ = TR.add_edges(add_truples);
      // add_truples.append(ret_seqs.first); 
      auto delete_end = std::chrono::high_resolution_clock::now();
      batchInsertEdge(del_pairs__.first, del_pairs__.second, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
      batchInsertEdge(add_pairs__.first, add_pairs__.second, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
      auto final_insert_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = final_insert_end - delete_start;
      
      std::cout << parlay::internal::init_num_workers() << ",";
      std::cout << graph_size << "," << ln << "," << mean << ",";
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
      std::cout << batch_insert_size << "," << diff.count() << ",";

      std::cout << std::endl;
    }
    else
    {
      batch_insertion = 1;
      long multiplier = 5;
      while(batch_insertion <= TG.get_num_dynamic_edges() )
      {
        double desired_prob = (static_cast<double>(batch_insertion)/TG.get_num_dynamic_edges());
        if(desired_prob > 1.0f)
          desired_prob = 1.0f;

        parlay::sequence<std::pair<long, long>> empty_pairs;
        long batch_insert_size = 0;
        auto del_pairs = TG.generateDeleteEdges(desired_prob);
        auto add_truples = TG.generateAddEdges(desired_prob);
        batch_insert_size = del_pairs.size() + add_truples.size();
        auto delete_start = std::chrono::high_resolution_clock::now();
        auto del_pairs__ = TR.delete_edges(del_pairs);
        auto add_pairs__ = TR.add_edges(add_truples);
        // add_truples.append(ret_seqs.first); 
        auto delete_end = std::chrono::high_resolution_clock::now();
        batchInsertEdge(del_pairs__.first, del_pairs__.second, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
        batchInsertEdge(add_pairs__.first, add_pairs__.second, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
        auto final_insert_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = final_insert_end - delete_start;
        
        std::cout << parlay::internal::init_num_workers() << ",";
        std::cout << graph_size << "," << ln << "," << mean << ",";
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
        std::cout << batch_insert_size << "," << diff.count() << ",";

        std::cout << std::endl;
        batch_insertion *= multiplier;
        if (multiplier == 2)
          multiplier = 5;
        else
          multiplier = 2;
      }
        double desired_prob = (static_cast<double>(batch_insertion)/TG.get_num_dynamic_edges());
        if(desired_prob > 1.0f)
          desired_prob = 1.0f;

        parlay::sequence<std::pair<long, long>> empty_pairs;
        long batch_insert_size = 0;
        auto del_pairs = TG.generateDeleteEdges(desired_prob);
        auto add_truples = TG.generateAddEdges(desired_prob);
        batch_insert_size = del_pairs.size() + add_truples.size();
        auto delete_start = std::chrono::high_resolution_clock::now();
        auto del_pairs__ = TR.delete_edges(del_pairs);
        auto add_pairs__ = TR.add_edges(add_truples);
        // add_truples.append(ret_seqs.first); 
        auto delete_end = std::chrono::high_resolution_clock::now();
        batchInsertEdge(del_pairs__.first, del_pairs__.second, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
        batchInsertEdge(add_pairs__.first, add_pairs__.second, clusters, (double) 0.0f, [] (double A, double B) {return A+B;}, randomized);
        auto final_insert_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = final_insert_end - delete_start;
        
        std::cout << parlay::internal::init_num_workers() << ",";
        std::cout << graph_size << "," << ln << "," << mean << ",";
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
        std::cout << batch_insert_size << "," << diff.count() << ",";

        std::cout << std::endl;
    }

    deleteRCtree(clusters); 

}