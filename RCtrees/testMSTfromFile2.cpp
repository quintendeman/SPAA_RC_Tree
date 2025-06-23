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

// Simple text progress bar
void print_progress(double fraction, double mbps) {
    const int width = 50;
    int pos = static_cast<int>(width * fraction);
    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(fraction * 100.0) << "% ("
              << mbps << " MB/s)\r";
    std::cout.flush();
}

void remove_duplicates(parlay::sequence<std::tuple<long,long,double>>& input) {
    parlay::sequence<std::atomic<char>> flags = parlay::tabulate<std::atomic<char>>(roundUpPow2(input.size() * 1.5f), [] (long i) {return 0;});
    auto vals = parlay::sequence<std::tuple<long, long, double>>(flags.size());

    parlay::parallel_for(0, input.size(), [&] (long i) {
        long& v = std::get<0>(input[i]);
        long& w = std::get<1>(input[i]);
        double& weight = std::get<2>(input[i]);

        if(w < v)
            std::swap(v, w);
        long expected_location = simple_hash(v, w) & (flags.size() - 1);
        while(true) {
            if(flags[expected_location] == 0) {
                long& v_ = std::get<0>(vals[expected_location]);
                long& w_ = std::get<1>(vals[expected_location]);
                double& weight_ = std::get<2>(vals[expected_location]);
                if(v_ == v && w_ == w)
                    return;
                char expected = 0;
                char desired = 1;
                bool success = flags[expected_location].compare_exchange_strong(expected, desired);
                if(success) {
                    v_ = w;
                    w_ = w;
                    weight_ = weight;
                    return;
                }
            }
            expected_location = (expected_location + 1) & (flags.size() - 1);
        }
    });

    auto indices = parlay::delayed_tabulate(flags.size(), [&] (long i){
        return i;
    });
    auto good_indices = parlay::filter(indices, [&] (long i) {
        return flags[indices[i]] == 1;
    });
    input = parlay::tabulate(good_indices.size(), [&] (long i) {
        return vals[good_indices[i]];
    });


}


std::pair<parlay::sequence<std::tuple<long, long, double>>, vertex>
read_wedge_from_file(const std::string &file_path) {
    std::ifstream file(file_path, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + file_path);
    }

    // Get total file size in bytes
    std::streamsize total_size = file.tellg();
    file.seekg(0, std::ios::beg);

    parlay::sequence<std::tuple<long, long, double>> edge_list;
    vertex max_vertex = 0;

    std::string line;
    std::streamsize processed = 0;
    auto last_time = std::chrono::steady_clock::now();
    std::streamsize last_processed = 0;

    long time = -1;

    while (std::getline(file, line)) {
        processed = file.tellg();
        if (processed < 0) processed = total_size; // on final line

        // Compute progress every 100k lines or so
        static size_t counter = 0;
        if (++counter % 100000 == 0) {
            double fraction = double(processed) / double(total_size);
            auto now = std::chrono::steady_clock::now();
            double dt = std::chrono::duration<double>(now - last_time).count();
            double dbytes = double(processed - last_processed);
            double mbps = (dbytes / (1024.0 * 1024.0)) / (dt > 0 ? dt : 1);
            // print_progress(fraction, mbps);
            last_time = now;
            last_processed = processed;
        }

        std::istringstream iss(line);
        std::string tag;
        long u, v;
        double w;
        if (!(iss >> u >> v )) continue;
        // std::cout << u << " " << v << " " << w << std::endl;
        // if (tag != "a") continue;
        // if (u >= v) continue;

        edge_list.emplace_back(u, v, (double) (u * v + u + v));
        max_vertex = std::max({max_vertex, u, v});
        assert(v >= 0 && w >= 0 && "indices should be positive");
    }

    // Final update to 100%
    // print_progress(1.0, 0.0);
    // std::cout << std::endl;

    // std::cout << "Duplicate removal reduced size from " << edge_list.size() << " to ";
    remove_duplicates(edge_list);
    // std::cout << edge_list.size() << std::endl;

    edge_list = parlay::random_shuffle(edge_list);

    parlay::parallel_for(0, edge_list.size(), [&] (long i) { //
        std::get<2>(edge_list[i]) = edge_list.size() * 2 + -1 * i;
    });

    return {edge_list, max_vertex};
}





int main(int argc, char* argv[]) {
    auto usage = "Usage: ./testMST.out [path_to_file]\nThe file must have data in the form of \"i j w\" per line where i and j are indices and w is an integer weight";

    srand(time(NULL));


    // vertex graph_size = rand() % max_rand_size; 
    

   
    assert(argc == 3 && "This expect only two arguments, the first is a path to a file and the second is the batch size");

    std::string path_to_file = std::string(argv[1]);

    long batch_size = std::stol(argv[2]);

    auto wedge_pair = read_wedge_from_file(path_to_file);
    // std::cout << "max node: " << wedge_pair.second << std::endl;
    // std::cout << "num edges: " << wedge_pair.first.size() << std::endl;

    vertex graph_size = wedge_pair.second + 1;  

    
    ternarizer<long, double> TR(graph_size, 0.0f);
    parlay::sequence<cluster<vertex, double>> clusters;
    
    parlay::sequence<std::tuple<long,long,double>> empty_edges;
    
    create_base_clusters(clusters, empty_edges, static_cast<long>(3), graph_size*extra_tern_node_factor);
    create_RC_tree(clusters, graph_size*extra_tern_node_factor, (double) 0.0f, [] (double A, double B) {return A+B;}, false);
    // incrementMST(clusters, wedge_pair.first, TR);

    // std::cout << "Created empty clusters" << std::endl;
    

    const long max_edges = wedge_pair.first.size();
    long current_edges = 0;
    parlay::internal::timer t;

    // assert(max_edges < graph_size);

    
    const long starting_insert_edges = 1000;


    // for(auto& wd : wedge_pair.first)
    //     std::cout << std::get<0>(wd) << " " << std::get<1>(wd) << " " << std::get<2>(wd) << std::endl;



    incrementMST(clusters, wedge_pair.first.subseq(0, starting_insert_edges), TR);
    
    parlay::parallel_for(0, wedge_pair.first.size(), [&] (long i) {
        std::get<2>(wedge_pair.first[i]) = -1 * i;
    });


    long times_done = 0;
    long max_times_done = 5;
    long warmup_trials = 0;
    // if(batch_size >= 1000) {
    //     max_times_done = 100;
    //     warmup_trials = 50;
    // }
    // if(batch_size >= 100000l) {
    //     max_times_done = 5;
    //     warmup_trials = 0;
    // }  if (batch_size >= 1000000l) {
    //     max_times_done = 1;
    //     warmup_trials = 0;
    // }

    // std::cout << "max_times_done " << max_times_done << std::endl;
    // std::cout << "warmup_trials " << warmup_trials << std::endl;

    // time measurement start
    
    double total_time = 0.0f;

    current_edges = starting_insert_edges;

    while(current_edges < max_edges) {
        // init start

        // std::cout << "CE " << current_edges << std::endl;

        long start_index = current_edges;
        long end_index = (current_edges + batch_size) >= max_edges ? max_edges : current_edges + batch_size;
        auto sub_seq = wedge_pair.first.subseq(start_index, end_index);
        // init end

        // comp start
        auto start = std::chrono::high_resolution_clock::now();

        incrementMST(clusters, sub_seq, TR);


        // comp end
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        if(times_done >= warmup_trials)
            total_time += elapsed.count(); 

        current_edges += batch_size;
        times_done++;
        if(times_done == max_times_done)
            break;
    }


   
    double time_per_edge = total_time / static_cast<double>((times_done - warmup_trials) * batch_size);
    // std::cout << "Time per edge: " << time_per_edge << " seconds" << std::endl;
    std::cout << parlay::internal::init_num_workers() << batch_size << "," << time_per_edge << std::endl;

    deleteRCtree(clusters);

    // std::cout << "Working with a graph of size " << graph_size << std::endl; 

    // TreeGen<long, double> TG(graph_size, 0.0f, 1.0f, 0.1, 1.1, exponential, true, 0);

    // TG.generateInitialEdges();

    // auto retedges = TG.getAllEdges();

    // ternarizer<long, double> TR(graph_size, 0.0f);

    // auto ternerized_initial_edges = TR.add_edges(retedges);
    
    // parlay::sequence<cluster<vertex, datatype>> clusters;

    // // // Measure creation time
    // // auto start_creation = std::chrono::high_resolution_clock::now();
    // create_base_clusters(clusters, ternerized_initial_edges, static_cast<long>(3), graph_size*extra_tern_node_factor);
    // create_RC_tree(clusters, graph_size*extra_tern_node_factor, static_cast<double>(0.0f), [] (double A, double B) {return A+B;}, false);

    // if(num_additions >= 0)
    // {

    //     std::cout << parlay::internal::init_num_workers() << ",";
    //     std::cout << graph_size << ",";
    //     std::cout << ln << "," << mean << ",";
    //         switch(distribution){
    //         case exponential:
    //             std::cout << "e";
    //             break;
    //         case geometric:
    //             std::cout << "g";
    //             break;
    //         case constant:
    //             std::cout << "c";
    //             break;
    //         case uniform:
    //             std::cout << "u";
    //             break;
    //         }
    //     std::cout << ",";
    //     std::cout << num_additions << ",";

    //     auto random_add_edges = generate_random_edges(TG.parents, TG.random_perm_map, num_additions);

    //     incrementMST(clusters, random_add_edges, TR);
    // }
    // else
    // {
    //     // sweep
    //     long multiplier = 5;
    //     num_additions = 1;

    //     while(num_additions <= graph_size)
    //     {
    //         std::cout << parlay::internal::init_num_workers() << ",";
    //         std::cout << graph_size << ",";
    //         std::cout << ln << "," << mean << ",";
    //             switch(distribution){
    //             case exponential:
    //                 std::cout << "e";
    //                 break;
    //             case geometric:
    //                 std::cout << "g";
    //                 break;
    //             case constant:
    //                 std::cout << "c";
    //                 break;
    //             case uniform:
    //                 std::cout << "u";
    //                 break;
    //             }
    //         std::cout << ",";
    //         std::cout << num_additions << ",";
    //         auto random_add_edges = generate_random_edges(TG.parents, TG.random_perm_map, num_additions);
    //         incrementMST(clusters, random_add_edges, TR);

    //         num_additions *= multiplier;
    //         if(multiplier == 5)
    //             multiplier = 2;
    //         else
    //             multiplier = 5;
    //     }

    // }
        



    // // if(graph_size <= 100)
    // //     printTree(clusters);

    // // if(graph_size <= 100)
    // //     printTree(clusters);

    // deleteRCtree(clusters);


    return 0;
}
