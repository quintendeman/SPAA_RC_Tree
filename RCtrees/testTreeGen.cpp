#include <iostream>
#include "treeGen.h"

#include <chrono>

#include "RC.h"
#include "RCdynamic.h"
#include "random_trees.h"
#include <ctime>    // For time()


using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

int main(int argc, char* argv[])
{
    std::srand(std::time(0));
    vertex graph_size;

    // std::cout << "Argc " << argc << std::endl;

    if(argc >= 2)
    {
        std::string arg = argv[1];
        if (arg.find_first_not_of("0123456789") == std::string::npos)
        {
            graph_size = std::stol(argv[1]);
        }
        else
        {
        graph_size = 100000000 + (rand() % 400000000l);    
        }
    }
    else
        graph_size = 100000000 + (rand() % 400000000l);    
        
    


    double min_weight = 0.0;
    double max_weight = 100.0f;
    // // graph_size+=1000;
    // graph_size = 100 + (rand() % 100);
    // graph_size = 30;

    std::cout << "Working with graph_size " << graph_size << std::endl;

    for(unsigned I = 0; I < 10; I ++)
    {
        std::cout << "Testing graph " << I+1 << std::endl;

        TreeGen<vertex, double> TG(graph_size, min_weight, max_weight, 0.9, 5, uniform);

        TG.generateInitialEdges();

        auto retedges = TG.getAllEdges();


        // std::cout << "OG edges " << std::endl;
        // for(auto& edg : retedges)
        //     std::cout << std::get<0>(edg) << " -- " << std::get<1>(edg) << ": " << std::get<2>(edg) << std::endl;
        // std::cout << std::endl;

        // std::cout << "original counts ";
        // for(unsigned int i = 0; i < TG.counts.size(); i++)
        // {
        //     std::cout << i << "/" << TG.counts[i] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << "original parents ";
        // for(unsigned int i = 0; i < TG.parents.size(); i++)
        // {
        //     std::cout << i << "->" << TG.parents[i] << " ";
        // }
        // std::cout << std::endl;
    

        TG.verifyParents();

        // return 0;
        

        // graph G;
        // G = convert_parents_to_graph(G, TG.parents);
        const vertex max_degree = 3;
        
        parlay::sequence<cluster<vertex, double>> clusters; 


        auto static_creation_start = std::chrono::high_resolution_clock::now(); //

        create_base_clusters(clusters, retedges, max_degree, graph_size);
        
        // printTree(clusters);
        double defretval = 0.0;

        create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);

        auto static_creation_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double>  dur = static_creation_end - static_creation_start;
        std::cout << "static elapsed time: " << dur.count() << " seconds\n";


        // return 0;

        // create_base_clusters(G, clusters, max_degree);

        // TG.clearWeights(0.0);



        // return 0;

        // std::cout << "deleted " << std::endl;
        // for(auto& pr : deledges)
        //     std::cout << pr.first << " " << pr.second << std::endl;
        // std::cout << std::endl;

        // std::cout << "after delete counts ";
        // for(auto& ct : TG.counts)
        //     std::cout << ct << " ";
        // std::cout << std::endl;

        std::cout << TG.interconnects.size() << " number of dynamic edges " << std::endl;

        for(unsigned int i = 0; i < 20; i++)
        {
            std::cout << "testing batchinsert " << i+1 << std::endl;
            double del_prob = (rand() % 100);
            del_prob/=100;
            double add_prob = (rand() % 100);
            add_prob/=100;
            
            auto deledges = TG.generateDeleteEdges(del_prob);
            TG.verifyParents();

            // std::cout << "deleted " << std::endl;
            // for(auto& pr : deledges)
            //     std::cout << pr.first << " " << pr.second << std::endl;
            // std::cout << std::endl;
            
            auto addedges = TG.generateAddEdges(add_prob);
            TG.verifyParents();

            // std::cout << "added " << std::endl;
            // for(auto& edg : addedges)
            //     std::cout << std::get<0>(edg) << " -- " << std::get<1>(edg) << ": " << std::get<2>(edg) << std::endl;
            // std::cout << std::endl;

            batchInsertEdge(deledges, addedges, clusters, (double) 0, [] (double A, double B){
                return A + B;
            });

            testPathQueryValid(clusters, TG.parents, TG.weights, TG.random_perm_map, graph_size);
        }

        deleteRCtree(clusters);

        

    }

    return 0;
}