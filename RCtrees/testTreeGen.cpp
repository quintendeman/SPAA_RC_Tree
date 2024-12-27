#include <iostream>
#include "treeGen.h"

#include "RC.h"
#include "RCdynamic.h"
#include <ctime>    // For time()


using vertex = long;
using utils = graph_utils<vertex>;
using graph = parlay::sequence<parlay::sequence<vertex> >;

int main()
{
    std::srand(std::time(0));

    double min_weight = 0.0;
    double max_weight = 100.0f;
    vertex graph_size;
    graph_size = 1000000 + (rand() % 10000000l);
    // // graph_size+=1000;
    // graph_size = 100 + (rand() % 100);
    // graph_size = 30;

    std::cout << "Working with graph_size " << graph_size << std::endl;
    TreeGen<vertex, double> TG(graph_size, min_weight, max_weight, 0.1, 40, exponential);

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
    

    graph G;
    G = convert_parents_to_graph(G, TG.parents);
    
    parlay::sequence<cluster<vertex, double>> clusters;
    
    const vertex max_degree = 3;

    create_base_clusters(G, clusters, max_degree);
    create_RC_tree(clusters, graph_size, false);

    TG.clearWeights(0.0);



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

    for(unsigned int i = 0; i < 10; i++)
    {
        double del_prob = (rand() % 100);
        del_prob/=100;
        double add_prob = (rand() % 100);
        add_prob/=100;
        
        auto deledges = TG.generateDeleteEdges(del_prob);
        TG.verifyParents();
        
        auto addedges = TG.generateAddEdges(add_prob);
        TG.verifyParents();

        batchInsertEdge(deledges, addedges, clusters, (double) 0, [] (double A, double B){
            return A + B;
        });

        testPathQueryValid(clusters, TG.parents, TG.weights);
    }

    deleteRCtree(clusters);

    return 0;
}