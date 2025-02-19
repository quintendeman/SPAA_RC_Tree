//testing functions for RC trees
//Pulling from RC.cpp because LCA code will use too 

#ifndef RCTEST_H
#define RCTEST_H


#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "cluster.h"
#include "subtree_query.h"


template <typename T, typename D>
void check_parents_children(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        auto cluster_ptr = &clusters[i];
        if(cluster_ptr->parent != nullptr)
        {
            bool in_parent = false;
            for(auto& par_child : cluster_ptr->parent->children)
            {
                if(par_child == cluster_ptr)
                    in_parent = true;
            }
            if(in_parent == false)
            {
                std::cout << red << "RC tree parents inconsistent!" << reset << std::endl;
                cluster_ptr->print();
                cluster_ptr->parent->print();
                exit(1);
            }
        }
        for(auto& child : cluster_ptr->children)
        {
            if(child == nullptr)
                continue;
            if(child->parent != cluster_ptr)
            {
                std::cout << red << "RC tree children inconsistent!" << reset << std::endl;
                cluster_ptr->print();
                child->print();
                exit(1);
            }
        }
    });
}

template<typename T, typename D>
void check_children_values(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        auto cluster_ptr = &clusters[i];
        D value_from_children = 0.0;
        for(auto& child : cluster_ptr->children)
        {
            if(child == nullptr)
                continue;
            value_from_children+=child->data;
        }
        if(value_from_children != cluster_ptr->data)
        {
            std::cout << red << "RC tree values inconsistent!" << reset << std::endl;
            cluster_ptr->print();
            exit(1);
        }
    });
}


//given an RC tree, makes sure that various invariants hold
//the parent of a child of v is v
//the child of a parent of v is v
//for addition) makes sure that the value of a parent is the sum of the values of the children
//specific type of subtree query works
//if print is true, print out debugging info, otherwise don't
template<typename vertex, typename datatype> 
void test_rc_valid(const parlay::sequence<vertex>& parents, parlay::sequence<cluster<vertex, datatype>>& clusters, bool print=true)
{

    check_parents_children(clusters);
    check_children_values(clusters);

    parlay::sequence<std::tuple<vertex, vertex, datatype>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, datatype>(i, parents[i], (datatype) ( 1.0));
    });
    //only keep nondegenerate edges (where the two endpoints are different)
    weighted_edges = parlay::filter(weighted_edges, [&] (auto edge) {
        return std::get<0>(edge) != std::get<1>(edge);
    });
    
    //change all of the weights in the tree to 1
    batchModifyEdgeWeights(weighted_edges, [] (datatype a, datatype b) {
        return a + b;
    }, clusters);

    auto random_index = rand() % parents.size();
    if (print)
    std::cout << "Root: " << random_index << " Parent: " << parents[random_index] << std::endl;

    auto manual_sum_val = manual_subtree_sum(&clusters[random_index], &clusters[parents[random_index]], clusters);
    if (print)
    std::cout << "Manual subtree sum: " << red << manual_sum_val << reset << std::endl;

    
    auto sub_ret_val = subtree_query(&clusters[random_index], &clusters[parents[random_index]], (datatype) 0.0, [] (datatype a, datatype b) {
        return a+b;
    });

    if (!isNearlyEqual(manual_sum_val,sub_ret_val)) {
        std::cout << "manual and RC tree subtree sums different, aborting " << manual_sum_val << " " << sub_ret_val << std::endl;
        exit(17);
    }
    if (print) {
        std::cout << "Subtree query returns: " << bold << red << sub_ret_val << reset << std::endl;
        std::cout << std::endl;

    }
   

    return;
}


/**
 * Make sure the base clusters are consistent i.e. that every neighbour of mine has me as their neighbour too
 */
template <typename T, typename D>
void check_consistency(parlay::sequence<node<T, D>*>& tree_nodes)
{
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        if((node_ptr->state & live) == 0)
        {
            return;
        }
        if(node_ptr->state & (unary_cluster | binary_cluster | base_edge))
        {
            return;
        }

        for(auto& ptr : node_ptr->adjacents)
        {
            if(ptr == nullptr || !(ptr->state & (base_edge | binary_cluster)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, ptr);
            if(get_other_side(other_node_ptr, ptr) != node_ptr)
            {
                node_ptr->cluster_ptr->print();
                ptr->cluster_ptr->print();
                other_node_ptr->cluster_ptr->print();
                std::cout << "Base clusters inconsistent" << std::endl;
                exit(1);
            }
        }
    });

    return;
}



template<typename T, typename D>
void check_counter(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        if(clusters[i].counter)
        {
            std::cout << red << "Counter wasn't reset" << reset << std::endl;
            exit(1);
        }
    });
}

#endif