#ifndef RC_TEST_H
#define RC_TEST_H
#include "parlay/sequence.h"
#include "cluster.h"

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
 
 
 template <typename T, typename D>
 void check_parents_children(parlay::sequence<cluster<T,D>>& clusters)
 {
     parlay::parallel_for(0, clusters.size(), [&] (T i) {
         auto cluster_ptr = &clusters[i];
         if(cluster_ptr->parent != nullptr)
         {
             bool in_parent = false;
             for(auto& par_child : cluster_ptr->parent.load()->children)
             {
                 if(par_child == cluster_ptr)
                     in_parent = true;
             }
             if(in_parent == false)
             {
                 std::cout << red << "RC tree parents inconsistent!" << reset << std::endl;
                 cluster_ptr->print();
                 cluster_ptr->parent.load()->print();
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
 #endif