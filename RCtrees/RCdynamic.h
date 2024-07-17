#ifndef RCDYNAMIC_H
#define RCDYNAMIC_H

#include <atomic>
#include <random>
#include <algorithm>
#include <iostream>
#include <mutex>
#include "cluster.h"
#include "RC.h"
#include <parlay/alloc.h>

template <typename edgetype>
struct edge_with_flag
{
    edgetype E;
    bool valid = true;
};


template<typename T, typename D>
bool first_condition(node<T,D>* node)
{
    return node->state & affected;
}

template<typename T, typename D>
bool second_condition(node<T,D>* node)
{
    return node->state & adjacency_changed;
}

template<typename T, typename D>
bool this_condition(node<T,D>* node)
{
    return node->state & all_contracted_affected;
}




template<typename T, typename D>
void batchInsertEdge( const parlay::sequence<std::pair<T, T>>& delete_edges, const parlay::sequence<std::tuple<T, T, D>>& add_edges, parlay::sequence<cluster<T, D>>& clusters)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    auto tree_nodes = parlay::flatten(parlay::tabulate(delete_edges.size() + add_edges.size(), [&] (T i) {
        T v, w;
        if(i < delete_edges.size())
        {
            v = delete_edges[i].first;
            w = delete_edges[i].second;
        }
        else
        {
            v = std::get<0>(add_edges[i-delete_edges.size()]);
            w = std::get<1>(add_edges[i-delete_edges.size()]);
        }
        parlay::sequence<node<T,D>*> ret_seq = {clusters[v].adjacency.get_head(), clusters[w].adjacency.get_head()};
        if(ret_seq[0] == nullptr || ret_seq[1] == nullptr) // TODO remove
        {
            std::cout << red << "This should never happen, nullptr on head?" << reset << std::endl;
            exit(1);
        }
        ret_seq[0]->tiebreak() = i*2;
        ret_seq[1]->tiebreak() = i*2 +1;
        return ret_seq;
    }));

    parlay::parallel_for(0, tree_nodes.size(), [&] (T i){
        if(tree_nodes[i]->tiebreak() != i)
            tree_nodes[i] = nullptr;
    });

    tree_nodes = parlay::filter(tree_nodes, [] (auto node_ptr){
        if(node_ptr != nullptr)
        {
            node_ptr->state |= affected | adjacency_changed | live;
            node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
            while(node_ptr->cluster_ptr->adjacency.size() > 1)
                node_ptr->cluster_ptr->adjacency.delete_tail();
        }
        return node_ptr != nullptr;
    });

    parlay::parallel_for(0, delete_edges.size(), [&] (T I){
        auto& v = delete_edges[I].first;
        auto& w = delete_edges[I].second;
        
        auto v_node_ptr = clusters[v].adjacency.get_head();
        auto w_node_ptr = clusters[w].adjacency.get_head();

        
        node<T,D>* edge_node_ptr = nullptr;
        for(auto& potential_edge_ptr : v_node_ptr->adjacents)
            if(potential_edge_ptr != nullptr && get_other_side(v_node_ptr, potential_edge_ptr) == w_node_ptr)
                edge_node_ptr = potential_edge_ptr;

        for(short i = 0; i < max_neighbours; i++)
        {
            if(v_node_ptr->adjacents[i] == edge_node_ptr)
            {
                v_node_ptr->adjacents[i] = nullptr;
            }
            if(w_node_ptr->adjacents[i] == edge_node_ptr)
            {
                w_node_ptr->adjacents[i] = nullptr;
            }
        }
        cluster_allocator::destroy(edge_node_ptr->cluster_ptr);
    });

    // auto add_edges_wflags = parlay::tabulate(add_edges.size(), [&] (T i) {
    //     edge_with_flag<std::tuple<T, T, D>> ret_edge;
    //     ret_edge.valid = true;
    //     ret_edge.E = add_edges[i];
    //     return std::move(ret_edge);
    // });



    return;
}

#endif