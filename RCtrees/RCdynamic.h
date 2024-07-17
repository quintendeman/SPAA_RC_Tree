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
    // auto tree_nodes = parlay::flatten(parlay::delayed_tabulate(delete_edges.size() + add_edges.size(), [&] (T i) {
    //     T v, w;
    //     if(i < delete_edges.size())
    //     {
    //         v = delete_edges[i].first;
    //         w = delete_edges[i].second;
    //     }
    //     else
    //     {
    //         v = std::get<0>(add_edges[i-delete_edges.size()]);
    //         w = std::get<1>(add_edges[i-delete_edges.size()]);
    //     }
    //     parlay::sequence<node<T,D>*> ret_seq = {clusters[v].adjacency.get_head(), clusters[w].adjacency.get_head()};
    //     ret_seq[0]->tiebreak = i*2;
    //     ret_seq[1]->tiebreak = i*2 +1;
    //     return ret_seq;
    // }));

    // parlay::parallel_for(0, tree_nodes.size(), [&] (T i){
    //     if(ret_seq[i]->tiebreak != i)
    //         ret_seq[i] = nullptr;
    // });

    // ret_seq = parlay::filter(ret_seq, [] (auto node_ptr){
    //     if(node_ptr != nullptr)
    //         node_ptr->state |=
    //     return node_ptr != nullptr;
    // });


    return;
}

#endif