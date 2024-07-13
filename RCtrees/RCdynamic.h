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

/**
 * Debug function, ideally we never have to use this
 */
template<typename T, typename D>
std::array<T, max_neighbours> getAdjacencyAtLevel(const T& cluster_idx, short contraction_level, const parlay::sequence<cluster<T, D>> &clusters)
{
    const cluster<T, D>* cluster_ptr = &clusters[cluster_idx];
    const short& current_level = cluster_ptr->contraction_time;
    if(contraction_level > current_level)
        contraction_level = current_level;

    std::array<T, max_neighbours> adjacency = cluster_ptr->initial_adjacency;

    for(short i = 0; i <= contraction_level; i++)
    {
        for(auto& v: adjacency)
        {
            if (v == -1)
                continue;
            
            while(cluster_ptr->isPtr(v) == -1)
            {
                // find the relevant edge
                auto relevant_node = &clusters[v];

                if(relevant_node->contraction_time > i)
                    break;
                // find other side of this -- must have been a compress
                cluster<T, D>* other_side = nullptr;
                if (cluster_ptr == relevant_node->get_parent())
                {
                    other_side = relevant_node->get_one_neighbour();
                }
                else
                {
                    other_side = relevant_node->get_parent();
                }
                v = other_side->index;
            }
        }
        bool finalized = true;
        for(const auto& v : adjacency)
            if(v != -1 && cluster_ptr->isPtr(v) == -1)
                finalized = false;
        if(finalized)
            break;
    }
    return adjacency;
}

template<typename T, typename D>
std::array<T, max_neighbours> getAdjacencyAtLevel(const T& cluster_idx, short contraction_level, std::array<T, max_neighbours> adjacency,  short initial_level, const parlay::sequence<cluster<T, D>> &clusters)
{
    const cluster<T, D>* cluster_ptr = &clusters[cluster_idx];
    const short& current_level = cluster_ptr->contraction_time;
    if(contraction_level > current_level)
        contraction_level = current_level;
    if(initial_level > contraction_level)
        initial_level = contraction_level;

    for(short i = initial_level + 1; i <= contraction_level; i++)
    {
        for(auto& v: adjacency)
        {
            if (v == -1)
                continue;
            
            while(cluster_ptr->isPtr(v) == -1)
            {
                // find the relevant edge
                auto relevant_node = &clusters[v];

                if(relevant_node->contraction_time > i)
                    break;
                // find other side of this -- must have been a compress
                cluster<T, D>* other_side = nullptr;
                if (cluster_ptr == relevant_node->get_parent())
                {
                    other_side = relevant_node->get_one_neighbour();
                }
                else
                {
                    other_side = relevant_node->get_parent();
                }
                v = other_side->index;
            }
        }
        bool finalized = true;
        for(const auto& v : adjacency)
            if(v != -1 && cluster_ptr->isPtr(v) == -1)
                finalized = false;
        if(finalized)
            break;
    }
    return adjacency;
}


/**
 * Assumes that "clusters" contains a bunch of ALREADY fully contracted trees i.e. a fully contracted
 * 
 * There must be no repeated edges in delete_edges
 * 
 * or in add_edges -- a repeated edge between delete_edges and add_edges is fine
 */

template<typename T, typename D>
void batchInsertEdge( const parlay::sequence<std::pair<T, T>>& delete_edges, const parlay::sequence<std::tuple<T, T, D>>& add_edges, parlay::sequence<cluster<T, D>>& clusters)
{

    int count = 0;
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    auto all_ptrs_pairs = parlay::delayed_tabulate(delete_edges.size() + add_edges.size(), [&] (T i) {
        T v, w;
        if(i < delete_edges.size())
        {
            v = delete_edges[i].first;
            w = delete_edges[i].second;
        }
        else
        {
            v = std::get<0>(add_edges[i - delete_edges.size()]);
            w = std::get<1>(add_edges[i - delete_edges.size()]);
        }
        
        parlay::short_sequence<cluster<T,D>*> ret_seq = {&clusters[v], &clusters[w]};
        clusters[v].counter = i*2;
        clusters[w].counter = i*2 + 1;
        return ret_seq;
    });
    auto all_ptrs = parlay::flatten(all_ptrs_pairs);

    parlay::parallel_for(0, all_ptrs.size(), [&] (T i) {
        if(all_ptrs[i]->counter != i)
            all_ptrs[i] = nullptr;
    });

    auto aff_ptrs = parlay::filter(all_ptrs, [] (const auto& ptr){
        if(ptr != nullptr)
        {
            ptr->counter = 0;
            ptr->state = affected | live;
        }
        return ptr != nullptr;
    });

    parlay::parallel_for(0, delete_edges.size(), [&] (T i) {
        const auto& edge = delete_edges[i];
        const auto& v = edge.first;
        const auto& w = edge.second;
        
        auto edgePtr = getEdge(v, w, clusters);
        if(edgePtr == nullptr) // TODO: REMOVE
        {
            std::cout << red << "Something went wrong, got nullptr" << reset << std::endl;
            exit(1);
        }
        cluster_allocator::destroy(edgePtr);
        
        auto& v_arr = *clusters[v].adjacency[0];
        for(short i = 0; i < clusters[v].size; i+=2)
        {
            if(v_arr[i/2] == w)
                clusters[v].types[i] = clusters[v].types[i+1] = deleted_type;
        }
        auto& w_arr = *clusters[w].adjacency[0];
        for(short i = 0; i < clusters[w].size; i+=2)
        {
            if(w_arr[i/2] == v)
                clusters[w].types[i] = clusters[w].types[i+1] = deleted_type;
        }
        return;
    });

    parlay::parallel_for(0, aff_ptrs.size(), [&] (T i) {
        auto& ptr = aff_ptrs[i];
        const auto& v = ptr->index;
        for(short i = 0; i < ptr->size; i++)
            if(ptr->types[i] == deleted_type)
            {
                ptr->types[i] = empty_type;
                ptr->ptrs[i] = nullptr;
            }
        auto& w_arr = *(ptr->adjacency[0]);
        for(short i = 0; i < ptr->size; i+=2)
        {
            auto& other_ptr = ptr->ptrs[i];
            auto& w = w_arr[i/2];
            if(other_ptr != nullptr && w != -1 && w != other_ptr->index)
            {
                ptr->ptrs[i] = &clusters[w];
                ptr->ptrs[i+1] = getEdge(v, w, clusters);
            }
        }
        ptr->counter = 0;
        ptr->add_alternate_level(affected | live, 0);
    });

    auto add_edges_wflags = parlay::tabulate(add_edges.size(), [&] (T i) {
        edge_with_flag<std::tuple<T, T, D>> ret_edge;
        ret_edge.valid = true;
        ret_edge.E = add_edges[i];
        return std::move(ret_edge);
    });

    auto reduce_edges = add_edges_wflags;
    
    bool first_time = true;
    do
    {
        parlay::parallel_for(0, reduce_edges.size(), [&] (T i) {
            auto& edge = reduce_edges[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);
            clusters[v].counter = i + 1;
        });
        parlay::parallel_for(0, reduce_edges.size(), [&] (T i) {
            auto& edge = reduce_edges[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);

            
            if(clusters[v].counter == (i+1))
            {
                reduce_edges[i].valid = false;

                auto newEdgePtr = cluster_allocator::alloc();
                newEdgePtr->state = live | base_edge;
                newEdgePtr->data = std::get<2>(edge.E);
                newEdgePtr->add_initial_neighbours(&clusters[v], &clusters[w]);
                clusters[v].add_neighbour(&clusters[w], newEdgePtr);
            
                // auto& w_arr = *(clusters[v].adjacency[0]);
                // for(short ii = 0; ii < w_arr.size(); ii++)
                // {
                //     if(w_arr[ii] == -1)
                //     {
                //         w_arr[ii] = w;
                //         break;
                //     }
                // }
                clusters[v].counter = 0;

            }

        });

        reduce_edges = parlay::filter(reduce_edges, [] (auto r) {return r.valid;});
    }while(reduce_edges.size() > 0);

    do
    {
        parlay::parallel_for(0, add_edges_wflags.size(), [&] (T i) {
            auto& edge = add_edges_wflags[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);
            clusters[w].counter = i + 1;
        });

        parlay::parallel_for(0, add_edges_wflags.size(), [&] (T i) {
            auto& edge = add_edges_wflags[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);
            if(clusters[w].counter == (i+1))
            {
                add_edges_wflags[i].valid = false;

                auto new_edge = getEdge(v, w, clusters);
                if(new_edge == nullptr) // TODO: REMOVE
                {
                    std::cout << red << "nullptr??" << reset << std::endl;
                    clusters[v].print();clusters[w].print();
                    exit(1);
                }
                clusters[w].add_neighbour(&clusters[v], new_edge);
                clusters[w].counter = 0;
            }
        });

        add_edges_wflags = parlay::filter(add_edges_wflags, [] (auto edge) {return edge.valid;});
    }while(add_edges_wflags.size() > 0);


    do
    {


        count++;
    }while(!(count > 1000));    

    if(count > 1000)
        std::cout << red << "[dynamic] definitely went into an infinite loop" << reset << std::endl; 

    return;
}

#endif