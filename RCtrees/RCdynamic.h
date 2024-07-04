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
 */

template<typename T, typename D>
void batchInsertEdge( const parlay::sequence<std::pair<T, T>>& delete_edges, const parlay::sequence<std::tuple<T, T, D>>& add_edges, parlay::sequence<cluster<T, D>>& clusters)
{

    int count = 0;

    // obtain initial modified edges

    do
    {


        count++;
    }while(!(count > 1000));    

    if(count > 1000)
        std::cout << red << "[dynamic] definitely went into an infinite loop" << reset << std::endl; 

    return;
}

#endif