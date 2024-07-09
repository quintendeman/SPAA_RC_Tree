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

    auto temp_aff_ptrs = parlay::tabulate(delete_edges.size() + add_edges.size(), [&] (T i) {
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
        parlay::sequence<cluster<T,D>*> ret_seq = {&clusters[v], &clusters[w]};
        return ret_seq;
    });

    auto aff_ptrs = parlay::flatten(temp_aff_ptrs);

    parlay::parallel_for(0, aff_ptrs.size(), [&] (T ii) {
        auto& ptr = aff_ptrs[ii];
        ptr->counter = ii;
    });

    
    parlay::parallel_for(0, aff_ptrs.size(), [&] (T ii) {
        auto& ptr = aff_ptrs[ii];
        if(ptr->counter != ii)
            ptr = nullptr;
    });


    // might as well use this oppurtunity to copy the adjacency list as well
    aff_ptrs = parlay::filter(aff_ptrs, [&] (auto ptr) {
        if(ptr == nullptr)
            return false;
        ptr->counter = 0;
        ptr->state = live | affected | adjacency_changed;
        for(auto i = 0; i < max_neighbours+1; i++)
            ptr->alternate_adjacencies.push_back(ptr->adjacencies[i]);
        return true;
    });

    parlay::sequence<cluster<T,D>*> all_ptrs = parlay::sequence<cluster<T,D>*>(aff_ptrs.size() * max_neighbours, nullptr);

    // revert affected ptrs to original state before clustering
    parlay::parallel_for(0, aff_ptrs.size(), [&] (T ii) {
        const auto& ptr = aff_ptrs[ii];
        for(short i = 1; i < max_neighbours+1; i++)
        {
            const auto& should_be_index = ptr->adjacencies[i];
            auto index_in_ptr = (i-1) * 2;
            if((ptr->ptrs[index_in_ptr] == nullptr && should_be_index == -1) || (ptr->ptrs[index_in_ptr]->index == should_be_index))
            {
                continue;
            }
            else
            {
                auto offset_in_all_ptrs = ii * max_neighbours + (i-1);
                all_ptrs[offset_in_all_ptrs] = getEdge(ptr->index, should_be_index, clusters);
            }
        }
    });
    // revert affected ptrs to original state 
    parlay::parallel_for(0, aff_ptrs.size(), [&] (T ii) {
        const auto& ptr = aff_ptrs[ii];
        for(short i = 1; i < max_neighbours+1; i++)
        {
            const auto& should_be_index = ptr->adjacencies[i];
            auto index_in_ptr = (i-1) * 2;
            if((ptr->ptrs[index_in_ptr] == nullptr && should_be_index == -1) || (ptr->ptrs[index_in_ptr]->index == should_be_index))
            {
                continue;
            }
            else
            {
                auto offset_in_all_ptrs = ii * max_neighbours + (i-1);
                ptr->ptrs[index_in_ptr] = &clusters[should_be_index];
                ptr->ptrs[index_in_ptr + 1] = all_ptrs[offset_in_all_ptrs];
                ptr->types[index_in_ptr] = neighbour_type;
                ptr->types[index_in_ptr+1] = edge_type;
            }
        }
    });

    if(clusters.size() <= 300)
        printTree(clusters);


    parlay::parallel_for(0, delete_edges.size(), [&] (T i) {
        const auto& v = delete_edges[i].first;
        const auto& w = delete_edges[i].second;
        auto edgePtr = getEdge(v, w, clusters);
        if(edgePtr == nullptr)
        {
            std::cout << red << "Empty edge" << reset << std::endl;
            exit(1);
        }
        // std::cout << "delete " << v << " -- " << w << std::endl;
        cluster_allocator::destroy(edgePtr);

        // remove v from alternate adjacency AND ptrs/types
        for(short i = 1; i < max_neighbours+1; i++)
        {
            if(clusters[v].alternate_adjacencies[i] == w)
            {
                clusters[v].alternate_adjacencies[i] = -1;
            }
            if (clusters[w].alternate_adjacencies[i] == v)
            {
                clusters[w].alternate_adjacencies[i] = -1;            
            }
        }
        for(short i = 0; i < clusters[v].size; i+=2)
        {
            if(clusters[v].ptrs[i] != nullptr && clusters[v].ptrs[i]->index == w)
            {
                clusters[v].types[i] = clusters[v].types[i+1] = deleted_type;
            }
            if(clusters[w].ptrs[i] != nullptr && clusters[w].ptrs[i]->index == v)
            {
                clusters[w].types[i] = clusters[w].types[i+1] = deleted_type;
            }
        }
    });

    parlay::parallel_for(0, aff_ptrs.size(), [&] (T i) {
        auto& cluster_ptr = aff_ptrs[i];
        for(short i = 0; i < cluster_ptr->size; i++)
        {
            auto& typ = cluster_ptr->types[i];
            if(typ == deleted_type)
            {
                cluster_ptr->ptrs[i] = nullptr;
            }
        }
    });

    if(clusters.size() <= 300)
        printTree(clusters);


    // ADDING -

    for(short cnt = 0; cnt < max_neighbours; cnt++)
    {
        parlay::parallel_for(0, add_edges.size(), [&] (T i) {
            auto v = std::get<0>(add_edges[i]);
            auto w = std::get<1>(add_edges[i]);
            if(w < v)
                std::swap(v, w);
            for(short i = 1; i < max_neighbours+1; i++)
            {
                if(clusters[v].alternate_adjacencies[i] == w)
                {
                    return;
                }
            }
            clusters[v].counter = i+1;
        });

        parlay::parallel_for(0, add_edges.size(), [&] (T ii){
            auto v = std::get<0>(add_edges[ii]);
            auto w = std::get<1>(add_edges[ii]);
            if(w < v)
                std::swap(v, w);
            if(clusters[v].counter == ii+1)
            {
                auto new_edge = cluster_allocator::alloc();
                new_edge->index = -1;
                new_edge->data = std::get<2>(add_edges[ii]);
                new_edge->state = live | base_edge;
                // std::cout << "Added " << v << " -- " << w << std::endl;
                for(auto i = 0; i < max_neighbours + 1; i++)
                {
                    if(i > 0 && clusters[v].alternate_adjacencies[i] == -1)
                    {
                        clusters[v].alternate_adjacencies[i] = w;
                        clusters[v].add_neighbour(&clusters[w], new_edge);
                        break;
                    }
                }

                clusters[v].counter = 0;
            }
        });
    }

    for(short cnt = 0; cnt < max_neighbours; cnt++)
    {
        parlay::parallel_for(0, add_edges.size(), [&] (T i) {
            auto v = std::get<0>(add_edges[i]);
            auto w = std::get<1>(add_edges[i]);
            if(w < v)
                std::swap(v, w);
            
            for(short i = 1; i < max_neighbours+1; i++)
            {
                if(clusters[w].alternate_adjacencies[i] == v)
                {
                    return;
                }
            }
            clusters[w].counter = i+1;
        });
        parlay::parallel_for(0, add_edges.size(), [&] (T ii){
            auto v = std::get<0>(add_edges[ii]);
            auto w = std::get<1>(add_edges[ii]);
            if(w < v)
                std::swap(v, w);
            if(clusters[w].counter == ii+1)
            {
                // find other edge
                cluster<T,D>* edge_ptr = nullptr;
                for(uint i = 0; i < clusters[v].size; i+=2)
                {
                    auto& ptr = clusters[v].ptrs[i];
                    if(ptr == &clusters[w])
                        edge_ptr = clusters[v].ptrs[i+1];
                }
                if(edge_ptr == nullptr)
                {
    
                    std::cout << red << "NULLPTR " << v << " -- " << w << reset << std::endl;
                    clusters[v].print();
                    clusters[w].print();

                    exit(1);
                }

                for(short i = 0; i < max_neighbours + 1; i++)
                {
                    if(i > 0 && clusters[w].alternate_adjacencies[i] == -1)
                    {
                        clusters[w].alternate_adjacencies[i] = v;                
                        clusters[w].add_neighbour(&clusters[v], edge_ptr);
                    }
                }
                clusters[w].counter = 0;
            }
        });
    }

    if(clusters.size() <= 300)
        printTree(clusters);

    // new loop
    // write it in w's list
    do
    {


        count++;
    }while(!(count > 1000));    

    if(count > 1000)
        std::cout << red << "[dynamic] definitely went into an infinite loop" << reset << std::endl; 

    return;
}

#endif