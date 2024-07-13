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

template<typename T, typename D>
void contract(cluster<T,D>* cluster_ptr)
{
    if(cluster_ptr->get_neighbour_count() == 0)
    {
        cluster_ptr->state&=(~live);
        cluster_ptr->state|=(nullary_cluster);
        cluster_ptr->state|=internal;
    }
    else if(cluster_ptr->get_neighbour_count() == 1)
    {
        cluster<T,D>* edge_ptr = nullptr;
        cluster<T,D>* other_side = nullptr;
        short neighbour_index = -1;
        
        for(short i = 0; i < cluster_ptr->size; i+=2)
        {
            if(cluster_ptr->types[i] & neighbour_type)
            {
                other_side = cluster_ptr->ptrs[i];
                edge_ptr = cluster_ptr->ptrs[i+1];
                neighbour_index = i;
            }
        }

        // now add these two as children and remove as neighbours if neighbours
        cluster_ptr->change_to_child(edge_ptr);
        other_side->change_to_child(cluster_ptr);

        edge_ptr->set_parent(cluster_ptr);
        cluster_ptr->set_parent(other_side);

        // // mark both of these as not live
        edge_ptr->state&=(~live);
        cluster_ptr->state&=(~live);

        cluster_ptr->state|=unary_cluster;
        cluster_ptr->state|=internal;
    }
    else 
    if (cluster_ptr->get_neighbour_count() == 2)
    {
        // find left and right vertices/nodes
        cluster<T,D>* left_edge_ptr = nullptr;
        cluster<T,D>* right_edge_ptr = nullptr;
        cluster<T,D>* left_node_ptr = nullptr;
        cluster<T,D>* right_node_ptr = nullptr;

        cluster_ptr->get_two_neighbours_edges(left_node_ptr, left_edge_ptr, right_node_ptr, right_edge_ptr);

        left_node_ptr->overwrite_neighbour(cluster_ptr, right_node_ptr, cluster_ptr);

        right_node_ptr->overwrite_neighbour(cluster_ptr, left_node_ptr, cluster_ptr);
        
        left_edge_ptr->set_parent(cluster_ptr);
        right_edge_ptr->set_parent(cluster_ptr);

        cluster_ptr->change_to_child(left_edge_ptr);
        cluster_ptr->change_to_child(right_edge_ptr);

        
        cluster_ptr->state&=(~live);

        cluster_ptr->state|=binary_cluster;
        cluster_ptr->state|=internal;
    }

}

/*
    Get two hop neighbours at that level
*/
template<typename T, typename D>
parlay::sequence<cluster<T,D>*> get_3dp1(const parlay::sequence<cluster<T,D>*>& aff_ptrs, parlay::sequence<cluster<T,D>>& clusters, const unsigned char level)
{
    parlay::sequence<cluster<T,D>*> frontier;
    
    frontier = parlay::flatten(
        parlay::delayed_tabulate(aff_ptrs.size(), [&] (T i) {
            auto& ptr = aff_ptrs[i];
            auto ret_seq = parlay::sequence<cluster<T,D>*>(max_neighbours, nullptr);

            if(ptr->alternate_adjacency.size() > 0)
            {
                auto& w_node = *ptr->alternate_adjacency[level];
                for(short i = 0; i < w_node.size(); i++)
                {
                    auto& w = w_node[i];
                    if(w != -1)
                        ret_seq[i] = &clusters[w];
                }
            }
            else
            {
                auto& w_node = *ptr->adjacency[level];
                for(short i = 0; i < w_node.size(); i++)
                {
                    auto& w = w_node[i];
                    if(w != -1)
                        ret_seq[i] = &clusters[w];
                }
            }
            return ret_seq;
        }));

    // TODO merge?
    frontier = parlay::flatten(
        parlay::delayed_tabulate(frontier.size(), [&] (T i) {
            auto& ptr = frontier[i];
            auto ret_seq = parlay::sequence<cluster<T,D>*>(max_neighbours + 1, nullptr);

            if(ptr == nullptr)
                return ret_seq;

            if(ptr->alternate_adjacency.size() > 0)
            {
                auto& w_node = *ptr->alternate_adjacency[level];
                for(short i = 0; i < w_node.size(); i++)
                {
                    auto& w = w_node[i];
                    if(w != -1)
                        ret_seq[i] = &clusters[w];
                }
            }
            else
            {
                auto& w_node = *ptr->adjacency[level];
                for(short i = 0; i < w_node.size(); i++)
                {
                    auto& w = w_node[i];
                    if(w != -1)
                        ret_seq[i] = &clusters[w];
                }
            }
            ret_seq[max_neighbours] = ptr;
            return ret_seq;
        }));
    
    parlay::parallel_for(0, frontier.size(), [&] (T i) {
        auto& ptr = frontier[i];
        if(ptr != nullptr)
            ptr->counter = i+1;
    });
    parlay::parallel_for(0, frontier.size(), [&] (T i) {
        auto& ptr = frontier[i];
        if(ptr != nullptr && ptr->counter != i+1)
            ptr = nullptr;
    });

    parlay::filter(frontier, [] (auto ptr){
        if(ptr != nullptr)
        {
            ptr->counter = 0;
            return true;
        }
        return false;
    });
    
    return frontier;
}

template<typename T, typename D>
bool first_condition(cluster<T,D>* ptr)
{
    
    return false;
}

template<typename T, typename D>
bool second_condition(cluster<T,D>* ptr)
{
    return false;
}

template<typename T, typename D>
bool third_condition(cluster<T,D>* ptr)
{
    return false;
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

    // TODO MERGE
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
                ptr->types[i] = neighbour_type;
                ptr->types[i+1] = edge_type;
            }
            else if (other_ptr != nullptr)
            {
                ptr->types[i] = neighbour_type;
                ptr->types[i+1] = edge_type;
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
            
                auto& w_arr = *(clusters[v].alternate_adjacency[0]);
                for(short ii = 0; ii < w_arr.size(); ii++)
                {
                    if(w_arr[ii] == -1)
                    {
                        w_arr[ii] = w;
                        break;
                    }
                }
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
                auto& v_arr = *(clusters[w].alternate_adjacency[0]);
                for(short ii = 0; ii < v_arr.size(); ii++)
                {
                    if(v_arr[ii] == -1)
                    {
                        v_arr[ii] = v;
                        break;
                    }
                }
                clusters[w].counter = 0;
            }
        });

        add_edges_wflags = parlay::filter(add_edges_wflags, [] (auto edge) {return edge.valid;});
    }while(add_edges_wflags.size() > 0);

    do
    {
        // Find MIS of update-eligible affected vertices in F/F'_{i-1}
        set_MIS(aff_ptrs);
        auto mis_set = parlay::filter(aff_ptrs, [] (auto ptr) {
            return ptr->state & IS_MIS_SET;;
        });

        parlay::parallel_for(0, mis_set.size(), [&] (T i) {
            auto& ptr = mis_set[i];
            contract(ptr);
            ptr->adjacency.adopt(&(ptr->alternate_adjacency)); // already contracted, we don't need to keep track of the adjacents
            ptr->state &= (~affected);
        });
        
        parlay::parallel_for(0, aff_ptrs.size(), [&] (T i){
            auto& ptr = aff_ptrs[i];
            if(ptr->state & affected)
                ptr->add_alternate_level(ptr->state, count);
        });

        aff_ptrs = get_3dp1(aff_ptrs, clusters, count+1);

        aff_ptrs = parlay::filter(aff_ptrs, [] (auto ptr) {
            return first_condition(ptr) || second_condition(ptr) || third_condition(ptr);
        });

        // set alternate adjacencies for new aff_ptrs?
        parlay::parallel_for(0, aff_ptrs.size(), [&] (T i) {
            auto& ptr = aff_ptrs[i];
            if(ptr->alternate_adjacency.size() == 0)
            {
                ptr->alternate_adjacency.copy_till_level(&ptr->adjacency, count+1);
            }
            for(short i = 0; i < ptr->size; i+=2)
            {
                auto& w_node = *ptr->alternate_adjacency.get_tail();
                auto& w = w_node[i/2];
                auto& other_ptr = ptr->ptrs[i];
                // get ptrs to all the edges?
                if(other_ptr != nullptr && w != -1 && other_ptr->index != w)
                {
                    ptr->ptrs[i+1] = getEdge(w, ptr->index, clusters);
                    if(ptr->ptrs[i+1] == nullptr) // TODO remove
                    {
                        std::cout << red << "nullptr" << reset << std::endl;
                        exit(1);
                    }
                    ptr->ptrs[i] = &clusters[w];
                    ptr->types[i] = neighbour_type;
                }
            }
        });

        // ensure ptrs haven't already contracted this session? done because we only care about live ones and the adjacency list is short
        count++;
    }while(!(count > 100) && aff_ptrs.size() > 0);    

    if(count > 100)
        std::cout << red << "[dynamic] definitely went into an infinite loop" << reset << std::endl; 

    return;
}

#endif