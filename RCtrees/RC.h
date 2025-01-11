#ifndef RC_H
#define RC_H

/**
 * Mainly deals with static RC tree creation
*/

#include <atomic>
#include <random>
#include <set>
#include <algorithm>
#include <iostream>
#include <mutex>
#include "../examples/helper/graph_utils.h"
#include <parlay/alloc.h>
#include "../examples/counting_sort.h"
#include <parlay/random.h>

#include "cluster.h"
#include "utils.h" //for bit tricks
#include "MIS.h"
#include "RC_test.h"



template <typename edgetype>
struct edge_with_flag
{
    edgetype E;
    bool valid = true;
};

template<typename T, typename D>
void finalize(node<T,D>* contracted_node)
{
    while(contracted_node->cluster_ptr->adjacency.get_tail() != contracted_node)
    {
        auto tail = contracted_node->cluster_ptr->adjacency.get_tail();
        if(tail->state & dont_finalize)
            break;
        contracted_node->cluster_ptr->adjacency.delete_tail();
    }

    return;
}



/**
 * MUST have at most 2 neighbours
 */
template<typename T, typename D>
void contract(node<T,D>* node_ptr, bool affect = false)
{
    
    if(node_ptr->get_num_neighbours_live() == 0)
    {
        node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        node_ptr->state |= contracts_this_round | nullary_cluster;
        node_ptr->state &= (~live);
        node_ptr->state &= (~affected);
        node_ptr->cluster_ptr->parent = nullptr;
        if(affect)
        {
            for(auto& ptr : node_ptr->adjacents)
                ptr = nullptr;
            finalize(node_ptr);
        }
    }
    else if(node_ptr->get_num_neighbours_live() == 1) // unary
    {
        node<T,D>* neighbour_node;
        node<T,D>* edge_node;
        for(auto& ptr : node_ptr->adjacents)
            if(ptr != nullptr && ptr->state & (base_edge | binary_cluster))
            {
                edge_node = ptr;
                break;
            }
        neighbour_node = get_other_side(node_ptr, edge_node);

        // do nothing for myself
        node_ptr->state &= (~live);
        node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        node_ptr->state |= (unary_cluster | contracts_this_round);

        node_ptr->cluster_ptr->parent = neighbour_node->cluster_ptr;
        for(short k = 0; k < max_neighbours; k++)
        {
            auto& ptr = neighbour_node->adjacents[k];
            if(ptr == edge_node)
            {
                neighbour_node->cluster_ptr->children[k] = node_ptr->cluster_ptr; 
                break;
            }
        }

        edge_node->cluster_ptr->parent = node_ptr->cluster_ptr;
        

        for(short k = 0; k < max_neighbours; k++)
        {
            auto& ptr = node_ptr->adjacents[k];
            if(ptr == edge_node)
            {
                node_ptr->cluster_ptr->children[k] = edge_node->cluster_ptr; 
                break;
            }
        }

        for(auto& ptr : neighbour_node->adjacents)
            if(ptr == edge_node)
                ptr = node_ptr;
        if(affect == true)
        {
            finalize(node_ptr);
            finalize(edge_node);
            // for(auto& ptr : node_ptr->adjacents)
            //     ptr = nullptr;
            neighbour_node->state |= (affected | adjacency_changed | live);
            node_ptr->state &= (~affected);
            node_ptr->state &= (~adjacency_changed);
            neighbour_node->state &= (~(binary_cluster | base_edge));
        }
    }
    else if (node_ptr->get_num_neighbours_live() == 2)
    {
        node<T,D>* left_edge;
        node<T,D>* left_node;
        node<T,D>* right_edge;
        node<T,D>* right_node;
        for(auto& ptr : node_ptr->adjacents)
            if(ptr != nullptr && ptr->state & (base_edge | binary_cluster))
            {
                left_edge = ptr;
                left_node = get_other_side(node_ptr, left_edge);
                break;
            }
        for(auto& ptr : node_ptr->adjacents)
            if(ptr != nullptr && ptr->state & (base_edge | binary_cluster) && ptr != left_edge)
            {
                right_edge = ptr;
                right_node = get_other_side(node_ptr, right_edge);
                break;
            }
        node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        node_ptr->state |= binary_cluster | contracts_this_round;
        node_ptr->state &= (~live);
        node_ptr->cluster_ptr->parent = nullptr;
        

        left_edge->cluster_ptr->parent = node_ptr->cluster_ptr;
        right_edge->cluster_ptr->parent = node_ptr->cluster_ptr;
        
        for (short k = 0; k < node_ptr->adjacents.size(); k++)
        {
            auto& ptr = node_ptr->adjacents[k];
            if (ptr == left_edge)
            {
                ptr = left_node;
                node_ptr->cluster_ptr->children[k] = left_edge->cluster_ptr;
            }
            if (ptr == right_edge)
            {
                node_ptr->cluster_ptr->children[k] = right_edge->cluster_ptr;
                ptr = right_node;
            }
        }

        left_node->set_ptr_from_list(left_edge,node_ptr);
        right_node->set_ptr_from_list(right_edge,node_ptr);


        if(affect == true)
        {
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != left_node && ptr != right_node)
                    ptr = nullptr;
            }
            right_node->state |= (affected | adjacency_changed | live);
            left_node->state |= (affected | adjacency_changed | live);
            right_node->state &= (~(binary_cluster | base_edge));
            left_node->state &= (~(binary_cluster | base_edge));
            node_ptr->state &= (~affected);
            node_ptr->state &= (~adjacency_changed);
            finalize(node_ptr);
            finalize(left_edge);
            finalize(right_edge);
        }
    }
}

/**
 * Given an symmetric graph, creates a set of clusters
 * In total, it creates n + m clusters in the array base_clusters
 * The first n are base_vertex clusters
 * it returns an array of n clusters corresponding to the original n vertices.
 * The edge clusters are present as edges.
 * No weights
*/
template <typename T, typename D>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T, D> > &base_clusters, const T max_size, D defretval = 0.00f)
{

    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    base_clusters = parlay::tabulate(G.size(), [] (T i) {
        cluster<T,D> base_cluster;
        return base_cluster;
    });
    // TODO merge with above!
    parlay::parallel_for(0, base_clusters.size(), [&] (T i){
        base_clusters[i].index = i;
        base_clusters[i].add_empty_level(live, 0);
    });


    parlay::parallel_for(0, base_clusters.size(), [&] (const T i){
        auto& clstr = base_clusters[i];
        clstr.state = base_vertex | live;
        clstr.data = defretval; //must reset base vertex clusters to default val as well as base edge clusters
        
        const auto& v= i;

        for(const auto& w : G[v])
        {
            if(w < v)
                continue;
            auto edge_cluster = cluster_allocator::alloc();
            edge_cluster->index = -1;
            edge_cluster->state = base_edge | live;
            edge_cluster->data = defretval;
            edge_cluster->max_weight_edge = edge_cluster;
            edge_cluster->add_empty_level(live | base_edge, 0);
            edge_cluster->add_ptr_to_highest_level(clstr.adjacency.get_tail());
            edge_cluster->add_ptr_to_highest_level(base_clusters[w].adjacency.get_tail());
            clstr.add_ptr_to_highest_level(edge_cluster->adjacency.get_tail());
        }
    });

    parlay::parallel_for(0, G.size(), [&] (T v) {
        auto cluster_ptr = &base_clusters[v];
        cluster_ptr->index = v;
        for(const auto& w : G[v])
        {
            if(v < w) // reversed
                continue;
            // find the edge that corresponds/joins to the other side
            const auto& node_ptr_arr = base_clusters[w].adjacency.get_tail()->adjacents;

            node<T,D>* edge_ptr = nullptr;
            bool found = false;
            for(auto& ptr : node_ptr_arr)
            {
                if(ptr != nullptr)
                {
                    const auto& final_adjacency_ptr_arr = ptr->adjacents;
                    for(const auto& v_ptr : final_adjacency_ptr_arr)
                    {
                        if(v_ptr != nullptr && v_ptr->cluster_ptr == cluster_ptr)
                        {    
                            edge_ptr = ptr;
                            found = true;
                            break;
                        }
                    }
                }
                if(found)
                    break;
            }
            base_clusters[v].add_ptr_to_highest_level(edge_ptr);
        }
    });

    return;
}

template<typename T, typename D>
void create_base_clusters(parlay::sequence<cluster<T,D>>& base_clusters, parlay::sequence<std::tuple<T,T,D>> initial_edges, const T max_degree, const T num_vertices)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;
    // parlay::sequence<cluster<T,D>> base_clusters 
    // auto& base_clusters = empty_clusters;
    base_clusters = parlay::tabulate(num_vertices, [] (T i) {
        cluster<T,D> base_cluster;
        // base_cluster.add_empty_level(live, 0);
        // return base_cluster;
        return base_cluster;
    });
    // TODO merge with above?
    parlay::parallel_for(0, base_clusters.size(), [&] (T i) {
        base_clusters[i].add_empty_level(live, 0);
        base_clusters[i].state = base_vertex | live;
        base_clusters[i].index = i;
        
    });

    auto edges_wflag = parlay::tabulate(initial_edges.size(), [&] (T i) {
        edge_with_flag<std::tuple<T, T, D>> ret_edge;
        ret_edge.valid = true;
        ret_edge.E = initial_edges[i];
        return std::move(ret_edge);
    });

    while(edges_wflag.size())
    {
        // std::cout << "left side add edge sizes " << edges_wflag.size() << std::endl;
        parlay::parallel_for(0, edges_wflag.size(), [&] (T i) {
            auto& edge = edges_wflag[i].E;
            T v = std::get<0>(edge);
            T w = std::get<1>(edge);
            D weight = std::get<2>(edge);
            assert(v >= 0 && v < num_vertices);
            assert(w >= 0 && w < num_vertices);

            if(v < w)
                std::swap(v, w); // have a consistent view
            base_clusters[v].tiebreak = i + 1; // 
        });

        parlay::parallel_for(0, edges_wflag.size(), [&] (T i) {
            auto& edge = edges_wflag[i].E;
            T v = std::get<0>(edge);
            T w = std::get<1>(edge);
            D weight = std::get<2>(edge);
            if(v < w)
                std::swap(v, w); // have a consistent view
            if(base_clusters[v].tiebreak != i + 1)
                return;
            edges_wflag[i].valid = false;

            // std::string print_string = std::to_string(v) + " -- " + std::to_string(w) + " tiebreak succeeded\n";
            // std::cout << print_string;

            auto edge_cluster = cluster_allocator::alloc();
            edge_cluster->index = -1;
            edge_cluster->state = base_edge | live;
            edge_cluster->data = weight;
            edge_cluster->max_weight_edge = edge_cluster;
            edge_cluster->add_empty_level(live | base_edge, 0);
            assert(base_clusters[v].adjacency.get_tail());
            assert(base_clusters[w].adjacency.get_tail());
            
            edge_cluster->add_ptr_to_highest_level(base_clusters[v].adjacency.get_tail());
            edge_cluster->add_ptr_to_highest_level(base_clusters[w].adjacency.get_tail());
            base_clusters[v].add_ptr_to_highest_level(edge_cluster->adjacency.get_tail());

        });

        edges_wflag = parlay::filter(edges_wflag, [] (auto ewf) {return ewf.valid;});

    }

    edges_wflag = parlay::tabulate(initial_edges.size(), [&] (T i) {
        edge_with_flag<std::tuple<T, T, D>> ret_edge;
        ret_edge.valid = true;
        ret_edge.E = initial_edges[i];
        return std::move(ret_edge);
    });

    while(edges_wflag.size())
    {
        // std::cout << "right side add edge sizes " << edges_wflag.size() << std::endl;
        parlay::parallel_for(0, edges_wflag.size(), [&] (T i) {
            auto& edge = edges_wflag[i].E;
            T v = std::get<0>(edge);
            T w = std::get<1>(edge);
            D weight = std::get<2>(edge);
            assert(v >= 0 && v < num_vertices);
            assert(w >= 0 && w < num_vertices);

            if(v < w)
                std::swap(v, w); // have a consistent view
            base_clusters[w].tiebreak = i + 1; // 
        });

        parlay::parallel_for(0, edges_wflag.size(), [&] (T i) {
            auto& edge = edges_wflag[i].E;
            T v = std::get<0>(edge);
            T w = std::get<1>(edge);
            D weight = std::get<2>(edge);
            if(v < w)
                std::swap(v, w); // have a consistent view
            if(base_clusters[w].tiebreak != i + 1)
                return;
            edges_wflag[i].valid = false;

            // std::string print_string = std::to_string(v) + " -- " + std::to_string(w) + " tiebreak succeeded\n";
            // std::cout << print_string;
            // std::cout << v << " -- " << w << " tiebreak succeeded" << std::endl;

            auto v_node_ptr = base_clusters[v].adjacency.get_tail();
            auto w_node_ptr = base_clusters[w].adjacency.get_tail();

            // base_clusters[v] should have myself
            for(auto& ptr : v_node_ptr->adjacents)
            {
                assert(ptr); // it should have already added me 
                if(get_other_side(v_node_ptr, ptr) == w_node_ptr)
                {
                    base_clusters[w].add_ptr_to_highest_level(ptr);
                    return;
                }
            }
            assert(false && "Shouldn't reach here, no joining edge found");
        });

        edges_wflag = parlay::filter(edges_wflag, [] (auto ewf) {return ewf.valid;});
    }

    // printTree(base_clusters);

    return;

}



/**
 * Make an exact copy of the current adjacency list, copying over everything (including the state), the only new thing will be the level.
 * When making a new round, pointers to edges will also get renewed, is that ideal? probably not but we don't have a choice do we 
 */
template<typename T, typename D>
void recreate_last_levels(parlay::sequence<node<T,D>*>& tree_nodes)
{
    
    
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        auto& cluster_ptr = node_ptr->cluster_ptr;
        auto& v = cluster_ptr->index;
        cluster_ptr->add_empty_level(node_ptr->state, node_ptr->contraction_level+1);

        for(short e = 0; e < node_ptr->adjacents.size(); e++)
        {
            auto edge_ptr = node_ptr->adjacents[e];
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
            auto& w = other_node_ptr->cluster_ptr->index;
            if(w < v)
                continue;
            edge_ptr->cluster_ptr->add_empty_level(edge_ptr->state, edge_ptr->contraction_level + 1); 
            node_ptr->next->adjacents[e] = edge_ptr->next;
            // edge_ptr->next->add_ptr(node_ptr->next);
            for(short k = 0; k < edge_ptr->adjacents.size(); k++)
            {
                if(edge_ptr->adjacents[k] == node_ptr)
                    edge_ptr->next->adjacents[k] = node_ptr->next;
            }
        }
        // for(auto& edge_ptr : node_ptr->adjacents)
        // {
        //     if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
        //         continue;
        //     auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
        //     auto& w = other_node_ptr->cluster_ptr->index;
        //     if(w < v)
        //         continue;
            
        //     // std::cout << "Binary added " << edge_ptr->index() << std::endl;

        //     edge_ptr->cluster_ptr->add_empty_level(edge_ptr->state, edge_ptr->contraction_level + 1); 
        //     cluster_ptr->add_ptr_to_highest_level(edge_ptr->next);
        //     edge_ptr->cluster_ptr->add_ptr_to_highest_level(node_ptr->next);
        // }
    });
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        auto& cluster_ptr = node_ptr->cluster_ptr;
        auto& v = cluster_ptr->index;

        for(short e = 0; e < node_ptr->adjacents.size(); e++)
        {
            auto edge_ptr = node_ptr->adjacents[e];
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
            auto& w = other_node_ptr->cluster_ptr->index;
            if(v < w)
                continue;
            node_ptr->next->adjacents[e] = edge_ptr->next;
            // edge_ptr->next->add_ptr(node_ptr->next);
            for(short k = 0; k < edge_ptr->adjacents.size(); k++)
            {
                if(edge_ptr->adjacents[k] == node_ptr)
                    edge_ptr->next->adjacents[k] = node_ptr->next;
            }
        }
    });

    parlay::parallel_for(0, tree_nodes.size(), [&] (T i){
        auto& node_ptr = tree_nodes[i];
        auto& cluster_ptr = node_ptr->cluster_ptr;
        for(short e = 0; e < node_ptr->adjacents.size(); e++)
        {
            auto edge_ptr = node_ptr->adjacents[e];
            if(edge_ptr == nullptr)
                continue;
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            if(edge_ptr->state & unary_cluster && !(edge_ptr->state & (binary_cluster | base_edge)))
            {
                edge_ptr->cluster_ptr->add_empty_level(unary_cluster, edge_ptr->contraction_level+1);
                edge_ptr->next->state = unary_cluster;
                node_ptr->next->adjacents[e] = edge_ptr->next;
            }
        }
        tree_nodes[i] = tree_nodes[i]->next;
    });
}


/**
 * The main workhorse, populates the set base_clusters with internal_clusters
 * As base_vertices aren't useful, it replaces them in place
 * 
 * The accumulation is simple for now -- each child points to the parent cluster
 * This can be used for connectivity tracking by starting at base_clusters[v] and base_clusters[w] and going to the parents
 * until they overlap
 * 
 * A more complex accumulator may be but currently the edges have no weights
 * This takes O(n) work and O(log n log n) span instead of O(log n log log n) span because the filter operatio is exact
 * not an approximate compaction
 * 
 * defretval is the default data stored in a cluster
*/

template <typename T, typename D, typename lambdafunc>
void create_RC_tree(parlay::sequence<cluster<T,D> > &base_clusters, T n, D defretval, lambdafunc assocfunc, bool randomized = false, bool print=false)
{


    auto tree_nodes = parlay::tabulate(base_clusters.size(), [&] (T i){
        return base_clusters[i].adjacency.get_tail();
    });

    // if(base_clusters.size() <= 100)
    //     printTree(base_clusters);

    recreate_last_levels(tree_nodes);

    auto count = 0;
    do
    {
        check_consistency(tree_nodes); // TODO remove
        // break;
        auto eligible = parlay::filter(tree_nodes, [] (auto node_ptr){
            return node_ptr->get_num_neighbours_live() <= 2;
        });

        
        set_MIS(eligible);
        auto candidates = parlay::filter(eligible, [] (auto node_ptr){
            return node_ptr->cluster_ptr->state & IS_MIS_SET;
        });

        check_mis(candidates); // TODO remove

        
        parlay::parallel_for(0, candidates.size(), [&] (T i){
            candidates[i]->cluster_ptr->first_contracted_node = candidates[i]->prev;
            contract(candidates[i]);
            accumulate(candidates[i], defretval, assocfunc);
        });

        
        tree_nodes = parlay::filter(tree_nodes, [] (auto node_ptr) {
            return node_ptr->state & live;
        });

        
        recreate_last_levels(tree_nodes);

        // if(base_clusters.size() <= 100)
        //     printTree(base_clusters);
    
        if (print)
        std::cout << "[" << count  << "]: " << bold << tree_nodes.size() << reset << std::endl;
        count++;
    }while(tree_nodes.size() > 0 && count < 100);

}

/**
 * Only frees the "floating" edge clusters
*/
template<typename T, typename D>
void deleteRCtree(parlay::sequence<cluster<T, D>> &base_clusters)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    parlay::parallel_for(0, base_clusters.size(), [&] (T i){
        auto& clstr = base_clusters[i];
        auto& node_ptr = clstr.first_contracted_node;
        if(node_ptr->state & binary_cluster)
        {
            node_ptr = node_ptr->prev;
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & base_edge)
                {
                    cluster_allocator::destroy(ptr->cluster_ptr);
                }
            }
        }
        else
        {
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & base_edge)
                {
                    cluster_allocator::destroy(ptr->cluster_ptr);
                }
            }
        }
    });
}

template<typename T, typename D>
void printTree(parlay::sequence<cluster<T, D>> &base_clusters, unsigned char level = -1)
{
    for(uint i = 0; i < base_clusters.size(); i++)
    {
        base_clusters[i].print(level);
    }   
}


#endif