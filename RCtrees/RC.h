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
#include "cluster.h"
#include "../examples/helper/graph_utils.h"
#include <parlay/alloc.h>



/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it

    Returns an array of parents such that the parents of index V would be parents[V]
*/
template <typename T>
parlay::sequence<T> generate_tree_graph(T num_elements)
{
    assert(num_elements > 0);

    parlay::sequence<T> dummy_initial_parents = parlay::tabulate(num_elements, [&] (T v) {return (T) 0;});

    std::mt19937 gen(std::random_device{}());
    

    parlay::parallel_for(0, num_elements, [&] (T v) {
        std::uniform_real_distribution<double> dis(0, 1);
        auto random_val = dis(gen);

        static const double anywhere_left_weight = 10;
        static const double immediate_left_weight = 30;
        static const double root_weight = 0.0;

        static const double anywhere_prob = (anywhere_left_weight/(anywhere_left_weight+immediate_left_weight+root_weight));
        static const double root_prob = anywhere_prob + (root_weight/(anywhere_left_weight+immediate_left_weight+root_weight));

        if (random_val <= anywhere_prob && v > 0)
        {
            std::uniform_int_distribution<> disint(0, v-1);

            dummy_initial_parents[v] = disint(gen);
        }
        else if (random_val < root_prob)
        {
            dummy_initial_parents[v] = v;
        }
        else
            dummy_initial_parents[v] = v == 0 ? 0 : v - 1;
    });  
    
    return dummy_initial_parents;
}

/*
    Converts the parents array into a symmetric graph
*/
template <typename graph, typename T>
graph convert_parents_to_graph(graph G, parlay::sequence<T> parents)
{
    parlay::sequence<T> vertices = parlay::tabulate(parents.size(), [&] (T v) {return v;});

    G = parlay::map(vertices, [&] (T v) {
        if(parents[v] == v) // root
        {
            parlay::sequence<T> empty;
            return empty;
        }
        parlay::sequence<T> temp = parlay::tabulate(1, [&] (T i) {return i;});
        temp[0] = parents[v];
        return temp;
    });

    G = graph_utils<T>::symmetrize(G);

    return G;
}

/**
 * Ensures that the degree is capped for the parents vector
 * i.e. not too many nodes have the same parent
 * Uses locking! Fortunately, we don't care about the performance for the graph generation too much
*/
template <typename T>
void degree_cap_parents(parlay::sequence<T> &parents, const T max_degree)
{
    parlay::sequence<std::atomic<T>> counts = parlay::tabulate(parents.size(), [] (size_t) {
       return std::atomic<T>(0); // Initialize each element with the value 0
    });


    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(counts[parents[v]].load() >= (max_degree))
        {
            parents[v] = v;
        }
        T parent_count = counts[parents[v]].fetch_add(1);
        if(parent_count >= (max_degree - 1))
        {
            parents[v] = v;
        }
    });

}


/**
 * Extracts a particular bit (counted from the right) from an element
*/
template <typename T>
inline bool extract_bit(T number, int offset_from_right)
{
    return (number >> offset_from_right) & 1;
}

/*
    Returns index of first different bit from the left
    Sets the value of bit in the boolean bit
    Sets the value of b
    sizeof(T) must be less than 16 bytes
*/
template <typename T>
inline char first_different_bit(const T a, const T b, bool* bit)
{
    T difference = a ^ b;
    char num_bits = sizeof(T) * 8;
    
    for(char i = num_bits-1; i >= 0; i--)
    {
        bool inspected_bit = extract_bit(difference, i);
        if(inspected_bit)
        {
            if(bit) *bit = extract_bit(b, i);
            return i;
        }
    }

    return -1;
}

/*
    returns a char with I_w and C_w(I_w) packed
    I_w is the index of the first different bit
    C_w is the value in neighbour w of this bit
    Also sets the different bit index in the char ptr different_bit_index

    Also, technically I waste the left-most bit in each char
*/

template <typename T>
static unsigned char get_single_colour_contribution(const T vcolour, const T wcolour, char* different_bit_index = NULL)
{
    bool wbit = false;
    char different_bit = first_different_bit(vcolour, wcolour, &wbit);
    char final_returned_character = (different_bit << 1) | wbit;
    if(different_bit_index) *different_bit_index = different_bit;
    return final_returned_character;
}





// /*
//     Given a set of cluster ptrs, colours them.
//     The clusters must have an initial valid colouring stored in their temp_colour folder
//     The clusters should be for "vertices" s.t. one hop from each vertex is a cluster representing an edge
//     and two hops away is a vertex representing a neighbouring node
// */
template<typename T, typename D>
void colour_clusters(parlay::sequence<cluster<T,D>*> clusters)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        unsigned long local_maximum = clusters[v]->get_default_colour();
        unsigned long local_minimum = clusters[v]->get_default_colour();

        for(uint i = 0; i < clusters[v]->size; i+=2)
        {

            if(clusters[v]->ptrs[i] == nullptr || (clusters[v]->types[i]&neighbour_type) == 0)
                continue;

            unsigned long compared_colour = clusters[v]->ptrs[i]->get_default_colour();

            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                local_minimum = compared_colour;   
        }

        if(local_maximum == clusters[v]->get_default_colour()) // This node is a local maximum, give it a unique colour
        {
            clusters[v]->colour = local_maximum_colour;
        }
        else if(local_minimum == clusters[v]->get_default_colour())
        {
            clusters[v]->colour = local_minimum_colour;
        }
        else
        {
            clusters[v]->colour = 2 + (get_single_colour_contribution(clusters[v]->get_default_colour(), local_maximum) / 2); // adding 2 and removing indicator bit
        }

    });

    return;

}

// /*
//     sets a boolean flag in the clusters indicating that they're part of MIS
//     also may change the boolean flag of some other clusters, only consider the clusters in this
//     These clusters must have a maximum degree of 2
// */
template<typename T, typename D>
void set_MIS(parlay::sequence<cluster<T,D>*> clusters, bool randomized = false)
{
    if(!randomized)
    {
        colour_clusters(clusters);

        parlay::parallel_for(0, clusters.size(), [&] (T v) {
            clusters[v]->state &= (~IS_MIS_SET);
            clusters[v]->set_neighbour_mis(false);
        });

        auto colours = parlay::tabulate(clusters.size(), [&] (T v) {
        return clusters[v]->colour;
        });

        auto vertices = parlay::tabulate(clusters.size(), [&] (T v) {
            return v;
        });

        auto result = vertices;

        parlay::sequence<unsigned long> offsets = counting_sort(vertices.begin(), vertices.end(), result.begin(), colours.begin(), 8 * sizeof(unsigned long));

        for(uint i = 0; i < offsets.size(); i++)
        {
            T start_index;
            if (i == 0)
                start_index = 0;
            else
                start_index = offsets[i-1];
            T end_index = offsets[i];

            parlay::parallel_for(start_index, end_index, [&] (T i) {
                T v = result[i];
                if(clusters[v]->get_neighbour_MIS() == true)
                {
                    clusters[v]->state &= (~IS_MIS_SET);
                    return;
                }
                clusters[v]->state |= IS_MIS_SET;
            });
        }
    }
    else
    {
        // set initial colour to zero
        parlay::parallel_for(0, clusters.size(), [&] (T v) {
            clusters[v]->set_neighbour_colour(0);
        });

        parlay::random_generator gen;
        std::uniform_int_distribution<T> dis(1, clusters.size());
        parlay::parallel_for(0, clusters.size(), [&] (T v) {
            auto r = gen[v];
            clusters[v]->colour = dis(r);
        });

        parlay::parallel_for(0, clusters.size(), [&] (T v) {
            bool flag = clusters[v]->is_max_neighbour_colour();
            if(flag)
                clusters[v]->state |= IS_MIS_SET;
            else
                clusters[v]->state &= (~IS_MIS_SET);
        });

    }

}


/**
 * Given an symmetric graph, creates a set of clusters
 * In total, it creates n + m clusters in the array base_clusters
 * The first n are base_vertex clusters
 * And the last m are base_edge clusters
 * These clusters are linked to point to each other appropriately
 * it returns an array of n clusters corresponding to the original n vertices.
 * The edge clusters are present as edges.
*/
template <typename T, typename D>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T, D> > &base_clusters, const T max_size, D defretval = 0.00f)
{

    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    T n = G.size();
    T total_num_clusters = n;

    base_clusters = parlay::tabulate(total_num_clusters, [&] (T v) {
        cluster<T,D> base_cluster;
        return base_cluster;
    });


    // Add all "outgoung" edges from one side
    parlay::parallel_for(0, n, [&] (T v) {
        auto& _cluster = base_clusters[v];
        _cluster.index = v;
        _cluster.state = base_vertex | live;
        
        for(const auto& w : G[v])
        {
            _cluster.add_initial_adjacency(w);
            if(w < v)
                continue;
            auto edge_cluster = cluster_allocator::alloc();
            edge_cluster->index = -1;
            edge_cluster->state = base_edge | live;
            edge_cluster->data = defretval;
            edge_cluster->add_initial_neighbours(&base_clusters[w], &base_clusters[v]);
            // std::cout << "[" << v << "] Edge added between " << v << "<->" << w << std::endl;
            _cluster.add_neighbour(&base_clusters[w], edge_cluster);
        }
    });

    // // Add incoming edges
    parlay::parallel_for(0, n, [&] (T v) {
        auto cluster_ptr = &base_clusters[v];
        cluster_ptr->index = v;

        for(const auto& w : G[v])
        {
            if(v < w) // reversed
                continue;
            auto other_node_ptr = &base_clusters[w];
            // find the edge that corresponds that joins to the other side
            for(short i = 0; i < cluster_ptr->size; i+=2)
            {
                if(other_node_ptr->ptrs[i] == cluster_ptr)
                {
                    auto& connecting_edge_ptr = other_node_ptr->ptrs[i+1];
                    cluster_ptr->add_neighbour(other_node_ptr, connecting_edge_ptr);
                    break;
                }
            }
        }

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

template <typename T, typename D>
void create_RC_tree(parlay::sequence<cluster<T,D> > &base_clusters, T n, bool randomized = false, D defretval = 0.00f)
{
    // std::cout << "create RC tree called" << std::endl;
    
    parlay::sequence<cluster<T,D>*> all_cluster_ptrs = parlay::tabulate(n, [&] (T v) {
        return &base_clusters[v];
    });

    parlay::sequence<cluster<T,D>*> forest, candidates;

    // Initially the forest of live nodes is all live nodes
    forest = parlay::filter(all_cluster_ptrs, [&] (cluster<T,D>* C) {
        return ((C->state & live));
    });

    bool first_time = true;

    // printTree(base_clusters);

    unsigned char contraction_round = 0;

    do
    {
        // Shrink the forst as live nodes decrease
        if(first_time == true)
        {
            
            first_time = false;
        }
        else
        {
            forest = parlay::filter(forest, [&] (cluster<T,D>* C) {
                return ((C->state & live));
            });
        }
        
        // std::cout << "forest.size(): " << forest.size() << std::endl;

        // Eligible nodes are those with 0, 1 or 2 neighbours
        auto eligible = parlay::filter(forest, [&] (cluster<T,D>* C) {
            return (C->get_neighbour_count() <= 2);
        });

        // std::cout << "eligible.size(): " << eligible.size() << std::endl;

        // Set the flag is_MIS amongst them
        set_MIS(eligible, randomized = randomized);

        // Filter out an MIS of eligible nodes
        candidates = parlay::filter(eligible, [&] (cluster<T,D>* C) {
            return (C->state & IS_MIS_SET);
        });

        // std::cout << "candidates.size(): " << candidates.size() << std::endl;

        // printTree(base_clusters);

        // do rake and compress
        parlay::parallel_for(0, candidates.size(), [&] (T v) {
            cluster<T,D>* cluster_ptr = candidates[v];
            cluster_ptr->data = defretval;
            // rake
            if(cluster_ptr->get_neighbour_count() == 0)
            {
                cluster_ptr->state&=(~live);
                cluster_ptr->state|=(nullary_cluster);
                cluster_ptr->state|=internal;
                cluster_ptr->contraction_time = contraction_round;
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
                cluster_ptr->contraction_time = contraction_round;
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

                left_edge_ptr->state&=(~live);
                right_edge_ptr->state&=(~live);

                cluster_ptr->state&=(~live);

                cluster_ptr->state|=binary_cluster;
                cluster_ptr->state|=internal;
                cluster_ptr->contraction_time = contraction_round;
            }
        });
        // printTree(base_clusters);
        // std::cout << "Forest size: " << forest.size() << std::endl;
        contraction_round++;
    }while(forest.size());

    // not necessary anynmore!
    // if(do_height == true)
    // set_heights(all_cluster_ptrs);

}





template <typename T, typename D, typename assocfunc>
D PathQuery( cluster<T, D>* v,  cluster<T, D>* w, const D& defretval, assocfunc func)
{
    
    if(v == w)
        return defretval;

    bool use_v = false;

    cluster<T,D>* prev_boundary_v_l = nullptr;
    cluster<T,D>* prev_boundary_v_r = nullptr;
    D prev_val_till_v_l = defretval;
    D prev_val_till_v_r = defretval;
    
    cluster<T,D>* prev_boundary_w_l = nullptr;
    cluster<T,D>* prev_boundary_w_r = nullptr;
    D prev_val_till_w_l = defretval;
    D prev_val_till_w_r = defretval;


    while(true)
    {
        // Ascend using the lower-height cluster. If either has to ascend above root, they are not connected.
        if(v == nullptr && w == nullptr && prev_boundary_v_l == nullptr && prev_boundary_v_r == nullptr)
            return defretval;   
        else if(v == w)
            break;
        else if (v == nullptr && w != nullptr)
            use_v = false;
        else if (w == nullptr && v != nullptr)
            use_v = true;
        else if(v->get_height() < w->get_height())
            use_v = true;
        else
            use_v = false;
        
        if(use_v)
        {
            cluster<T,D>* &v_cluster = v;

            // We are at the start
            if (prev_val_till_v_l == defretval && prev_val_till_v_r == defretval) 
            {
                v_cluster->find_boundary_vertices(prev_boundary_v_l, prev_val_till_v_l, prev_boundary_v_r, prev_val_till_v_r, defretval);
                
            } 
            else 
            {
                D lval;
                D rval;
                cluster<T,D>* l;
                cluster<T,D>* r;

                v_cluster->find_boundary_vertices(l, lval, r, rval, defretval);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_v_l == v && prev_boundary_v_r == v) 
                {
                    prev_boundary_v_l = l;
                    prev_boundary_v_r = r;
                    
                    if(l != r) // ascending into unary
                    {
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    }
                    else
                        for (uint i = 0; i < v_cluster->size; i++) // any children which are edges
                        {
                            if (v_cluster->types[i] & child_type && v_cluster->ptrs[i]->state & (binary_cluster | base_edge)) 
                            {
                                prev_val_till_v_l = func(prev_val_till_v_l, v_cluster->ptrs[i]->data);
                                prev_val_till_v_r = func(prev_val_till_v_r, v_cluster->ptrs[i]->data);
                            }
                        }
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if(prev_boundary_v_r == v)
                    {
                        prev_val_till_v_r = prev_val_till_v_l;
                    }
                    else
                        prev_val_till_v_l = prev_val_till_v_r;

                    prev_boundary_v_l = prev_boundary_v_r = l;
                }
                // Binary to binary
                else 
                {
                    if (prev_boundary_v_l == l) 
                    {
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    } 
                    else if (prev_boundary_v_r == l) 
                    {
                        std::swap(prev_boundary_v_r, prev_boundary_v_l);
                        std::swap(prev_val_till_v_l, prev_val_till_v_r);
                        
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    
                    } 
                    else if (prev_boundary_v_l == r) 
                    {
                        std::swap(r, l);
                        std::swap(rval, lval);    
                        
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    } 
                    else // prev_boundary_v_r == r 
                    {     
                        prev_boundary_v_l = l;
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                    }
                }
            }
            // std::cout << bright_magenta << v->index << "[" << (short) v->get_height() << "]" << " ";
            // std::cout << prev_boundary_v_l->index << "(" << prev_val_till_v_l << ") " << prev_boundary_v_r->index << "(" << prev_val_till_v_r << ") ";
            // auto& for_printing = v;
            // if(for_printing->state & unary_cluster)
            //     std::cout << " unary ";
            // else if (for_printing->state & binary_cluster)
            //     std::cout << " binary ";
            // else
            //     std::cout << " nullary ";
            // std::cout << reset << std::endl;
            v = v->get_parent();
        }
        else
        {
            cluster<T,D>* &w_cluster = w;

            // We are at the start
            if (prev_val_till_w_l == defretval && prev_val_till_w_r == defretval) 
            {
                w_cluster->find_boundary_vertices(prev_boundary_w_l, prev_val_till_w_l, prev_boundary_w_r, prev_val_till_w_r, defretval);
            } 
            else 
            {
                D lval, rval;
                cluster<T,D>* l;
                cluster<T,D>* r;
                w_cluster->find_boundary_vertices(l, lval, r, rval, defretval);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_w_l == w && prev_boundary_w_r == w) 
                {
                    
                    prev_boundary_w_l = l;
                    prev_boundary_w_r = r;

                    if(l != r)
                    {
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    }
                    
                    else
                        for (uint i = 0; i < w_cluster->size; i++) // any children which are edges
                        {
                            if (w_cluster->types[i] & child_type && w_cluster->ptrs[i]->state & (binary_cluster | base_edge)) 
                            {
                                prev_val_till_w_l = func(prev_val_till_w_l, w_cluster->ptrs[i]->data);
                                prev_val_till_w_r = func(prev_val_till_w_r, w_cluster->ptrs[i]->data);
                            }
                        }
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if(prev_boundary_w_r == w)
                    {
                        prev_val_till_w_r = prev_val_till_w_l;
                    }
                    else
                        prev_val_till_w_l = prev_val_till_w_r;

                    prev_boundary_w_l = prev_boundary_w_r = l;
                }
                // Binary
                else 
                {
                    if (prev_boundary_w_l == l) 
                    {
                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    } 
                    else if (prev_boundary_w_r == l) 
                    {
                        
                        std::swap(prev_boundary_w_l, prev_boundary_w_r);
                        std::swap(prev_val_till_w_l, prev_val_till_w_r);

                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);

                    } 
                    else if (prev_boundary_w_l == r) 
                    {
                        std::swap(l, r);
                        std::swap(lval, rval);

                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    } 
                    else // prev_boundary_w_r == r 
                    { 
                        prev_boundary_w_l = l;
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                    }
                }
            }

            // std::cout << bright_green << w->index << "[" << (short) w->get_height() << "]" << " ";
            // std::cout << prev_boundary_w_l->index << "(" << prev_val_till_w_l << ") " << prev_boundary_w_r->index << "(" << prev_val_till_w_r << ") ";
            // auto& for_printing = w;
            // if(for_printing->state & unary_cluster)
            //     std::cout << " unary ";
            // else if (for_printing->state & binary_cluster)
            //     std::cout << " binary ";
            // else
            //     std::cout << " nullary ";
            // std::cout << reset << std::endl;
            w = w->get_parent();
        }
        
    }
    std::cout << std::endl;

    if(prev_boundary_v_l == nullptr && prev_boundary_w_l == nullptr)
    {
        return defretval;
    }
    if(prev_boundary_v_l == nullptr)
    {
        if(w == prev_boundary_w_r)
            return prev_val_till_w_r;
        else
            return prev_val_till_w_l;   
    }
    if(prev_boundary_w_l == nullptr)
    {
        if(v == prev_boundary_v_l)
            return prev_val_till_v_l;
        else
            return prev_val_till_v_r;
    }
    
    D v_contrib, w_contrib;
    
    if(v == prev_boundary_v_l)
        v_contrib = prev_val_till_v_l;
    else
        v_contrib = prev_val_till_v_r;

    if(w == prev_boundary_w_l)
        w_contrib = prev_val_till_w_l;
    else
        w_contrib = prev_val_till_w_r;

    return func(v_contrib, w_contrib);
   
}

/**
 * Given a pair [v,w], returns the  (if it exists)
 * Exploits the fact that an edge [v,w] must be the child of v or w in an RC tree
*/
template<typename T, typename D>
cluster<T, D>* getEdge(T v, T w, parlay::sequence<cluster<T, D>>& clusters)
{

    if(v == w)
        return nullptr;
    

    // check v
    auto& cluster_v = clusters[v];
    for(uint i = 0; i < cluster_v.size; i+=2)
    {
        if(cluster_v.ptrs[i] != nullptr && cluster_v.ptrs[i]->index == w && cluster_v.ptrs[i+1]->state & base_edge)
            return cluster_v.ptrs[i+1];
    }

    // check w
    auto& cluster_w = clusters[w];
    for(uint i = 0; i < cluster_w.size; i+=2)
    {
        if(cluster_w.ptrs[i] != nullptr && cluster_w.ptrs[i]->index == v && cluster_w.ptrs[i+1]->state & base_edge)
            return cluster_w.ptrs[i+1];
    }

    std::cout << red << v << " - " << w <<": edge not found" << reset << std::endl;

    return nullptr;
}

template<typename T, typename D, typename assocfunc>
void batchModifyEdgeWeights(const parlay::sequence<std::tuple<T, T, D>>& edges, assocfunc func, parlay::sequence<cluster<T, D>>& clusters, const D& defretval = 0.00f)
{

    // increment counter and return ptrs
    parlay::sequence<cluster<T,D>*> edge_ptrs = parlay::tabulate(edges.size(), [&] (const T i) {

        const std::tuple<T, T, D> edge_tuple = edges[i];
        const auto& weight = std::get<2>(edge_tuple);
        const auto& v = std::get<1>(edge_tuple);
        const auto& w = std::get<0>(edge_tuple);

        cluster<T,D>* cluster_ptr = getEdge(v, w, clusters);

        if(cluster_ptr == nullptr)
            return cluster_ptr;
        cluster_ptr->data = weight;

        auto retptr = cluster_ptr;


        while(cluster_ptr != nullptr)
        {
            auto old_value = cluster_ptr->counter.fetch_add(1);
            if(old_value > 0) 
                break;
            cluster_ptr = cluster_ptr->get_parent();
        }

        return retptr;

    });

    // std::cout << bright_magenta;
    // for(const auto& e : edges)
    //     std::cout << std::get<0>(e) << "<->" << std::get<1>(e) << ":" << std::get<2>(e) << " ";
    // std::cout << reset << std::endl << std::endl;



    parlay::parallel_for(0, edge_ptrs.size(), [&] (T i) {
        if(edge_ptrs[i] == nullptr)
            return;

        auto cluster_ptr = edge_ptrs[i];

        // decrement counter
        while(cluster_ptr != nullptr)
        {
            auto ret_val = cluster_ptr->counter.fetch_add(-1);
            if(ret_val > 1)
                break;
            // prepare children
            bool first_child = true;
            if(!(cluster_ptr->state & base_edge))
            {
                D final_val = defretval; 
                for(uint i = 0; i < cluster_ptr->size; i++)
                {
                    if(cluster_ptr->ptrs[i] != nullptr && cluster_ptr->types[i] & child_type && (cluster_ptr->ptrs[i]->state & (binary_cluster | base_edge)))
                    {
                        if(first_child == true)
                        {
                            first_child = false;
                            final_val = cluster_ptr->ptrs[i]->data;
                        }
                        else
                        {
                            final_val = func(final_val, cluster_ptr->ptrs[i]->data);
                        }
                    }
                }
                cluster_ptr->data = final_val;
            }
            cluster_ptr = cluster_ptr->get_parent();
        }

    });
    
    
    // for(const auto& ep : edge_ptrs)
    //     if(ep)
    //         std::cout << ep->data << " ";
    // std::cout << std::endl;
    

    return;
}



/**
 * Only frees the "floating" edge clusters
*/
template<typename T, typename D>
void deleteRCtree(parlay::sequence<cluster<T, D>> &base_clusters)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    parlay::parallel_for(0, base_clusters.size(), [&] (T v) {
        auto cluster_ptr = &base_clusters[v];
        
        for(short i = 0; i < cluster_ptr->size; i+=1)
        {
            cluster<T,D>*& potential_ptr = cluster_ptr->ptrs[i];
            short potential_type = cluster_ptr->types[i];
            
            if((potential_ptr != nullptr) && (potential_ptr->state & base_edge) && (potential_type & (child_type)))
            {
                cluster_allocator::destroy(potential_ptr);
            }
        }

    });
}

template<typename T, typename D>
void printTree(parlay::sequence<cluster<T, D>> &base_clusters)
{
    for(uint i = 0; i < base_clusters.size(); i++)
    {
        base_clusters[i].print();
    }   
}

#endif