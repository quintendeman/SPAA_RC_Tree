#include <atomic>
#include <random>
#include <set>
#include <iostream>
#include <mutex>
#include "cluster.h"
#include "../examples/helper/graph_utils.h"


#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"




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
        std::uniform_real_distribution<> dis(0, 1);
        static const float nullary_weight = 5.0;
        static const float anywhere_weight = 500;
        static const float adjacent_weight = 50;
        static const float nullary = (nullary_weight)/(nullary_weight+anywhere_weight+adjacent_weight);
        static const float anywhere = (anywhere_weight)/(nullary_weight+anywhere_weight+adjacent_weight);
        auto random_val = dis(gen);
        if (random_val <= anywhere && v > 0)
        {
            std::uniform_int_distribution<T> disint(0, v-1);

            dummy_initial_parents[v] = disint(gen);
        }
        else if (random_val < (anywhere + nullary))
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
static inline bool extract_bit(T number, int offset_from_right)
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
template<typename T>
void colour_clusters(parlay::sequence<T>& cluster_indices, parlay::sequence<cluster<T>>& all_clusters)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, cluster_indices.size(), [&] (T v) {
        cluster<T>& my_cluster = all_clusters[cluster_indices[v]];

        unsigned long local_maximum = my_cluster.get_default_colour();
        unsigned long local_minimum = my_cluster.get_default_colour();

        for(uint i = 0; i < my_cluster.size; i+=2)
        {

            if(my_cluster.indices[i] == -1 || (my_cluster.types[i]&neighbour_type) == 0)
                continue;

            unsigned long compared_colour = all_clusters[my_cluster.indices[i]].get_default_colour();

            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                local_minimum = compared_colour;   
        }

        if(local_maximum == my_cluster.get_default_colour()) // This node is a local maximum, give it a unique colour
        {
            my_cluster.colour = local_maximum_colour;
        }
        else if(local_minimum == my_cluster.get_default_colour())
        {
            my_cluster.colour = local_minimum_colour;
        }
        else
        {
            my_cluster.colour = 2 + (get_single_colour_contribution(my_cluster.get_default_colour(), local_maximum) / 2); // adding 2 and removing indicator bit
        }

    });

    return;

}

// /*
//     sets a boolean flag in the clusters indicating that they're part of MIS
//     also may change the boolean flag of some other clusters, only consider the clusters in this
//     These clusters must have a maximum degree of 2
// */
template<typename T>
void set_MIS(parlay::sequence<T>& cluster_indices,  parlay::sequence<cluster<T>>& all_clusters, bool randomized = false)
{
    if(!randomized)
    {
        colour_clusters(cluster_indices, all_clusters);


        auto colours = parlay::tabulate(cluster_indices.size(), [&] (T v) {
            all_clusters[cluster_indices[v]].set_MIS(false);
            all_clusters[cluster_indices[v]].state |= is_candidate;
            return all_clusters[cluster_indices[v]].colour;
        });

        auto vertices = parlay::tabulate(cluster_indices.size(), [&] (T v) {
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
                if(all_clusters[cluster_indices[v]].get_neighbour_MIS(all_clusters) == true)
                {
                    all_clusters[cluster_indices[v]].set_MIS(false);
                    return;
                }
                all_clusters[cluster_indices[v]].set_MIS(true);
            });
        }

        parlay::parallel_for(0, cluster_indices.size(), [&] (T v) {
            all_clusters[cluster_indices[v]].state &= (~is_candidate);
        });
    }
    else
    {
        // set as is_candidate
        parlay::parallel_for(0, cluster_indices.size(), [&] (T v) {
            all_clusters[cluster_indices[v]].state |= is_candidate;
        });

        parlay::random_generator gen;
        std::uniform_int_distribution<T> dis(1, cluster_indices.size());
        parlay::parallel_for(0, cluster_indices.size(), [&] (T v) {
            auto r = gen[v];
            all_clusters[cluster_indices[v]].colour = dis(r);
        });

        parlay::parallel_for(0, cluster_indices.size(), [&] (T v) {
            all_clusters[cluster_indices[v]].set_MIS(all_clusters[cluster_indices[v]].is_max_neighbour_colour(all_clusters));
        });

        parlay::parallel_for(0, cluster_indices.size(), [&] (T v) {
            all_clusters[cluster_indices[v]].state &= (~is_candidate);
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
template <typename T>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T> > &base_clusters, const T max_size)
{

    T n = G.size();
    T total_num_clusters = n + n-1;

    base_clusters = parlay::tabulate(total_num_clusters, [&] (T v) {
        cluster<T> base_cluster;
        return base_cluster;
    });

    parlay::parallel_for(0, n, [&] (T v) {
        cluster<T>& cluster = base_clusters[v];
        cluster.index = v;
        cluster.state = base_vertex | live;

        for(short i = 0; i < G[v].size(); i++) //TODO: Assumes ordering i.e. root is towards zero
        {
            if(G[v][i] < v)
            {
                base_clusters[v+n].add_initial_neighbours(G[v][i], v);
                base_clusters[v+n].state = base_edge;
                base_clusters[v+n].index = v+n;
                base_clusters[v].add_neighbour(G[v][i], v+n);
            }
            else
            {
                base_clusters[v].add_neighbour(G[v][i], G[v][i]+n);
            }
        }

    });

}

template <typename T>
void set_heights(parlay::sequence<T> &all_cluster_indices, parlay::sequence<cluster<T>>& all_clusters)
{

    parlay::parallel_for(0, all_cluster_indices.size(), [&] (T v) {
        auto cluster = &all_clusters[all_cluster_indices[v]];
        T my_height = 0;
        T parent_index = cluster->get_parent();
        cluster->height = 0;

        while(parent_index != -1)
        {
            bool swapped = false;
            while(swapped == false)
            {                
                T height = cluster->height.load();
                if(height >= my_height)
                {
                    return;
                }
                swapped = cluster->height.compare_exchange_strong(height, my_height);
            }            

            my_height++;
            parent_index = cluster->get_parent();
            cluster = &all_clusters[parent_index];
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
*/

template <typename T>
void create_RC_tree(parlay::sequence<cluster<T> > &base_clusters, T n, bool do_height = true, bool randomized = false, bool print_forest = false)
{
    // std::cout << "create RC tree called" << std::endl;
    
    parlay::sequence<T> all_cluster_indices = parlay::tabulate(n, [&] (T v) {
        return v;
    });

    parlay::sequence<T> forest, candidates;

    // Initially the forest of live nodes is all live nodes
    forest = parlay::filter(all_cluster_indices, [&] (T C) {
        
        return base_clusters[C].state & live;
    });

    auto last_forest_size = 0;
    auto count = 0;
    bool first_time = true;

    do
    {
       

        // printTree(base_clusters);
        // std::cout << "\n\n\n\n\n\n" << std::endl;
        // Shrink the forst as live nodes decrease
        if(first_time == true)
        {
            
            first_time = false;
        }
        else
        {
            forest = parlay::filter(forest, [&] (T C) {
                return base_clusters[C].state & live;
            });
            if(forest.size() == last_forest_size)
            {
                count++;
            }
            assert(count < 100 && "INFINITE LOOP, FOREST SIZE CONSTANT");
            last_forest_size = forest.size();
        }
        
        // std::cout << "forest.size(): " << forest.size() << std::endl;

        // Eligible nodes are those with 0, 1 or 2 neighbours
        auto eligible = parlay::filter(forest, [&] (T C) {
            return base_clusters[C].get_neighbour_count() <= 2;
        });

        // std::cout << "eligible.size(): " << eligible.size() << std::endl;

        // Set the flag is_MIS amongst them
        set_MIS(eligible, base_clusters, randomized = randomized);

        // Filter out an MIS of eligible nodes
        candidates = parlay::filter(eligible, [&] (T C) {
            return base_clusters[C].state & is_MIS;
        });

        // std::cout << "candidates.size(): " << candidates.size() << std::endl;

        
        

        // do rake and compress
        parlay::parallel_for(0, candidates.size(), [&] (T v) {
            T cluster_index = candidates[v];
            // base_clusters[cluster_index];
            // rake
            if(base_clusters[cluster_index].get_neighbour_count() == 0)
            {
                base_clusters[cluster_index].state&=(~live);
                base_clusters[cluster_index].state|=(nullary_cluster);
                base_clusters[cluster_index].state|=internal;
            }
            if(base_clusters[cluster_index].get_neighbour_count() == 1)
            {
                T edge_idx = -1;
                T other_side = -1;
                short neighbour_index = -1;
                
                for(short i = 0; i < base_clusters[cluster_index].size; i+=2) // TODO: replace with a get_first_neighbour functions?
                {
                    if(base_clusters[cluster_index].types[i] & neighbour_type)
                    {
                        other_side = base_clusters[cluster_index].indices[i];
                        edge_idx = base_clusters[cluster_index].indices[i+1];
                        neighbour_index = i;
                    }
                }
                // assert(edge_idx != -1 && other_side != -1);

                // now add these two as children and remove as neighbours if neighbours
                auto child1 = base_clusters[other_side].change_to_child(edge_idx);
                auto child2 = base_clusters[other_side].change_to_child(cluster_index);
                
                // assert(child1 != -1 && child2 != -1);

                // // make other_side be the parent for both
                base_clusters[edge_idx].set_parent(cluster_index);
                base_clusters[cluster_index].set_parent(other_side);

                //
                base_clusters[cluster_index].state&=(~live);
                base_clusters[cluster_index].state|=unary_cluster;
                base_clusters[cluster_index].state|=internal;
            }
            if (base_clusters[cluster_index].get_neighbour_count() == 2)
            {
                // find left and right vertices/nodes
                T left_edge_idx = -1;
                T right_edge_idx = -1;
                T left_node_idx = -1;
                T right_node_idx = -1;

                base_clusters[cluster_index].get_two_neighbours_edges(left_node_idx, left_edge_idx, right_node_idx, right_edge_idx);

                // assert(left_edge_idx != -1);
                // assert(right_edge_idx != -1);
                
                base_clusters[left_node_idx].overwrite_neighbour(cluster_index, right_node_idx, cluster_index);

                base_clusters[right_node_idx].overwrite_neighbour(cluster_index, left_node_idx, cluster_index);
                
                base_clusters[left_edge_idx].set_parent(cluster_index);
                base_clusters[right_edge_idx].set_parent(cluster_index);

                auto child1 = base_clusters[cluster_index].change_to_child(left_edge_idx);
                auto child2 = base_clusters[cluster_index].change_to_child(right_edge_idx);

                

                // assert(child1 != -1 && child2 != -1);

                // base_clusters[left_edge_idx].state&=(~live);
                // base_clusters[right_edge_idx].state&=(~live);

                base_clusters[cluster_index].state&=(~live);

                base_clusters[cluster_index].state|=binary_cluster;
                base_clusters[cluster_index].state|=internal;
                
            }
        });

        // std::cout << "forest: " << forest.size() << std::endl;

    }while(forest.size());

    if(do_height == true)
        set_heights(all_cluster_indices, base_clusters);

}

template <typename T, typename assocfunc>
T PathQuery(T v, T w, T defretval, parlay::sequence<cluster<T> > &base_clusters, assocfunc func)
{
    
    bool use_v = false;

    T prev_boundary_v_l = base_clusters[v].get_parent();
    T prev_boundary_v_r = prev_boundary_v_l;
    T prev_val_till_v_l = defretval;
    T prev_val_till_v_r = defretval;
    T prev_v = v;


    T prev_boundary_w_l = base_clusters[w].get_parent();
    T prev_boundary_w_r = prev_boundary_w_l;
    T prev_val_till_w_l = defretval;
    T prev_val_till_w_r = defretval;



    while(true)
    {
        // Ascend using the lower-height cluster. If either has to ascend above root, they are not connected.
        if(v == -1 || w == -1)
            return defretval;   
        else if(v == w)
            break;
        else if(base_clusters[v].get_height() < base_clusters[w].get_height())
            use_v = true;
        else
            use_v = false;
        


        if(use_v)
        {
            cluster<T>& v_cluster = base_clusters[v]; 

            // We are at the start
            if(prev_val_till_v_l == defretval && prev_val_till_v_r == defretval)
            {   
                v_cluster.find_boundary_vertices(prev_boundary_v_l, prev_val_till_v_l, prev_boundary_v_r, prev_val_till_v_r, defretval, base_clusters);
            }
            else
            {
                T l, lval, r, rval;
                v_cluster.find_boundary_vertices(l, lval, r, rval, defretval, base_clusters);
                    
                // covers both unary to unary
                // and unary to binary case
                if(prev_boundary_v_l == v && prev_boundary_v_r == v) 
                {
                    prev_boundary_v_l = l;
                    prev_boundary_v_r = r;
                    prev_val_till_v_l = func(prev_val_till_v_l, lval);
                    prev_val_till_v_r = func(prev_val_till_v_r, rval);

                    for(uint i = 0; i < v_cluster.size; i++)
                    {
                        if(v_cluster.types[i] & child_type)
                        {   
                            prev_val_till_v_l = func(prev_val_till_v_l, base_clusters[v_cluster.indices[i]].data);
                            prev_val_till_v_r = func(prev_val_till_v_r, base_clusters[v_cluster.indices[i]].data);
                        }
                    }
                }
                // covers binary ascending to:
                // unary/nullary
                else if(l == r)
                {
                    if(prev_boundary_v_l == v)
                    {
                        prev_boundary_v_l = l;
                        prev_val_till_v_l = func(prev_boundary_v_l, lval);
                    }
                    else
                    {
                        prev_boundary_v_l = r;
                        prev_val_till_v_r = func(prev_boundary_v_r, rval);
                    }
                    for(uint i = 0; i < v_cluster.size; i++)
                    {
                        if(v_cluster.types[i] & child_type && v_cluster.indices[i] != l)
                        {   
                            prev_val_till_v_l = func(prev_val_till_v_l, base_clusters[v_cluster.indices[i]].data);
                            prev_val_till_v_r = func(prev_val_till_v_r, base_clusters[v_cluster.indices[i]].data);
                        }
                    }
                    
                }
                // binary
                else
                {
                    if(prev_boundary_v_l == l)
                    {
                        prev_boundary_v_l = l; // not necessary but helps me thinks
                        prev_val_till_v_l = prev_boundary_v_l;
                        prev_boundary_v_r = r;
                        prev_boundary_v_r = func(prev_boundary_v_r, rval);
                        prev_boundary_v_r = func(prev_boundary_v_r, lval);
                    }
                    else if(prev_boundary_v_r == l)
                    {
                        prev_boundary_v_r = l; // not necessary but helps me thinks
                        prev_val_till_v_r = prev_boundary_v_r;
                        prev_boundary_v_l = r;
                        prev_boundary_v_l = func(prev_boundary_v_l, rval);
                        prev_boundary_v_l = func(prev_boundary_v_l, lval);
                        
                    }
                    else if(prev_boundary_v_l == r)
                    {
                        prev_boundary_v_l = r; // not necessary but helps me thinks
                        prev_val_till_v_l = prev_boundary_v_l;
                        prev_boundary_v_r = l;
                        prev_boundary_v_r = func(prev_boundary_v_r, rval);
                        prev_boundary_v_r = func(prev_boundary_v_r, lval);
                    }
                    else /*prev_boundary_v_r == r*/
                    {
                        prev_boundary_v_r = r; // not necessary but helps me thinks
                        prev_val_till_v_r = prev_boundary_v_r;
                        prev_boundary_v_l = l;
                        prev_boundary_v_l = func(prev_boundary_v_l, rval);
                        prev_boundary_v_l = func(prev_boundary_v_l, lval);
                    }
                }
            }
            
            v = base_clusters[v].get_parent();
        }
        else
        {
            cluster<T>& w_cluster = base_clusters[w];

            // We are at the start
            if (prev_val_till_w_l == defretval && prev_val_till_w_r == defretval) 
            {
                w_cluster.find_boundary_vertices(prev_boundary_w_l, prev_val_till_w_l, prev_boundary_w_r, prev_val_till_w_r, defretval, base_clusters);
            } 
            else 
            {
                T l, lval, r, rval;
                w_cluster.find_boundary_vertices(l, lval, r, rval, defretval, base_clusters);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_w_l == w && prev_boundary_w_r == w) 
                {
                    prev_boundary_w_l = l;
                    prev_boundary_w_r = r;
                    prev_val_till_w_l = func(prev_val_till_w_l, lval);
                    prev_val_till_w_r = func(prev_val_till_w_r, rval);

                    for (uint i = 0; i < w_cluster.size; i++) 
                    {
                        if (w_cluster.types[i] & child_type) 
                        {
                            prev_val_till_w_l = func(prev_val_till_w_l, base_clusters[w_cluster.indices[i]].data);
                            prev_val_till_w_r = func(prev_val_till_w_r, base_clusters[w_cluster.indices[i]].data);
                        }
                    }
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if (prev_boundary_w_l == w) 
                    {
                        prev_boundary_w_l = l;
                        prev_val_till_w_l = func(prev_boundary_w_l, lval);
                    } 
                    else 
                    {
                        prev_boundary_w_l = r;
                        prev_val_till_w_r = func(prev_boundary_w_r, rval);
                    }
                    for (uint i = 0; i < w_cluster.size; i++) 
                    {
                        if (w_cluster.types[i] & child_type && w_cluster.indices[i] != l) 
                        {
                            prev_val_till_w_l = func(prev_val_till_w_l, base_clusters[w_cluster.indices[i]].data);
                            prev_val_till_w_r = func(prev_val_till_w_r, base_clusters[w_cluster.indices[i]].data);
                        }
                    }
                }
                // Binary
                else 
                {
                    if (prev_boundary_w_l == l) 
                    {
                        prev_boundary_w_l = l; // Not necessary but helps me think
                        prev_val_till_w_l = prev_val_till_w_l;
                        
                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, lval);
                    } 
                    else if (prev_boundary_w_r == l) 
                    {
                        prev_boundary_w_r = l; // Not necessary but helps me think
                        prev_boundary_v_l = r;

                        T temp_r = prev_val_till_v_r;
                        
                        prev_val_till_v_r = func(prev_val_till_v_l, rval);
                        prev_val_till_v_l = temp_r;
                    
                    } 
                    else if (prev_boundary_w_l == r) 
                    {
                        prev_boundary_w_l = r; // Not necessary but helps me think
                        prev_boundary_v_r = l;

                        T temp_l = prev_val_till_v_l;
                        
                        prev_val_till_v_l = func(prev_val_till_v_r, lval);
                        prev_val_till_v_r = temp_l;
                    } 
                    else // prev_boundary_w_r == r 
                    { 
                        prev_boundary_w_r = r; // Not necessary but helps me think
                        prev_val_till_w_r = prev_val_till_w_r;
                        
                        prev_boundary_w_l = l;
                        prev_val_till_w_l = func(prev_val_till_w_l, rval);
                    }
                }
            }

            w = base_clusters[w].get_parent();
        }

        
    }

    if(prev_boundary_v_l == defretval && prev_boundary_w_l == defretval)
    {
        return defretval;
    }
    if(prev_boundary_v_l == defretval)
    {
        if(w == prev_boundary_w_r)
            return prev_val_till_w_r;
        else
            return prev_val_till_w_l;   
    }
    if(prev_boundary_w_l == defretval)
    {
        if(v == prev_boundary_v_l)
            return prev_val_till_v_l;
        else
            return prev_val_till_v_r;
    }
    
    T v_contrib, w_contrib;
    
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

template<typename T, typename assocfunc>
T subtreeQuery(T root, T u, T defretval, parlay::sequence<cluster<T>> all_clusters, assocfunc func)
{
    // rewrite this, I will have to think about this form the POV of clusters expanding,
    // NOT RC trees!

    return (T) 0;
}



template<typename T>
void printTree(parlay::sequence<T> cluster_indices, parlay::sequence<cluster<T> > &base_clusters)
{
    for(uint i = 0; i < cluster_indices.size(); i++)
    {
        base_clusters[cluster_indices[i]].print(base_clusters);
    }   
}

template<typename T>
void printTree(parlay::sequence<cluster<T> > &base_clusters)
{
    for(uint i = 0; i < base_clusters.size(); i++)
    {
        base_clusters[i].print(base_clusters);
    }   
}