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
        auto random_val = dis(gen);
        if (random_val <= 0.9 && v > 0)
        {
            std::uniform_int_distribution<> disint(0, v-1);

            dummy_initial_parents[v] = disint(gen);
        }
        else if (random_val < 0.905)
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
void colour_clusters(parlay::sequence<cluster<T>*> clusters)
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
template<typename T>
void set_MIS(parlay::sequence<cluster<T>*> clusters, bool randomized = false)
{
    if(!randomized)
    {
        colour_clusters(clusters);

        parlay::parallel_for(0, clusters.size(), [&] (T v) {
            clusters[v]->is_MIS = false;
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
                    clusters[v]->is_MIS = false;
                    return;
                }
                clusters[v]->is_MIS = true;
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
            clusters[v]->is_MIS=clusters[v]->is_max_neighbour_colour();
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
        auto cluster = &base_clusters[v];
        cluster->index = v;
        cluster->state = base_vertex | live;

        for(short i = 0; i < G[v].size(); i++) //TODO: Assumes ordering i.e. root is towards zero
        {
            if(G[v][i] < v)
            {
                base_clusters[v+n].add_initial_neighbours(&base_clusters[G[v][i]], cluster);
                base_clusters[v+n].state = base_edge;
                base_clusters[v+n].index = v+n;
                base_clusters[v].add_neighbour(&base_clusters[G[v][i]], &base_clusters[v+n]);
            }
            else
            {
                base_clusters[v].add_neighbour(&base_clusters[G[v][i]], &base_clusters[G[v][i]+n]);
            }
        }

    });

}

template <typename T>
void set_heights(parlay::sequence<cluster<T>*> &all_cluster_ptrs)
{

    parlay::parallel_for(0, all_cluster_ptrs.size(), [&] (T v) {
        auto cluster = all_cluster_ptrs[v];
        T my_height = 0;

        while(cluster != nullptr)
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
            cluster = cluster->get_parent();
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
void create_RC_tree(parlay::sequence<cluster<T> > &base_clusters, T n, bool do_height = true, bool randomized = false)
{
    // std::cout << "create RC tree called" << std::endl;
    
    parlay::sequence<cluster<T>*> all_cluster_ptrs = parlay::tabulate(n, [&] (T v) {
        return &base_clusters[v];
    });

    parlay::sequence<cluster<T>*> forest, candidates;

    // Initially the forest of live nodes is all live nodes
    forest = parlay::filter(all_cluster_ptrs, [&] (cluster<T>* C) {
        return ((C->state & live));
    });

    bool first_time = true;

    do
    {
        // Shrink the forst as live nodes decrease
        if(first_time == true)
        {
            
            first_time = false;
        }
        else
        {
            forest = parlay::filter(forest, [&] (cluster<T>* C) {
                return ((C->state & live));
            });
        }
        
        // std::cout << "forest.size(): " << forest.size() << std::endl;

        // Eligible nodes are those with 0, 1 or 2 neighbours
        auto eligible = parlay::filter(forest, [&] (cluster<T>* C) {
            return (C->get_neighbour_count() <= 2);
        });

        // std::cout << "eligible.size(): " << eligible.size() << std::endl;

        // Set the flag is_MIS amongst them
        set_MIS(eligible, randomized = randomized);

        // Filter out an MIS of eligible nodes
        candidates = parlay::filter(eligible, [&] (cluster<T>* C) {
            return (C->is_MIS);
        });

        // std::cout << "candidates.size(): " << candidates.size() << std::endl;

        // printTree(base_clusters);

        // do rake and compress
        parlay::parallel_for(0, candidates.size(), [&] (T v) {
            cluster<T>* cluster_ptr = candidates[v];
            // rake
            if(cluster_ptr->get_neighbour_count() == 0)
            {
                cluster_ptr->state&=(~live);
                cluster_ptr->state|=(nullary_cluster);
                cluster_ptr->state|=internal;
            }
            if(cluster_ptr->get_neighbour_count() == 1)
            {
                cluster<T>* edge_ptr = nullptr;
                cluster<T>* other_side = nullptr;
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
                other_side->change_to_child(edge_ptr);
                other_side->change_to_child(cluster_ptr);

                // // make other_side be the parent for both
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
                cluster<T>* left_edge_ptr = nullptr;
                cluster<T>* right_edge_ptr = nullptr;
                cluster<T>* left_node_ptr = nullptr;
                cluster<T>* right_node_ptr = nullptr;

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
                
            }
        });
    }while(candidates.size());

    if(do_height == true)
        set_heights(all_cluster_ptrs);

}




template<typename T>
void printTree(parlay::sequence<cluster<T> > &base_clusters)
{
    for(uint i = 0; i < base_clusters.size(); i++)
    {
        base_clusters[i].print();
    }   
}