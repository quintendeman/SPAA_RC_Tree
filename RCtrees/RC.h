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
#include "cluster.h"


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

        static const double anywhere_left_weight = 15;
        static const double immediate_left_weight = 30;
        static const double root_weight = 0.0; /* warning, does not represet probability of a null cluster as degree capping may create more forests */

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
        if(v == parents[v])
            return;
        T parent_count = counts[parents[v]].fetch_add(1);
        if(parent_count < (max_degree - 1))
        {
            return;
        }
        else
            parents[v] = v;
    });

}

template<typename T, typename D>
void degree_cap_add_edge(parlay::sequence<T> &parents, const T max_degree, parlay::sequence<std::tuple<T, T, D> >& tuples)
{
    parlay::sequence<std::atomic<T>> counts = parlay::tabulate(parents.size(), [] (size_t) {
       return std::atomic<T>(0); // Initialize each element with the value 0
    });

    std::vector<std::mutex> mutexes(counts.size());

    parlay::parallel_for(0, counts.size(), [&] (T chld) {
        if(parents[chld] != chld)
        {
            counts[parents[chld]].fetch_add(1);
        }
    });

    tuples = parlay::filter(tuples, [&] (auto tple) {
        auto& child = std::get<0>(tple);
        auto& parent = std::get<1>(tple);
        if(child == parent)
            return false;
        bool ret_val = false;

        mutexes[parent].lock();
        mutexes[child].lock();
        if(counts[child] < max_neighbours)
        {
            if(counts[parent] < max_neighbours - 2)
            {
                counts[parent]++;
                parents[child] = parent;
                ret_val = true;
            }
            else if(counts[parent] == max_neighbours - 1 && parents[parent] == parent)
            {
                counts[parent]++;
                parents[child] = parent;
                ret_val = true;        
            }
        }
        mutexes[parent].unlock();
        mutexes[child].unlock();
        return ret_val;
    });
    

    return;
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
void colour_nodes(parlay::sequence<node<T,D>*> tree_nodes, parlay::sequence<cluster<T,D>>& clusters)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, tree_nodes.size(), [&] (T I) {
        auto& cluster_ptr = clusters[tree_nodes[I]->index];
        unsigned long local_maximum = cluster_ptr->get_default_colour();
        unsigned long local_minimum = local_maximum;
        unsigned long my_colour = local_maximum;

        auto& node = *tree_nodes[I];

        for(uint i = 0; i < tree_nodes[I]->size(); i++)
        {
            // get other side
            const auto& node_ptr = node[i];
            if(node_ptr == nullptr || node_ptr->state & (binary_cluster | base_edge))
                continue;
            cluster<T,D>* other_ptr = nullptr;
            
            for(const auto& pot_other_ptr : node->adjacents)
                if(pot_other_ptr != nullptr && pot_other_ptr->index != cluster_ptr->index)
                    other_ptr = clusters[pot_other_ptr->index];
            if(other_ptr == nullptr)
                continue;

            unsigned long compared_colour = other_ptr->get_default_colour();
            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                local_minimum = compared_colour;   
        }
        if(local_maximum == my_colour) // This node is a local maximum, give it a unique colour
        {
            cluster_ptr->colour = local_maximum_colour;
        }
        else if(local_minimum == my_colour)
        {
            cluster_ptr->colour = local_minimum_colour;
        }
        else
        {
            cluster_ptr->colour = 2 + (get_single_colour_contribution(my_colour, local_maximum) / 2); // adding 2 and removing indicator bit
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
void set_MIS(parlay::sequence<node<T,D>*> tree_nodes, parlay::sequence<cluster<T, D> > &base_clusters, bool randomized = false)
{
    colour_nodes(tree_nodes, base_clusters);

}


/**
 * Given an symmetric graph, creates a set of clusters
 * In total, it creates n + m clusters in the array base_clusters
 * The first n are base_vertex clusters
 * it returns an array of n clusters corresponding to the original n vertices.
 * The edge clusters are present as edges.
*/
template <typename T, typename D>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T, D> > &base_clusters, const T max_size, D defretval = 0.00f)
{

    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    base_clusters = parlay::tabulate(G.size(), [] (T i) {
        cluster<T,D> base_cluster;
        return base_cluster;
    });

    parlay::parallel_for(0, base_clusters.size(), [&] (T i){
        base_clusters[i].index = i;
        base_clusters[i].add_empty_level(live, 0);
    });


    parlay::parallel_for(0, base_clusters.size(), [&] (const T i){
        auto& clstr = base_clusters[i];
        clstr.state = base_vertex | live;
        
        const auto& v= i;

        for(const auto& w : G[v])
        {
            if(w < v)
                continue;
            auto edge_cluster = cluster_allocator::alloc();
            edge_cluster->index = -1;
            edge_cluster->state = base_edge | live;
            edge_cluster->data = defretval;
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
    
}





template <typename T, typename D, typename assocfunc>
D PathQuery( cluster<T, D>* v,  cluster<T, D>* w, const D& defretval, assocfunc func)
{
    return defretval;
}

/**
 * Given a pair [v,w], returns the  (if it exists)
 * Exploits the fact that an edge [v,w] must be the child of v or w in an RC tree
*/
template<typename T, typename D>
cluster<T, D>* getEdge(const T v, const T w, parlay::sequence<cluster<T, D>>& clusters)
{

    if(v == w)
        return nullptr;

    return nullptr;
}

template<typename T, typename D, typename assocfunc>
void batchModifyEdgeWeights(const parlay::sequence<std::tuple<T, T, D>>& edges, assocfunc func, parlay::sequence<cluster<T, D>>& clusters, const D& defretval = 0.00f)
{   

    return;
}



/**
 * Only frees the "floating" edge clusters
*/
template<typename T, typename D>
void deleteRCtree(parlay::sequence<cluster<T, D>> &base_clusters)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

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