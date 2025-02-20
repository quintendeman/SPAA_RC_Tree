#ifndef MIS_H
#define MIS_H
#include "parlay/sequence.h"

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
void colour_nodes(parlay::sequence<node<T,D>*> tree_nodes)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, tree_nodes.size(), [&] (T I) {
        auto& cluster_ptr = tree_nodes[I]->cluster_ptr;
        unsigned long local_maximum = cluster_ptr->get_default_colour();
        unsigned long local_minimum = local_maximum;
        unsigned long my_colour = local_maximum;

        auto& node_ptr = tree_nodes[I];

        for(auto& edge_ptr : node_ptr->adjacents)
        {
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            cluster<T,D>* other_ptr = get_other_side(node_ptr, edge_ptr)->cluster_ptr;
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
void set_MIS(parlay::sequence<node<T, D>*>& tree_nodes, bool use_tree_nodes = false, bool randomized = true)
{
    
    if(!randomized)
    {

        auto cluster_ptrs = parlay::delayed_tabulate(tree_nodes.size(), [&] (T i) {
            return tree_nodes[i]->cluster_ptr;
        });
        colour_nodes(tree_nodes);

        parlay::parallel_for(0, cluster_ptrs.size(), [&] (T v) {
            cluster_ptrs[v]->state &= (~IS_MIS_SET);
            if(use_tree_nodes)
                cluster_ptrs[v]->set_neighbour_mis(false, tree_nodes[v]);
            else
                cluster_ptrs[v]->set_neighbour_mis(false);
        });

        auto colours = parlay::delayed_tabulate(cluster_ptrs.size(), [&] (T v) {
            return cluster_ptrs[v]->colour.load();
        });

        auto vertices = parlay::tabulate(cluster_ptrs.size(), [] (T v) {
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
                if(use_tree_nodes)
                {
                    if(cluster_ptrs[v]->get_neighbour_mis(tree_nodes[v]) == true)
                    {
                        cluster_ptrs[v]->state &= (~IS_MIS_SET);
                        return;
                    } 
                }
                else if(cluster_ptrs[v]->get_neighbour_mis() == true)
                {
                    cluster_ptrs[v]->state &= (~IS_MIS_SET);
                    return;
                }
                cluster_ptrs[v]->state |= IS_MIS_SET;
            });
        }
    }
    else // randomized
    {
        //
        auto cluster_ptrs = parlay::delayed_tabulate(tree_nodes.size(), [&] (T i) {
            return tree_nodes[i]->cluster_ptr;
        });
        
        parlay::random_generator gen;
        std::uniform_int_distribution<T> dis(1, std::numeric_limits<T>::max()-1); //changed to max to account for variable types
        

        parlay::parallel_for(0, cluster_ptrs.size(), [&] (T v) {
            cluster_ptrs[v]->state &= (~IS_MIS_SET);
            auto r = gen[v];
            cluster_ptrs[v]->colour = dis(r);
            // if(use_tree_nodes)
            // {
            //     // auto& my_node = tree_nodes[v];
            //     // for(auto& edge_ptr : my_node->adjacents)
            //     // {
            //     //     if(edge_ptr == nullptr)
            //     //         continue;
            //     //     if(edge_ptr->state & (base_edge | binary_cluster))
            //     //     {
            //     //         auto ngbr_node = get_other_side(my_node, edge_ptr);
            //     //         auto r = gen[v];
            //     //         ngbr_node->colour = static_cast<T>(0);
            //     //     }
            //     // }
            // }
            // else
            // {
            //     // auto& my_node = cluster_ptrs[v]->adjacency.get_tail();
            //     // for(auto& edge_ptr : my_node->adjacents)
            //     // {
            //     //     if(edge_ptr == nullptr)
            //     //         continue;
            //     //     if(edge_ptr->state & (base_edge | binary_cluster))
            //     //     {
            //     //         auto ngbr_node = get_other_side(my_node, edge_ptr);
            //     //         ngbr_node->colour = static_cast<T>(0);
            //     //     }
            //     // }
            // }
        });

        parlay::parallel_for(0, cluster_ptrs.size(), [&] (T v) {
            bool in_mis = true;
            T my_colour = cluster_ptrs[v]->colour.load();
            if(use_tree_nodes)
            {
                auto& my_node = tree_nodes[v];
                for(auto& edge_ptr : my_node->adjacents)
                {
                    if(edge_ptr == nullptr)
                        continue;
                    if(edge_ptr->state & (base_edge | binary_cluster))
                    {
                        auto ngbr_node = get_other_side(my_node, edge_ptr);
                        if(ngbr_node->cluster_ptr->colour > my_colour)
                            in_mis = false;
                    }
                }
            }
            else
            {
                auto my_node = cluster_ptrs[v]->adjacency.get_tail();
                for(auto& edge_ptr : my_node->adjacents)
                {
                    if(edge_ptr == nullptr)
                        continue;
                    if(edge_ptr->state & (base_edge | binary_cluster))
                    {
                        auto ngbr_node = get_other_side(my_node, edge_ptr);
                        if(ngbr_node->cluster_ptr->colour >= my_colour)
                            in_mis = false;
                    }
                }
            }
            if(in_mis)
                cluster_ptrs[v]->state |= IS_MIS_SET;
        });


    }

}


template <typename T, typename D>
void check_mis(parlay::sequence<node<T, D>*>& tree_nodes)
{
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        for(auto& ptr : node_ptr->adjacents)
        {
            if(ptr == nullptr || !(ptr->state & (base_edge | binary_cluster)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, ptr);
            if(other_node_ptr->cluster_ptr->state & IS_MIS_SET && node_ptr->cluster_ptr->state & IS_MIS_SET )
            {
                std::cout << red << "Neighbours are not MIS!" << reset << std::endl;
                exit(1);
            }
        }
    });

    return;
}

#endif