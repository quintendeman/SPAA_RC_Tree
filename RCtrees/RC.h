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
#include "../examples/counting_sort.h"
#include "utils.h"

static const char PRINT_QUERY = 0;

template <typename edgetype>
struct edge_with_flag
{
    edgetype E;
    bool valid = true;
};


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
        std::uniform_int_distribution<T> dis(1, 1l << 60l); // range should be enough
        

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
    node_ptr->cluster_ptr->colour = static_cast<T>(0);
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

        for (short k = 0; k < left_node->adjacents.size(); k++)
        {
            auto& ptr = left_node->adjacents[k];
            if (ptr == left_edge)
            {
                ptr = node_ptr;
            }
        }

        for (short k = 0; k < right_node->adjacents.size(); k++)
        {
            auto& ptr = right_node->adjacents[k];
            if (ptr == right_edge)
            {
                ptr = node_ptr;
            }
        }


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
    // Don't! it crashes the program :)
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
            auto edge_cluster = cluster_allocator::create();
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

            auto edge_cluster = cluster_allocator::create();
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
            assert(false && "Shouldn't reach here, no joining edge founf");
        });

        edges_wflag = parlay::filter(edges_wflag, [] (auto ewf) {return ewf.valid;});
    }

    // printTree(base_clusters);

    return;

}

/**
 * Make sure the base clusters are consistent i.e. that every neighbour of mine has me as their neighbour too
 */
template <typename T, typename D>
void check_consistency(parlay::sequence<node<T, D>*>& tree_nodes)
{
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        if((node_ptr->state & live) == 0)
        {
            return;
        }
        if(node_ptr->state & (unary_cluster | binary_cluster | base_edge))
        {
            return;
        }

        for(auto& ptr : node_ptr->adjacents)
        {
            if(ptr == nullptr || !(ptr->state & (base_edge | binary_cluster)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, ptr);
            if(get_other_side(other_node_ptr, ptr) != node_ptr)
            {
                node_ptr->cluster_ptr->print();
                ptr->cluster_ptr->print();
                other_node_ptr->cluster_ptr->print();
                std::cout << "Base clusters inconsistent" << std::endl;
                exit(1);
            }
        }
    });

    return;
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

template <typename T, typename D>
void check_parents_children(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        auto cluster_ptr = &clusters[i];
        if(cluster_ptr->parent != nullptr)
        {
            bool in_parent = false;
            for(auto& par_child : cluster_ptr->parent.load()->children)
            {
                if(par_child == cluster_ptr)
                    in_parent = true;
            }
            if(in_parent == false)
            {
                std::cout << red << "RC tree parents inconsistent!" << reset << std::endl;
                cluster_ptr->print();
                cluster_ptr->parent.load()->print();
                exit(1);
            }
        }
        for(auto& child : cluster_ptr->children)
        {
            if(child == nullptr)
                continue;
            if(child->parent != cluster_ptr)
            {
                std::cout << red << "RC tree children inconsistent!" << reset << std::endl;
                cluster_ptr->print();
                child->print();
                exit(1);
            }
        }
    });
}

template<typename T, typename D>
void check_children_values(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        auto cluster_ptr = &clusters[i];
        D value_from_children = 0.0;
        for(auto& child : cluster_ptr->children)
        {
            if(child == nullptr)
                continue;
            value_from_children+=child->data;
        }
        if(value_from_children != cluster_ptr->data)
        {
            std::cout << red << "RC tree values inconsistent!" << reset << std::endl;
            cluster_ptr->print();
            exit(1);
        }
    });
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
void create_RC_tree(parlay::sequence<cluster<T,D> > &base_clusters, T n, D defretval, lambdafunc assocfunc, bool randomized = false)
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
        // check_consistency(tree_nodes); // TODO remove
        // break;
        auto eligible = parlay::filter(tree_nodes, [] (auto node_ptr){
            return node_ptr->get_num_neighbours_live() <= 2;
        });

        
        set_MIS(eligible, false, randomized);
        auto candidates = parlay::filter(eligible, [] (auto node_ptr){
            return node_ptr->cluster_ptr->state & IS_MIS_SET;
        });

        // check_mis(candidates); // TODO remove

        
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
    
        
        // std::cout << "[" << count  << "]: " << bold << tree_nodes.size() << reset << std::endl;
        count++;
    }while(tree_nodes.size() > 0 && count < 100);

    assert(count < 100 && "Did your tree go into an infinite loop?");

}

template<typename T, typename D>
node<T,D>* get_edge(T v, T w, parlay::sequence<cluster<T,D>>& clusters, unsigned char level = 0)
{
    if(v == w)
        return nullptr;
    auto& v_cluster = clusters[v];
    auto v_node_ptr = v_cluster.adjacency.get_head();
    for(auto& edge_ptr : v_node_ptr->adjacents)
    {
        if(edge_ptr == nullptr)
            continue;
        auto w_node_ptr = get_other_side(v_node_ptr, edge_ptr);
        if(w_node_ptr->cluster_ptr->index == w)
            return edge_ptr;
    }

    auto& w_cluster = clusters[w];
    auto w_node_ptr = w_cluster.adjacency.get_head();
    for(auto& edge_ptr : w_node_ptr->adjacents)
    {
        if(edge_ptr == nullptr)
            continue;
        v_node_ptr = get_other_side(w_node_ptr, edge_ptr);
        if(v_node_ptr->cluster_ptr->index == v)
            return edge_ptr;
    }
    return nullptr;
}




template<typename T, typename D, typename assocfunc>
void batchModifyEdgeWeights(const parlay::sequence<std::tuple<T, T, D>>& edges, assocfunc func, parlay::sequence<cluster<T, D>>& clusters, const D& defretval = 0.00f)
{   

    // TODO merge  with up pass
    auto edge_ptrs = parlay::tabulate(edges.size(), [&] (T i){
        T v = std::get<0>(edges[i]);
        T w = std::get<1>(edges[i]);
        D data = std::get<2>(edges[i]);

        auto edge_ptr = get_edge(v, w, clusters);
        // if(edge_ptr == nullptr) // TODO REMOVE
        // {
        //     std::cout << red << "Got a nullptr" << reset << std::endl;
        //     exit(1);
        // }
        edge_ptr->cluster_ptr->data = data;

        auto ret_ptr = edge_ptr->cluster_ptr;
        auto cluster_ptr = ret_ptr;

        cluster_ptr->counter = 0;

        while(cluster_ptr != nullptr)
        {
            if(cluster_ptr->counter.fetch_add(1))
                break;
            cluster_ptr = cluster_ptr->get_parent();
        }


        return ret_ptr;
    });

    parlay::parallel_for(0, edge_ptrs.size(), [&] (T i) {
        if(edge_ptrs[i] == nullptr)
            return;

        auto cluster_ptr = edge_ptrs[i]; // base edge, technically not a cluster dw.. depends on which paper you look at haha
        
        auto start = edge_ptrs[i];

        while(cluster_ptr != nullptr)
        {
            if(cluster_ptr->counter.fetch_add(-1) > 1)
                break;
            if(cluster_ptr != start)
            {
                D final_val = defretval; 
                cluster<T,D>* final_heaviest_edge = nullptr;
                node<T,D>* node_ptr = cluster_ptr->first_contracted_node;
                bool first_child = true;

                for(auto& child : cluster_ptr->children)
                {
                    if(child == nullptr)
                        continue;
                    if((child->adjacency.get_head()->state & base_edge) || (child->first_contracted_node->next->state & (binary_cluster | base_edge)))
                    {
                        if(first_child)
                        {
                            final_val = child->data;
                            final_heaviest_edge = child->max_weight_edge;
                            first_child = false;
                        }
                        else
                        {
                            if(child->data > final_val)
                                final_heaviest_edge = child->max_weight_edge;
                            final_val = func(final_val, child->data);
                            
                        }
                    }
                }              
                cluster_ptr->data = final_val;
                cluster_ptr->max_weight_edge = final_heaviest_edge;
            }
            cluster_ptr = cluster_ptr->get_parent();
        }
    });

    return;
}

template<typename T, typename D>
void check_counter(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        if(clusters[i].counter)
        {
            std::cout << red << "Counter wasn't reset" << reset << std::endl;
            exit(1);
        }
    });
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

            D lval_print,rval_print = defretval; // TODO remove

            // We are at the start probably can add a bunch of nullptr checks too TODO
            if (prev_val_till_v_l == defretval && prev_val_till_v_r == defretval) 
            {
                v_cluster->find_boundary_vertices(prev_boundary_v_l, prev_val_till_v_l, prev_boundary_v_r, prev_val_till_v_r, defretval);
                lval_print = prev_val_till_v_l;
                rval_print = prev_val_till_v_r;
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
                    lval_print = prev_val_till_w_l;
                    rval_print = prev_val_till_w_r;
                    
                    if(l != r) // ascending into unary
                    {
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                        
                    }
                    else
                    {
                        auto v_node_ptr = v_cluster->first_contracted_node;
                        for(auto& ptr : v_node_ptr->adjacents)
                        {
                            if(ptr != nullptr && ptr->state & (binary_cluster | base_edge))
                            {
                                prev_val_till_v_l = func(prev_val_till_v_l, ptr->cluster_ptr->data);
                                prev_val_till_v_r = func(prev_val_till_v_r, ptr->cluster_ptr->data);

                            }
                        }
                    }
                    prev_boundary_v_l = l;
                    prev_boundary_v_r = r;
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
                lval_print = lval;
                rval_print = rval;
            }
            // std::cout << bright_magenta << v->index << "[" << (short) v->get_height() << "]" << " ";
            // std::cout << prev_boundary_v_l->index << "(" << prev_val_till_v_l << ") " << prev_boundary_v_r->index << "(" << prev_val_till_v_r << ") ";
            // std::cout << "lval,rval=" << lval_print << "," << rval_print << " ";
            // auto& for_printing = v;
            // if(for_printing->adjacency.get_tail()->state & unary_cluster)
            //     std::cout << " unary ";
            // else if (for_printing->adjacency.get_tail()->state & binary_cluster)
            //     std::cout << " binary ";
            // else
            //     std::cout << " nullary ";
            // std::cout << reset << std::endl;
            v = v->get_parent();
        }
        else
        {
            cluster<T,D>* &w_cluster = w;
            D lval_print,rval_print = defretval;

            // We are at the start
            if (prev_val_till_w_l == defretval && prev_val_till_w_r == defretval) 
            {
                w_cluster->find_boundary_vertices(prev_boundary_w_l, prev_val_till_w_l, prev_boundary_w_r, prev_val_till_w_r, defretval);
                lval_print = prev_val_till_w_l;
                rval_print = prev_val_till_w_r;
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

                    if(l != r)
                    {
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    }
                    else
                    {
                        auto w_node_ptr = w_cluster->first_contracted_node;
                        for(auto& ptr : w_node_ptr->adjacents)
                        {
                            if(ptr != nullptr && ptr->state & (binary_cluster | base_edge))
                            {
                                prev_val_till_w_l = func(prev_val_till_w_l, ptr->cluster_ptr->data);
                                prev_val_till_w_r = func(prev_val_till_w_r, ptr->cluster_ptr->data);

                            }
                        }
                    }
                    prev_boundary_w_l = l;
                    prev_boundary_w_r = r;
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
                lval_print = lval;
                rval_print = rval;
            }
            // std::cout << bright_green << w->index << "[" << (short) w->get_height() << "]" << " ";
            // std::cout << prev_boundary_w_l->index << "(" << prev_val_till_w_l << ") " << prev_boundary_w_r->index << "(" << prev_val_till_w_r << ") ";
            // std::cout << "lval,rval=" << lval_print << "," << rval_print << " ";
            // auto& for_printing = w;
            // if(for_printing->adjacency.get_tail()->state & unary_cluster)
            //     std::cout << " unary ";
            // else if (for_printing->adjacency.get_tail()->state & binary_cluster)
            //     std::cout << " binary ";
            // else
            //     std::cout << " nullary ";
            // std::cout << reset << std::endl;
            w = w->get_parent();
        }
        
    }
    // std::cout << std::endl;

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

template<typename T, typename D>
void testPathQueryValid(parlay::sequence<cluster<T,D>>& clusters, parlay::sequence<T>& parents, parlay::sequence<D>& weights, parlay::sequence<T>& random_perm_map, long base_size)
{
    parlay::parallel_for(0, 100000, [&] (T i) {
        T start_index = rand() % base_size;
        T index = start_index;
        T end_index;
        // doing sum!
        D manual_result = 0.0;
        while(parents[index] != index)
        {
            manual_result = manual_result + weights[index];
            index = parents[index];
        }
        end_index = index;

        D pathQueryResult = PathQuery(&clusters[random_perm_map[start_index]], &clusters[random_perm_map[end_index]], (D) 0.0, [] (D a, D b) {return a + b;});

        auto almost_equal = [] (D A, D B) -> bool {
        const double tolerance = 0.001;
            if(A == B)
                return true;
            if(A > B)
                if((A-B) < tolerance)
                    return true;
            if(B > A)
                if((B-A) < tolerance)
                    return true;
            return false;
        };

        if(!almost_equal(pathQueryResult, manual_result))
        {
            std::cout << "Path query from " << start_index << " to " << end_index << " invalid! " << std::endl;
            std::cout << pathQueryResult << " != " << manual_result << std::endl;
            assert("RC Path query and manual path query should be same" && false);
        }
       
    });

}


#endif