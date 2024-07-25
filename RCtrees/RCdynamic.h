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


template<typename T, typename D>
bool first_condition(node<T,D>*& node_ptr)
{
    return node_ptr->state & affected;
}

template<typename T, typename D>
bool second_condition(node<T,D>*& node_ptr)
{
    return node_ptr->state & adjacency_changed;
}

template<typename T, typename D>
bool third_condition(node<T,D>*& node_ptr)
{
    if(node_ptr->cluster_ptr->first_contracted_node == node_ptr)
        return false;
    // check neighbours
    for(auto& edge_ptr : node_ptr->adjacents)
    {
        if(edge_ptr == nullptr)
            continue;
        if(edge_ptr->state & (base_edge | binary_cluster))
        {
            auto other_node = get_other_side(node_ptr, edge_ptr);
            if(!(first_condition(other_node) || second_condition(other_node)))
                if(other_node->cluster_ptr->first_contracted_node == other_node)
                    return false;
        }
    }

    return true;
}

template<typename T, typename D>
parlay::sequence<node<T,D>*> get_3dp1(parlay::sequence<node<T,D>*>& initial_nodes)
{
    parlay::sequence<node<T,D>*> frontier = parlay::flatten(parlay::tabulate(initial_nodes.size(), [&] (T I) {
        auto& node_ptr = initial_nodes[I];
        parlay::sequence<node<T,D>*> ret_seq = parlay::sequence<node<T,D>*>(max_neighbours, nullptr);
        for(short i = 0; i < max_neighbours; i++)
        {
            if(node_ptr->adjacents[i] == nullptr || node_ptr->adjacents[i]->state & (binary_cluster | base_edge) == 0)
                continue;
            ret_seq[i] = get_other_side(node_ptr, node_ptr->adjacents[i]);
        }
        return ret_seq;
    }));

    frontier = parlay::flatten(parlay::tabulate(frontier.size(), [&] (T I) {
        auto& node_ptr = frontier[I];
        parlay::sequence<node<T,D>*> ret_seq = parlay::sequence<node<T,D>*>(max_neighbours + 1, nullptr);
        if(node_ptr != nullptr)
            for(short i = 0; i < max_neighbours; i++)
            {
                if(node_ptr->adjacents[i] == nullptr || node_ptr->adjacents[i]->state & (binary_cluster | base_edge) == 0)
                    continue;
                ret_seq[i] = get_other_side(node_ptr, node_ptr->adjacents[i]);
            }
        ret_seq[max_neighbours] = node_ptr;
        return ret_seq;
    }));

    return tiebreak(frontier);
}

template<typename T, typename D>
bool is_update_eligible(node<T,D>* node_ptr)
{
    for(auto& edge_ptr : node_ptr->adjacents)
    {
        if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster|base_edge)))
            continue;
        auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
        if((other_node_ptr->state & affected) == 0 && other_node_ptr->cluster_ptr->first_contracted_node == other_node_ptr)
            return false;
    }

    return true;
}

/** 
 * Filter and remove redundant pointers and nullptrs
*/

template<typename T, typename D>
parlay::sequence<node<T,D>*> tiebreak(parlay::sequence<node<T,D>*>& input_nodes)
{
    parlay::parallel_for(0, input_nodes.size(), [&] (T i){
        auto& node_ptr = input_nodes[i];
        if(node_ptr  == nullptr)
            return;
        node_ptr->cluster_ptr->tiebreak = i;
    });
    parlay::parallel_for(0, input_nodes.size(), [&] (T i){
        auto& node_ptr = input_nodes[i];
        if(node_ptr  != nullptr && node_ptr->cluster_ptr->tiebreak != i)
            node_ptr = nullptr;
    });
    return parlay::filter(input_nodes, [] (auto node_ptr) {return node_ptr != nullptr;});
}

// does essentially the same thing as recreate last levels but with additional checks if a node has already been allocated in a new layer or not
template<typename T, typename D>
void create_decompressed_affected(parlay::sequence<node<T,D>*>& affected_nodes)
{
    parlay::parallel_for(0, affected_nodes.size(), [&] (T i) {
        auto& aff_node = affected_nodes[i];
        if(aff_node->next == nullptr)
            aff_node->cluster_ptr->add_empty_level(aff_node->state & (~(binary_cluster | unary_cluster | nullary_cluster)) | live, aff_node->contraction_level+1);
    });

    parlay::parallel_for(0, affected_nodes.size(), [&] (T i) {
        auto& aff_node = affected_nodes[i];
        if(aff_node->next == nullptr)
            aff_node->cluster_ptr->add_empty_level(aff_node->state & (~(binary_cluster | unary_cluster | nullary_cluster)) | live, aff_node->contraction_level+1);
    });

    // decontract rakes
    parlay::parallel_for(0, affected_nodes.size(), [&] (T i){
        auto& aff_node = affected_nodes[i];
        if(!is_update_eligible(aff_node)) // if it isn't update eligible, it must be alive in the next round
            return;
        if(aff_node->next == nullptr)
            aff_node->cluster_ptr->add_empty_level(aff_node->state, aff_node->contraction_level+1);

        // for each of the vertices, check if they still exist
        // if it raked i.e. the other side exists, revert that
        for (short n = 0; n < aff_node->adjacents.size(); ++n) 
        {
            auto Li_edge_ptr = aff_node->adjacents[n];
            if (Li_edge_ptr == nullptr || !(Li_edge_ptr->state & (base_edge | binary_cluster))) {
                continue;
            }
            auto Li_neighbour_node = get_other_side(aff_node, Li_edge_ptr);

            auto current_edge_ptr = Li_edge_ptr->next;
            auto current_neighbour_node = Li_neighbour_node->next;
            auto current_aff_node = aff_node->next;

            current_aff_node->state &= (~(binary_cluster | unary_cluster | nullary_cluster));

            // check if it exists
            for (short m = 0; m < current_neighbour_node->adjacents.size(); ++m) {
                auto cur_edge = current_neighbour_node->adjacents[m];
                if (cur_edge == current_edge_ptr) {
                    break;
                } else if (cur_edge == current_aff_node) { // need to revert on both sides
                    current_aff_node->state |= affected | adjacency_changed;
                    current_neighbour_node->state |= affected | adjacency_changed;
                    current_neighbour_node->adjacents[m] = current_edge_ptr; // we now point to edge
                    for (short p = 0; p < current_aff_node->adjacents.size(); ++p) {
                        if (current_aff_node->adjacents[p] == current_neighbour_node) {
                            current_aff_node->adjacents[p] = current_edge_ptr;
                            break;
                        }
                    }
                }
            }
        }

        
    });

    // if an edge doesn't exist, both sides must be considered affected already 
    parlay::parallel_for(0, affected_nodes.size(), [&] (T i){
        auto& aff_node = affected_nodes[i];
        if(aff_node->next == nullptr)
            aff_node->cluster_ptr->add_empty_level(aff_node->state, aff_node->contraction_level+1);
        aff_node->next->state |= affected;

        // for each of the vertices, check if they still exist
        // if it raked i.e. the other side exists, revert that
        for (short n = 0; n < aff_node->adjacents.size(); ++n) {
            auto Li_edge_ptr = aff_node->adjacents[n];
            if (Li_edge_ptr == nullptr || !(Li_edge_ptr->state & (base_edge | binary_cluster))) {
                continue;
            }
            auto Li_neighbour_node = get_other_side(aff_node, Li_edge_ptr);

            auto current_edge_ptr = Li_edge_ptr->next;
            auto current_neighbour_node = Li_neighbour_node->next;
            auto current_aff_node = aff_node->next;

            bool exists = false;
            // does it exist in the next stage?
            
            for (short m = 0; m < current_aff_node->adjacents.size(); ++m) {
                auto cur_edge = current_aff_node->adjacents[m];
                if (cur_edge != nullptr && (cur_edge == current_edge_ptr || cur_edge == current_neighbour_node)) {
                    exists = true;
                }
            }
            if (!exists) {
                auto& v = aff_node->cluster_ptr->index;
                auto& w = Li_neighbour_node->cluster_ptr->index;
                if (v < w) {
                    if (current_edge_ptr == nullptr) {    
                        Li_edge_ptr->cluster_ptr->add_empty_level(Li_edge_ptr->state, Li_edge_ptr->contraction_level + 1);
                        current_edge_ptr = Li_edge_ptr->next;
                    }
                    current_edge_ptr->state |= dont_finalize;
                    current_edge_ptr->cluster_ptr->asked_to_increment = Li_edge_ptr->contraction_level;
                    current_edge_ptr->adjacents.fill(nullptr);
                    current_edge_ptr->add_ptr(current_neighbour_node);
                    current_edge_ptr->add_ptr(current_aff_node);
                    current_aff_node->add_ptr(current_edge_ptr);
                    current_aff_node->state |= affected | adjacency_changed;
                }
            }
        }
        
    });

    // if an edge doesn't exist, both sides must be considered affected already 
    parlay::parallel_for(0, affected_nodes.size(), [&] (T i){
        auto& aff_node = affected_nodes[i];
        if(aff_node->next == nullptr)
            aff_node->cluster_ptr->add_empty_level(aff_node->state, aff_node->contraction_level+1);

        // for each of the vertices, check if they still exist
        // if it raked i.e. the other side exists, revert that
        for (short n = 0; n < aff_node->adjacents.size(); ++n) {
            auto Li_edge_ptr = aff_node->adjacents[n];
            if (Li_edge_ptr == nullptr || !(Li_edge_ptr->state & (base_edge | binary_cluster))) {
                continue;
            }
            auto Li_neighbour_node = get_other_side(aff_node, Li_edge_ptr);

            auto current_edge_ptr = Li_edge_ptr->next;
            auto current_neighbour_node = Li_neighbour_node->next;
            auto current_aff_node = aff_node->next;

            bool exists = false;
            // does it exist in the next stage?
            for (short m = 0; m < current_aff_node->adjacents.size(); ++m) {
                auto cur_edge = current_aff_node->adjacents[m];
                if (cur_edge != nullptr && (cur_edge == current_edge_ptr || cur_edge == current_neighbour_node)) {
                    exists = true;
                }
            }

            if (!exists) {
                auto& v = aff_node->cluster_ptr->index;
                auto& w = Li_neighbour_node->cluster_ptr->index;
                if (!(v < w)) {
                    current_aff_node->add_ptr(current_edge_ptr);
                }
                current_aff_node->state |= affected | adjacency_changed;
            }
        }
   
    });

    parlay::parallel_for(0, affected_nodes.size(), [&] (T i){
        auto& aff_node = affected_nodes[i];
        auto& current_aff_node = aff_node->next;

        for (short n = 0; n < current_aff_node->adjacents.size(); ++n) {
            auto new_edge = current_aff_node->adjacents[n];
            if (new_edge == nullptr) {
                continue;
            }
            bool existed = false;

            for (short m = 0; m < aff_node->adjacents.size(); ++m) {
                auto old_edge = aff_node->adjacents[m];
                if (old_edge == nullptr) {
                    continue;
                }
                if (old_edge->next == new_edge) {
                    existed = true;
                    break;
                }
                auto old_other_side = get_other_side(aff_node, old_edge);
                if (old_other_side == nullptr) {
                    continue;
                }
                if (old_other_side->next == new_edge) {
                    existed = true;
                    break;
                }
            }
            if (!existed) {
                current_aff_node->adjacents[n] = nullptr;
            }
        }
    });

    parlay::parallel_for(0, affected_nodes.size(), [&] (T i){
        auto& aff_node = affected_nodes[i];
        auto& current_aff_node = aff_node->next;

        for (short n = 0; n < aff_node->adjacents.size(); ++n) {
            auto old_edge = aff_node->adjacents[n];
            if (old_edge == nullptr) {
                continue;
            }
            if (old_edge->state & unary_cluster && old_edge->next == nullptr) {
                old_edge->state |= dont_finalize;
                old_edge->cluster_ptr->add_empty_level(unary_cluster, old_edge->contraction_level + 1);
                old_edge->next->state = unary_cluster;
                old_edge->next->adjacents.fill(nullptr);
                current_aff_node->add_ptr(old_edge->next);
                // std::cout << "Adding unary for " << old_edge->index() << std::endl;
            } else if (old_edge->state & unary_cluster) {
                old_edge->state |= dont_finalize;
                old_edge->next->adjacents.fill(nullptr);
                old_edge->next->state = unary_cluster;
                current_aff_node->add_ptr(old_edge->next);
                // std::cout << "Modifying unary for " << old_edge->index() << " with adjacency size " << (int) old_edge->cluster_ptr->adjacency.size() << std::endl;
                finalize(old_edge->next);
            }
        }
    });

    return;
}

template<typename T, typename D>
void finalize(node<T,D>* contracted_node)
{
    if(contracted_node->state & dont_finalize)
        return;
    contracted_node->cluster_ptr->finalize_time = contracted_node->contraction_level;
    // std::cout << contracted_node->index() << " went from " << (int) contracted_node->cluster_ptr->adjacency.size() << " ";
    while(contracted_node->cluster_ptr->adjacency.get_tail() != contracted_node)
        {
            auto tail = contracted_node->cluster_ptr->adjacency.get_tail();
            if(tail->state & dont_finalize)
                break;

            for(short i1 = 0; i1 < tail->adjacents.size(); i1++)
            {
                auto one_hop_neighbour = tail->adjacents[i1];
                if(one_hop_neighbour == nullptr)
                    continue;
                for(short i2 = 0; i2 < one_hop_neighbour->adjacents.size(); i2++)
                {
                    auto two_hop_neighbour = one_hop_neighbour->adjacents[i2];
                    if(two_hop_neighbour == nullptr)
                        continue;
                    for(short j = 0; j < two_hop_neighbour->adjacents.size(); j++)
                    {
                        auto back_ptr = two_hop_neighbour->adjacents[j];
                        if(back_ptr == tail)
                            two_hop_neighbour->adjacents[j] = nullptr;
                    }
                    
                    if(two_hop_neighbour == tail)
                        one_hop_neighbour->adjacents[i2] = nullptr;
                }

            }
            tail->adjacents.fill(nullptr);
            contracted_node->cluster_ptr->adjacency.delete_tail();
        }
    // std::cout << " to " << (int) contracted_node->cluster_ptr->adjacency.size() << std::endl;
        return;

}

template<typename T, typename D>
void batchInsertEdge( const parlay::sequence<std::pair<T, T>>& delete_edges, const parlay::sequence<std::tuple<T, T, D>>& add_edges, parlay::sequence<cluster<T, D>>& clusters)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;
    auto tree_nodes = parlay::flatten(parlay::tabulate(delete_edges.size() + add_edges.size(), [&] (T i) {
        T v, w;
        if(i < delete_edges.size())
        {
            v = delete_edges[i].first;
            w = delete_edges[i].second;
            // std::cout << v << " -- " << w << " deleted" << std::endl;
        }
        else
        {
            v = std::get<0>(add_edges[i-delete_edges.size()]);
            w = std::get<1>(add_edges[i-delete_edges.size()]);
            // std::cout << v << " -- " << w << " added" << std::endl;
        }
        parlay::sequence<node<T,D>*> ret_seq = {clusters[v].adjacency.get_head(), clusters[w].adjacency.get_head()};
        if(ret_seq[0] == nullptr || ret_seq[1] == nullptr) // TODO remove
        {
            std::cout << red << "This should never happen, nullptr on head?" << reset << std::endl;
            exit(1);
        }
        ret_seq[0]->tiebreak() = i*2;
        ret_seq[1]->tiebreak() = i*2 +1;
        return ret_seq;
    }));

    if(clusters.size() <= 100) // TODO remove
    {
        for(uint i = 0; i < delete_edges.size() + add_edges.size(); i++)
        {
            T v, w;
            if(i < delete_edges.size())
            {
                v = delete_edges[i].first;
                w = delete_edges[i].second;
                std::cout << v << " -- " << w << " deleted" << std::endl;
            }
            else
            {
                v = std::get<0>(add_edges[i-delete_edges.size()]);
                w = std::get<1>(add_edges[i-delete_edges.size()]);
                std::cout << v << " -- " << w << " added" << std::endl;
            }
        }
    }

    parlay::parallel_for(0, tree_nodes.size(), [&] (T i){
        if(tree_nodes[i]->tiebreak() != i)
            tree_nodes[i] = nullptr;
    });




    tree_nodes = parlay::filter(tree_nodes, [] (auto node_ptr){
        if(node_ptr != nullptr)
        {
            node_ptr->state |= affected | adjacency_changed | live;
            node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        }
        return node_ptr != nullptr;
    });

    auto frontier = get_3dp1(tree_nodes);

    parlay::parallel_for(0, delete_edges.size(), [&] (T I){
        auto& v = delete_edges[I].first;
        auto& w = delete_edges[I].second;
        
        cluster<T,D>* cluster_edge_ptr = nullptr;

        auto v_node_ptr = clusters[v].adjacency.get_head();
       
        for (short i = 0; i < v_node_ptr->adjacents.size(); ++i) {
            auto edge_ptr = v_node_ptr->adjacents[i];
            if (edge_ptr != nullptr && edge_ptr->state & base_edge &&
                get_other_side(v_node_ptr, edge_ptr) != nullptr &&
                get_other_side(v_node_ptr, edge_ptr)->cluster_ptr->index == w) {
                if (v_node_ptr == clusters[v].adjacency.get_head()) {
                    cluster_edge_ptr = edge_ptr->cluster_ptr;
                }
            }
        }
        
        while(v_node_ptr != nullptr) {
            for(short i = 0; i < v_node_ptr->adjacents.size(); ++i) {
                auto edge_ptr = v_node_ptr->adjacents[i];
                if(edge_ptr != nullptr && (edge_ptr->cluster_ptr == cluster_edge_ptr || edge_ptr->cluster_ptr == &clusters[w]))
                    v_node_ptr->adjacents[i] = nullptr;
            }
            v_node_ptr = v_node_ptr->next;
        }

        auto w_node_ptr = clusters[w].adjacency.get_head();
        while(w_node_ptr != nullptr) {
            for(short i = 0; i < w_node_ptr->adjacents.size(); ++i) {
                auto edge_ptr = w_node_ptr->adjacents[i];
                if(edge_ptr != nullptr && (edge_ptr->cluster_ptr == cluster_edge_ptr || edge_ptr->cluster_ptr == &clusters[v]))
                    w_node_ptr->adjacents[i] = nullptr;
            }
            w_node_ptr = w_node_ptr->next;
        }
        cluster_allocator::destroy(cluster_edge_ptr);
    });

    auto add_edges_wflags = parlay::tabulate(add_edges.size(), [&] (T i) {
        edge_with_flag<std::tuple<T, T, D>> ret_edge;
        ret_edge.valid = true;
        ret_edge.E = add_edges[i];
        return std::move(ret_edge);
    });

    auto reduce_edges = add_edges_wflags;

    do
    {
        parlay::parallel_for(0, reduce_edges.size(), [&] (T i) {
            auto& edge = reduce_edges[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);
            clusters[v].tiebreak = i + 1;
        });

        parlay::parallel_for(0, reduce_edges.size(), [&] (T i) {
            auto& edge = reduce_edges[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);

            auto v_node_ptr = clusters[v].adjacency.get_head();
            auto w_node_ptr = clusters[w].adjacency.get_head();

            if(clusters[v].tiebreak == (i+1))
            {
                // does it already exist?
                node<T,D>* edge_node = nullptr;
                for(auto& edge_ptr : v_node_ptr->adjacents)
                {
                    if(edge_ptr != nullptr && get_other_side(v_node_ptr, edge_ptr)->cluster_ptr->index == w)
                        edge_node = edge_ptr;
                }
                if(edge_node == nullptr)
                {
                    auto newEdgeClstrPtr = cluster_allocator::alloc();
                    newEdgeClstrPtr->add_empty_level(base_edge, 0);
                    newEdgeClstrPtr->data = std::get<2>(edge.E);
                    newEdgeClstrPtr->index = -1;
                    edge_node = newEdgeClstrPtr->adjacency.get_head();
                    newEdgeClstrPtr->add_ptr_to_highest_level(v_node_ptr);
                    newEdgeClstrPtr->add_ptr_to_highest_level(w_node_ptr);
                    v_node_ptr->add_ptr(edge_node);
                }
                else
                {
                    edge_node->cluster_ptr->data = std::get<2>(edge.E);
                }
                edge.valid = false;
            }

        });
        reduce_edges = parlay::filter(reduce_edges, [] (auto r){
            return r.valid;
        });
    }while(reduce_edges.size() > 0);

    do
    {
        parlay::parallel_for(0, add_edges_wflags.size(), [&] (T i) {
            auto& edge = add_edges_wflags[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);
            clusters[w].tiebreak = i + 1;
        });

        parlay::parallel_for(0, add_edges_wflags.size(), [&] (T i) {
            auto& edge = add_edges_wflags[i];
            auto& v = std::get<0>(edge.E);
            auto& w = std::get<1>(edge.E);
            if(v < w)
                std::swap(v, w);

            auto v_node_ptr = clusters[v].adjacency.get_head();
            auto w_node_ptr = clusters[w].adjacency.get_head();

            if(clusters[w].tiebreak == (i+1))
            {
                // find other node_ptr
                for(auto& potential_edge_ptr : v_node_ptr->adjacents)
                {
                    if(potential_edge_ptr != nullptr && get_other_side(v_node_ptr, potential_edge_ptr) == w_node_ptr)
                    {  
                        w_node_ptr->add_ptr(potential_edge_ptr); 
                        break;
                    }
                }
                edge.valid = false;
            }

        });
        add_edges_wflags = parlay::filter(add_edges_wflags, [] (auto r){
            return r.valid;
        });
    }while(add_edges_wflags.size() > 0);

    
    // auto frontier = get_3dp1(tree_nodes); // NOT to be confused with the frontier in the paper

    std::cout << "endpoints size is " << tree_nodes.size() << std::endl;
    std::cout << "frontier size is " << frontier.size() << std::endl;
    
    frontier = parlay::filter(frontier, [] (auto node_ptr){
        if(first_condition(node_ptr) || second_condition(node_ptr) || third_condition(node_ptr))
        {
            return true;
        }
        return false;
    });
    parlay::parallel_for(0, frontier.size(), [&] (T i) {
        frontier[i]->state |= affected;
    });

    std::cout << "Which gets reduced to " << frontier.size() << std::endl;
    
    // create_decompressed_affected(frontier);

    parlay::sequence<node<T,D>*> mis_set;

    unsigned short count = 0;
    static const short max_count = 10;
    do
    {
        auto checker_nodes = parlay::tabulate(clusters.size(), [&] (T i){
            return clusters[i].adjacency[count];
        });
        check_consistency(checker_nodes);

        if(clusters.size() <= 100)
        {
            std::cout << "frontier: ";
            for(auto& ptr : frontier)
               std::cout <<  ptr->cluster_ptr->index << " ";
            std::cout << std::endl;
            std::cout << "At level " << (int) frontier[0]->contraction_level << std::endl;
            printTree(clusters, count+1);
        }
        auto new_frontier = get_3dp1(frontier);

        create_decompressed_affected(frontier);



        if(clusters.size() <= 1000)
        {
            std::cout << "\"finalize\" " << mis_set.size() << "nodes: ";
            for(auto& ptr : mis_set)
                std::cout << ptr->cluster_ptr->index << " ";
            std::cout << std::endl;
        }
        parlay::parallel_for(0, mis_set.size(), [&] (T i){
            finalize(mis_set[i]->next);
            auto& check_ptr = mis_set[i]; // TODO remove
        });

        
        std::cout << "checking level " << count+1 << " after decompression and finalizing" << std::endl;
        checker_nodes = parlay::tabulate(clusters.size(), [&] (T i){
            return clusters[i].adjacency[count+1];
        });
        check_consistency(checker_nodes);

        auto update_eligible_set = parlay::filter(frontier, [] (auto node_ptr) {
            return is_update_eligible(node_ptr);
        });

        update_eligible_set = parlay::filter(update_eligible_set, [] (auto node_ptr) {
            return node_ptr->get_num_neighbours_live() <= 2;
        });

        set_MIS(update_eligible_set, true);
        mis_set = parlay::filter(update_eligible_set, [] (auto node_ptr){
            return node_ptr->cluster_ptr->state & IS_MIS_SET;
        });

        if(clusters.size() <= 100)
        {
            printTree(clusters, count+2);
        }
        if(clusters.size() <= 1000)
        {
            std::cout << "contracted " << mis_set.size() << "nodes: ";
            for(auto& ptr : mis_set)
                std::cout << ptr->cluster_ptr->index << " ";
            std::cout << std::endl;
        }

        parlay::parallel_for(0, mis_set.size(), [&] (T i){
            mis_set[i]->cluster_ptr->first_contracted_node = mis_set[i];
            if(mis_set[i]->next == nullptr)
            {
                std::cout << "shouldn't be here" << std::endl;
                exit(1);
            }
            contract(mis_set[i]->next, true);
        });

        std::cout << "checking level " << count+1 << " after contracting" << std::endl;
        checker_nodes = parlay::tabulate(clusters.size(), [&] (T i){
            return clusters[i].adjacency[count+1];
        });
        check_consistency(checker_nodes);

        frontier = new_frontier;

        parlay::parallel_for(0, frontier.size(), [&] (T i){
            if(frontier[i]->next == nullptr) // TODO remove
            {
                std::cout << "Should never happen!" << std::endl;
                exit(1);
            }
            frontier[i] = frontier[i]->next;
        });

        frontier = parlay::filter(frontier, [] (auto node_ptr) {
            return node_ptr->state & affected;
        });
        
        frontier = parlay::filter(frontier, [] (auto node_ptr){
            if(first_condition(node_ptr) || second_condition(node_ptr) || third_condition(node_ptr))
            {
                return true;
            }
            return false;
        });
        parlay::parallel_for(0, frontier.size(), [&] (T i) {
            frontier[i]->state |= affected;
        });

        std::cout << "Frontier size " << frontier.size() << std::endl;        

        count++;
    }while(count < max_count && frontier.size());

    std::cout << "[dynamic] exited" << std::endl;

    if(count == max_count)
    {
        std::cout << red << "[dynamic] Definitely went into infinite loop (didn't terminate for " << max_count << " iters)" << reset << std::endl;
    }

    return;
}

#endif