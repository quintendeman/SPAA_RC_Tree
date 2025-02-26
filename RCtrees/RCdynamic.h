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
#include "RC_test.h"

static const char PRINT_DYNAMIC = 0;




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
        parlay::sequence<node<T,D>*> ret_seq = parlay::sequence<node<T,D>*>(max_neighbours + 1, nullptr);
        for(short i = 0; i < max_neighbours; i++)
        {
            if(node_ptr->adjacents[i] == nullptr || ((node_ptr->adjacents[i]->state & (binary_cluster | base_edge)) == 0))
                continue;
            ret_seq[i] = get_other_side(node_ptr, node_ptr->adjacents[i]);
        }
        ret_seq[max_neighbours] = node_ptr;
        return ret_seq;
    }));

    frontier = parlay::flatten(parlay::tabulate(frontier.size(), [&] (T I) {
        auto& node_ptr = frontier[I];
        parlay::sequence<node<T,D>*> ret_seq = parlay::sequence<node<T,D>*>(max_neighbours + 1, nullptr);
        if(node_ptr != nullptr)
            for(short i = 0; i < max_neighbours; i++)
            {
                if(node_ptr->adjacents[i] == nullptr || ((node_ptr->adjacents[i]->state & (binary_cluster | base_edge)) == 0))
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
    parlay::parallel_for(0, affected_nodes.size(), [&] (T i){
        auto& aff_node = affected_nodes[i];
        if(aff_node->next == nullptr)
            aff_node->cluster_ptr->add_empty_level((aff_node->state & (~(binary_cluster | unary_cluster | nullary_cluster))) | live, aff_node->contraction_level+1);
        
        for(short e = 0; e < aff_node->adjacents.size(); e++)
        {
            if(aff_node->adjacents[e] == nullptr)
                aff_node->next->adjacents[e] = nullptr;
            else if(aff_node->adjacents[e]->state & unary_cluster)
                aff_node->next->adjacents[e] = nullptr;    
            else
            {
                if(aff_node->adjacents[e]->next != aff_node->next->adjacents[e])
                {
                    auto old_other_side = get_other_side(aff_node,aff_node->adjacents[e]);
                    if(old_other_side != nullptr && old_other_side->next != aff_node->next->adjacents[e])
                    {
                        if (affected_nodes.size() <= 100 && PRINT_DYNAMIC) {
                            std::cout << "Deleted " 
                                    << aff_node->index() << " " 
                                    << aff_node->adjacents[e]->index() << " " 
                                    << old_other_side->index() << " " 
                            
                                    << std::endl;
                        }
                        aff_node->next->adjacents[e] = nullptr;
                        
                    }
                    else
                        continue;
                }
            }
        }
    });

    // Did I originally rake into something?
    parlay::parallel_for(0, affected_nodes.size(), [&] (T I){
        auto& aff_node = affected_nodes[I];
        if(!is_update_eligible(aff_node))
            return;
        if((aff_node->next->state & unary_cluster) == 0)
            return;
        if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
            std::cout << "rake should be reversed: " << aff_node->index() << std::endl;

        for(short i = 0; i < aff_node->adjacents.size(); i++)
        {
            const auto old_edge = aff_node->adjacents[i];
            const auto new_edge = aff_node->next->adjacents[i];
            if(old_edge == nullptr || new_edge == nullptr)
                continue;

            // did I rake into someone?
            if(old_edge->state & (binary_cluster | base_edge))
            {
                auto old_neighbour = get_other_side(aff_node, old_edge);
                auto new_neighbour = old_neighbour->next;
                
                for(short o = 0; o < old_neighbour->adjacents.size(); o++)
                {
                    if(old_neighbour->adjacents[o] == old_edge && new_neighbour->adjacents[o] == aff_node->next)
                    {
                        if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                            std::cout << "rake got reversed: " << aff_node->index() << std::endl;
                        
                        new_neighbour->adjacents[o] = new_edge;
                        new_neighbour->state |= affected | live;
                    }
                }
            }
            else
                continue;

        }
    });

    // Did I compress into something?
    parlay::parallel_for(0, affected_nodes.size(), [&] (T I){
        auto& aff_node = affected_nodes[I];
        if(!is_update_eligible(aff_node))
            return;
        
        for(short i = 0; i < aff_node->adjacents.size(); i++)
        {
            const auto old_edge = aff_node->adjacents[i];
            const auto new_edge = aff_node->next->adjacents[i];            
            if(old_edge == nullptr)
                continue;
            if(new_edge == nullptr)
                continue;

            if(old_edge->state & (binary_cluster | base_edge))
            {
                if(aff_node->next->state & binary_cluster)
                {
                    const auto old_neighbour = get_other_side(aff_node, old_edge);
                    short k;
                    for(k = 0; k < old_neighbour->adjacents.size(); k++)
                    {
                        if(old_neighbour->adjacents[k] == old_edge)
                            break;
                    }
                    old_neighbour->next->adjacents[k] = old_edge->next;
                    aff_node->next->adjacents[i] = old_edge->next;
                    if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                        std::cout << "compress reversed: " << aff_node->index() << " " << old_edge->index() << "->" << old_neighbour->index() << std::endl;
                    old_neighbour->next->state |= adjacency_changed | affected;
                }
            }
        }
    });

    // is an edge missing?

    parlay::parallel_for(0, affected_nodes.size(), [&] (T I){
        auto& aff_node = affected_nodes[I];
            for(short i = 0; i < aff_node->adjacents.size(); i++)
            {
                const auto old_edge = aff_node->adjacents[i];
                const auto new_edge = aff_node->next->adjacents[i];           
                 
                if(old_edge == nullptr)
                    continue;
                if((old_edge->state & (binary_cluster | base_edge)) && (new_edge == nullptr ))
                {
                    auto old_neighbour = get_other_side(aff_node, old_edge);

                    if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                        std::cout << "eligible for binary? " << aff_node->index() << " " << old_edge->index() << " " << old_neighbour->index() << std::endl;
                    old_edge->state |= one_sided;
                }
            }
        });

    parlay::parallel_for(0, affected_nodes.size(), [&] (T I){
        auto& aff_node = affected_nodes[I];
        for(short i = 0; i < aff_node->adjacents.size(); i++)
        {
            const auto old_edge = aff_node->adjacents[i];
            const auto new_edge = aff_node->next->adjacents[i];            
            if(old_edge == nullptr)
                continue;
            if(old_edge->state & one_sided)
            {
                auto old_neighbour = get_other_side(aff_node, old_edge);
                auto v = aff_node->index();
                auto w = old_neighbour->index();

                if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                    std::cout << "eligiblerrr for binary? " << aff_node->index() << " " << old_edge->index() << " " << old_neighbour->index() << std::endl;
                

                if(w < v)
                    continue;
                if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                    std::cout << "Binary added one sided added: " << aff_node->index() << " " << old_edge->index() << " " << old_neighbour->index() << std::endl;
                if(old_edge->next == nullptr)
                {
                    if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                        std::cout << "Binary new constructed " << std::endl;
                    if(old_edge->state & binary_cluster)
                        old_edge->cluster_ptr->add_empty_level(binary_cluster, old_edge->contraction_level + 1);
                    else
                        old_edge->cluster_ptr->add_empty_level(base_edge, old_edge->contraction_level + 1);
                }
                auto newly_created_edge = old_edge->next;
                newly_created_edge->state = old_edge->state;
                newly_created_edge->adjacents.fill(nullptr);
                // newly_created_edge->add_ptr(aff_node->next);
                for(short k = 0; k < old_edge->adjacents.size(); k++)
                {
                    if(old_edge->adjacents[k] == aff_node)
                        old_edge->next->adjacents[k] = aff_node->next;
                }
                aff_node->next->adjacents[i] = newly_created_edge;
                old_edge->state |= one_sided;
            }
        }
    });

    // is an edge missing? -- other side
    parlay::parallel_for(0, affected_nodes.size(), [&] (T I){
        auto& aff_node = affected_nodes[I];
        for(short i = 0; i < aff_node->adjacents.size(); i++)
        {
            const auto old_edge = aff_node->adjacents[i];
            const auto new_edge = aff_node->next->adjacents[i];            
            if(old_edge == nullptr)
                continue;
            // ((old_edge->state & (binary_cluster | base_edge)) && (new_edge == nullptr || new_edge != old_edge->next)) || 
            if((old_edge->state & one_sided)) //TODO or it is changed?
            {
                auto old_neighbour = get_other_side(aff_node, old_edge);
                auto v = aff_node->index();
                auto w = old_neighbour->index();
                if(v < w)
                    continue;
                if(affected_nodes.size() <= 100 && PRINT_DYNAMIC)
                    std::cout << "Binary other sided added: " << aff_node->index() << " " << old_edge->index() << " " << old_neighbour->index() <<  std::endl;
                auto newly_created_edge = old_edge->next;
                // newly_created_edge->add_ptr(aff_node->next);
                for(short k = 0; k < old_edge->adjacents.size(); k++)
                {
                    if(old_edge->adjacents[k] == aff_node)
                        old_edge->next->adjacents[k] = aff_node->next;
                }
                aff_node->next->adjacents[i] = newly_created_edge;
                old_edge->state &= ~one_sided;
            }
        }
        aff_node->next->state = live | affected;
        
    });
    return;
}


template<typename T, typename D>
bool unParent(cluster<T,D>* cluster_ptr)
{
    // if(cluster_ptr->parent == nullptr)
    //     return false;
    cluster<T,D>* expected = cluster_ptr->parent.load();
    cluster<T,D>* old_parent_value = expected;
    if(expected == nullptr)
        return false;
    cluster<T,D>* new_value = nullptr;
    bool success = cluster_ptr->parent.compare_exchange_strong(expected, new_value);
    if(!success)
        return false;

    for(auto& parents_child : expected->children)
        if(parents_child == cluster_ptr)
        {
            parents_child = nullptr;
            return true;
            
        }
    return false;
}

template<typename T, typename D>
void recursive_unParent(cluster<T,D>* cluster_ptr)
{
    auto old_cluster_ptr = cluster_ptr;
    while(cluster_ptr != nullptr)
    {
        old_cluster_ptr = cluster_ptr;
        cluster_ptr = cluster_ptr->parent;
        bool retval = unParent(old_cluster_ptr);
        if(retval == false)
            break;
        // old_cluster_ptr->parent = nullptr;
    }
}

template<typename T, typename D>
void remove_affection_back(node<T,D>* node_ptr)
{
    do
    {
        node_ptr->state &= (~affected);
        node_ptr = node_ptr->prev;
    }while(node_ptr != nullptr && (node_ptr->state & affected));
}

template<typename T, typename D, typename lambdafunc>
void batchInsertEdge( const parlay::sequence<std::pair<T, T>>& delete_edges, const parlay::sequence<std::tuple<T, T, D>>& add_edges, parlay::sequence<cluster<T, D>>& clusters, D defretval, lambdafunc func, bool randomized = false, const bool do_path_query = true)
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
        // if(ret_seq[0] == nullptr || ret_seq[1] == nullptr) // TODO remove
        // {
        //     std::cout << red << "This should never happen, nullptr on head?" << reset << std::endl;
        //     exit(1);
        // }
        ret_seq[0]->cluster_ptr->tiebreak = i*2;
        ret_seq[1]->cluster_ptr->tiebreak = i*2 +1;
        return ret_seq;
    }));

    if(clusters.size() <= 100 && PRINT_DYNAMIC) // TODO remove
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
        if(tree_nodes[i]->cluster_ptr->tiebreak != i)
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
        // unParent(cluster_edge_ptr);
        recursive_unParent(cluster_edge_ptr);
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
                    auto newEdgeClstrPtr = cluster_allocator::create();
                    newEdgeClstrPtr->add_empty_level(base_edge, 0);
                    newEdgeClstrPtr->data = std::get<2>(edge.E);
                    newEdgeClstrPtr->max_weight_edge = newEdgeClstrPtr;
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

    // std::cout << "endpoints size is " << tree_nodes.size() << std::endl;
    // std::cout << "frontier size is " << frontier.size() << std::endl;
    
    frontier = parlay::filter(frontier, [] (auto node_ptr){
        if(first_condition(node_ptr) || second_condition(node_ptr) || third_condition(node_ptr))
        {
            // unParent(node_ptr->cluster_ptr);
            recursive_unParent(node_ptr->cluster_ptr);
            node_ptr->state |= affected;
            return true;
        }
        return false;
    });
    
    // parlay::parallel_for(0, frontier.size(), [&] (T i) {
    //         frontier[i]->state |= affected;
    // });

    // std::cout << "Which gets reduced to " << frontier.size() << std::endl;
    
    // create_decompressed_affected(frontier);
    parlay::sequence<node<T,D>*> mis_set;

    unsigned short count = 0;
    static const short max_count = 200;
    do
    {
        if(PRINT_DYNAMIC || false) // TODO yknow
        {       
            // std::cout << "checking level at the start of round " << count << std::endl;
           auto checker_nodes = parlay::tabulate(clusters.size(), [&] (T i){
                return clusters[i].adjacency[count];
            });
            check_consistency(checker_nodes);
        }

        if(clusters.size() <= 100 && PRINT_DYNAMIC)
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

        if(clusters.size() <= 100 && PRINT_DYNAMIC)
        {
            std::cout << "\"finalize\" " << mis_set.size() << "nodes: ";
            for(auto& ptr : mis_set)
                std::cout << ptr->cluster_ptr->index << " ";
            std::cout << std::endl;
        }
        

        if(PRINT_DYNAMIC || false) // TODO yknow
        {
            // std::cout << "checking level " << count+1 << " after decompression " << std::endl;
            auto checker_nodes = parlay::tabulate(clusters.size(), [&] (T i){
                return clusters[i].adjacency[count+1];
            });
            check_consistency(checker_nodes);
        }

        auto update_eligible_set = parlay::filter(frontier, [] (auto node_ptr) {
            if(is_update_eligible(node_ptr))
            {
                node_ptr->state |= update_eligible;
                return true;
            }
            return false;
        });

        update_eligible_set = parlay::filter(update_eligible_set, [] (auto node_ptr) {
            return node_ptr->get_num_neighbours_live() <= 2;
        });

        set_MIS(update_eligible_set, true, randomized);
        mis_set = parlay::filter(update_eligible_set, [] (auto node_ptr){
            return node_ptr->cluster_ptr->state & IS_MIS_SET;
        });

        if(clusters.size() <= 100 && PRINT_DYNAMIC)
        {
            printTree(clusters, count+2);
        }
        if(clusters.size() <= 100 && PRINT_DYNAMIC)
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
            remove_affection_back(mis_set[i]->next);
            accumulate(mis_set[i], defretval, func, do_path_query);
        });

        if(PRINT_DYNAMIC || false) // TODO yknow
        {
            // std::cout << "checking level " << count+1 << " after contracting" << std::endl;
            auto checker_nodes = parlay::tabulate(clusters.size(), [&] (T i){
                return clusters[i].adjacency[count+1];
            });
            check_consistency(checker_nodes);
        }

        frontier = new_frontier;

        parlay::parallel_for(0, frontier.size(), [&] (T i){
            // if(frontier[i]->next == nullptr) // TODO remove
            // {
            //     std::cout << "Should never happen!" << std::endl;
            //     std::cout << "frontier[i] = " << frontier[i]->index() << std::endl;
            //     exit(1);
            // }
            frontier[i] = frontier[i]->next;
        });

        frontier = parlay::filter(frontier, [] (auto node_ptr) {
            return node_ptr->state & affected;
        });

        
        
        frontier = parlay::filter(frontier, [] (auto node_ptr){
            if(first_condition(node_ptr) || second_condition(node_ptr) || third_condition(node_ptr))
            {
                // unParent(node_ptr->cluster_ptr);
                recursive_unParent(node_ptr->cluster_ptr);
                node_ptr->state |= affected;
                return true;
            }
            return false;
        });

        // TODO merge with above
        // parlay::parallel_for(0, frontier.size(), [&] (T i) {
        //     frontier[i]->state |= affected;
        // });
        if(PRINT_DYNAMIC)
            std::cout << "Finalizing contracted nodes " << std::endl;
        parlay::parallel_for(0, mis_set.size(), [&] (T i){
            finalize(mis_set[i]->next);
        });

        

        // std::cout << "Dynamic Frontier ["  << count << "]  size " << frontier.size() << std::endl;        

        count++;
    }while(count < max_count && frontier.size());

    // std::cout << "[dynamic] exited" << std::endl;

    // check_parents_children(clusters); // todo remove

    if(count == max_count)
    {
        std::cout << red << "[dynamic] Definitely went into infinite loop (didn't terminate for " << max_count << " iters)" << reset << std::endl;
    }

    return;
}

#endif