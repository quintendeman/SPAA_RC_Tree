#ifndef CLUSTER_H
#define CLUSTER_H


#include <atomic>
#include <algorithm>
#include "adjacency_linked_list.h"
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include <atomic>

const char* reset = "\033[0m";
const char* black = "\033[30m";
const char* red = "\033[31m";
const char* green = "\033[32m";
const char* yellow = "\033[33m";
const char* blue = "\033[34m";
const char* magenta = "\033[35m";
const char* cyan = "\033[36m";
const char* white = "\033[37m";
const char* bright_black = "\033[90m";
const char* bright_red = "\033[91m";
const char* bright_green = "\033[92m";
const char* bright_yellow = "\033[93m";
const char* bright_blue = "\033[94m";
const char* bright_magenta = "\033[95m";
const char* bright_cyan = "\033[96m";
const char* bright_white = "\033[97m";
const char* bold = "\033[1m";
const char* dim = "\033[2m";
const char* italic = "\033[3m";
const char* underline = "\033[4m";
const char* blink = "\033[5m";
const char* reverse = "\033[7m";
const char* hidden = "\033[8m";
const char* backspace = "\033[D";


/*
    This represents a cluster in an RC tree.
    All these variables might be excessive but these flags and such are necessary since will be relying on pointer chasing
*/
template<typename T, typename D>
class adjacency_list;

template<typename T, typename D>
class node;

template<typename T, typename D>
std::pair<node<T,D>*,node<T,D>*> get_boundary_vertices_dead(node<T,D>* vertex);

template <typename T, typename D>
struct cluster
{
public:
    T index;
    adjacency_list<T,D> adjacency;
    std::atomic<T> counter;
    std::atomic<T> tiebreak;
    std::atomic<T>& colour = tiebreak; // does this work
    std::atomic<cluster<T,D>*> parent = nullptr; // TODO convert to atomic and modify RCdynamic.h too to prevent obscure scenario of extra work (but should still be n)
    std::array<cluster<T,D>*, max_neighbours> children;
    node<T,D>* first_contracted_node = nullptr;
    D data;

    cluster<T,D>* max_weight_edge = nullptr;
    std::array<std::pair<T,D>, max_neighbours> CGGraphEntry;
    std::array<cluster<T,D>*, max_neighbours> CGGraphEdgePtr;
    std::array<std::atomic<bool>, max_neighbours> CGEntryValid;
    T vertex_count = 0;

    // subtree stuff
    std::array<D, max_neighbours> partial_sums;
    std::array<std::atomic<char>, max_neighbours> partial_sum_complete;


    // short finalize_time = 0; // TODO remove
    // short asked_to_increment = 0; // TODO remove
    int state;

    cluster(void) : counter(0), colour(tiebreak)
    {
        for(short i = 0; i < max_neighbours; i++)
            this->children[i] = nullptr;
    }

    cluster(const cluster& other) :
        counter(other.counter.load())
    {
        for(short i = 0; i < max_neighbours; i++)
            this->children[i] = other.children[i];
    }
    
    short get_height(void)
    {
        if(this->adjacency.size() > 0 && this->first_contracted_node != nullptr)
        {
            return this->first_contracted_node->contraction_level;
        }
        assert(true && "Shouldn't ever get the height of an uncontracted node");
        return -1;
    }

    cluster<T,D>* get_parent()
    {
        return this->parent;
    }

    void insertCGentry(T other_index, D value, cluster<T,D>* corresponding_edge)
    {
        // std::cout << "Inserting " << this->index << " -- " << other_index << " " << value << "/" << corresponding_edge->data << "@" << corresponding_edge->index << std::endl;
        for(unsigned int i = 0; i < max_neighbours; i++)
        {
            bool expected = false;
            bool desired = true;
            bool success = this->CGEntryValid[i].compare_exchange_strong(expected, desired);
            if(success)
            {
                CGGraphEntry[i].first = other_index;
                CGGraphEntry[i].second = value;
                CGGraphEdgePtr[i] = corresponding_edge;
                return;
            }
        }
    }

    void clearCGgraphEntry(T other_index)
    {
        for(unsigned int i = 0; i < max_neighbours; i++)
        {
            if(CGGraphEntry[i].first == other_index)
                CGEntryValid[i] = false;
        }
    }

    void clearCGALL()
    {
        for(unsigned int i = 0; i < max_neighbours; i++)
        {
            CGEntryValid[i] = false;
        }
    }

    short getCGNumEdges()
    {
        short retnum = 0;
        for(unsigned int i = 0; i < max_neighbours; i++)
        {
            if((CGEntryValid[i] == true) && (CGGraphEntry[i].first > this->index))
                retnum++;
        }
        return retnum;
    }

    //default color is the numerical value of the pointer
    unsigned long get_default_colour(void)
    {
        return reinterpret_cast<unsigned long>(this);
    }

    // Filled with null pointers and no state
    void add_empty_level(void)
    {
        std::array<node<T,D>*, max_neighbours> arr = { nullptr };
        if(this->adjacency.size())
            this->adjacency.add_level(arr, empty_type, this->adjacency.get_tail()->contraction_level, this);
        else
            this->adjacency.add_level(arr, empty_type, 0, this);
    }

    void add_empty_level(int state, unsigned char level)
    {
        std::array<node<T,D>*, max_neighbours> arr = { nullptr };
        this->adjacency.add_level(arr, state, level, this);
    }

    void add_ptr_to_highest_level(node<T,D>* nodeptr)
    {
        auto& w_ptr_arr = this->adjacency.get_tail()->adjacents;
        for(auto& w_ptr : w_ptr_arr)
            if(w_ptr == nullptr)
            {
                w_ptr = nodeptr;
                break;
            }
    }

    void replicate_last_level(int state, unsigned char level)
    {
        if(this->adjacency.size() == 0)
            this->add_empty_level(state, level);
        else
        {
            auto& w_ptr_arr = this->adjacency.get_tail()->adjacents;
            this->adjacency.add_level(w_ptr_arr, state, level, this->index);
        }
    }
    void replicate_last_level(void)
    {
        this->replicate_last_level(empty_type, 1);
    }

    void set_neighbour_mis(bool IS_MIS, node<T,D>* node_ptr = nullptr)
    {
        node<T,D>* this_node = node_ptr;
        if(node_ptr == nullptr)
            this_node = this->adjacency.get_tail();
        for(auto& edge_ptr : this_node->adjacents)
        {
            if(edge_ptr != nullptr && edge_ptr->state & (base_edge | binary_cluster))
            {
                auto other_ptr = get_other_side(this_node, edge_ptr);
                if(IS_MIS == true)
                    other_ptr->cluster_ptr->state |= IS_MIS_SET;
                else
                    other_ptr->cluster_ptr->state &= (~IS_MIS_SET);
            }
        }
    }

    bool get_neighbour_mis(node<T,D>* node_ptr = nullptr)
    {
        node<T,D>* this_node = node_ptr;
        if(node_ptr == nullptr)
            this_node = this->adjacency.get_tail();
        for(auto& edge_ptr : this_node->adjacents)
        {
            if(edge_ptr != nullptr && edge_ptr->state & (base_edge | binary_cluster))
            {
                auto other_ptr = get_other_side(this_node, edge_ptr);
                if(other_ptr->cluster_ptr->state & IS_MIS_SET)
                    return true;
            }
        }
        return false;
    }

    void find_endpoints(cluster<T,D>*& left, cluster<T,D>*& right)
    {
        left = right = nullptr;
        if(this->adjacency.get_head()->state & base_edge)
        {
            bool first = true;
            for(auto& adjacent : this->adjacency.get_head()->adjacents)
            {
                if(adjacent == nullptr)
                    continue;
                if(first)
                {
                    first = false;
                    left = adjacent->cluster_ptr;
                }
                else
                    right = adjacent->cluster_ptr;
            }
        }
        else
        {
            std::cout << "NOT A BASE EDGE!" << std::endl;
            exit(1);
        }
    }

    void find_boundary_vertices(cluster<T,D>*& left, D& lval, cluster<T,D>*& right, D& rval, D defretval, short level = -1) //TODO convert to return a struct, this is irritating
    {
        left = right = nullptr;
        lval = rval = defretval;
        node<T,D>* node_ptr = nullptr;
        if(level == -1)
            node_ptr = this->first_contracted_node;
        else
            node_ptr = this->adjacency[level];
        
        if(node_ptr->next->state & nullary_cluster)
            return;
        else if (node_ptr->next->state & unary_cluster)
        {
            // find the one edge going to the other side
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & (base_edge | binary_cluster))
                {
                    left = get_other_side(node_ptr, ptr)->cluster_ptr;
                    lval = ptr->cluster_ptr->data;
                    right = left;
                    rval = lval;
                    return;
                }
            }
            std::cout << "SHOULD NEVER REACH HERE" << std::endl;
            exit(1);
        }
        else if (node_ptr->next->state & binary_cluster)
        {
            if(level == -1)
                node_ptr = this->first_contracted_node;
            node<T,D>* left_edge = nullptr;
            // find the one edge going to the other side
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & (base_edge | binary_cluster))
                {
                    left = get_other_side(node_ptr, ptr)->cluster_ptr;
                    lval = ptr->cluster_ptr->data;
                    left_edge = ptr;
                    break;
                }
            }
            // find the one edge going to the other side
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & (base_edge | binary_cluster) && ptr != left_edge)
                {
                    right = get_other_side(node_ptr, ptr)->cluster_ptr;
                    rval = ptr->cluster_ptr->data;
                    return;
                }
            }
        }
        else
        {
            std::cout << "IT SHOULD NEVER REACH HERE" << std::endl;
            exit(1);
        }
    }
    
    void print(unsigned char level = -1)
    {
        if(this->state & is_marked_endpoint)
            std::cout << italic << red << this->index << " " << reset;
        else if(this->state & is_marked)
            std::cout << bright_red << this->index << " " << reset;
        else
            std::cout << this->index << " " << reset;
        std::cout << bright_white <<  this->get_height() << " " << reset;
        std::cout << bright_green << this->data << " " << reset;
        std::cout << bright_magenta << this->counter << " " << reset;
        std::cout << bright_yellow;
        if(this->parent.load())
            std::cout << this->parent.load()->index << " ";
        else
            std::cout << "nl ";
        std::cout << reset;

        if(this->max_weight_edge)
            std::cout << max_weight_edge->index << "/" << max_weight_edge->data;
        else
            std::cout << "nl/-";
        std::cout << " ";

        std::cout << red;
        for(short i = 0; i < max_neighbours; i++)
            std::cout << this->partial_sums[i] << "/" << (int) this->partial_sum_complete[i] << " ";
        

        std::cout << bold << blue << "[ ";
        for(auto& child : this->children)
        {
            if(child != nullptr)
                std::cout << child->index << "/" << child->data << " ";
        }
        std::cout << "] " << reset;

        // for(auto i = 0; i < (this->adjacency.get_head()->state & affected ? 2 :  this->adjacency.size()); i++)
        for(volatile auto i = 0; i < this->adjacency.size(); i++)
        {
            if(level != -1 && i == level)
                break;
            const auto& node_ptr_arr = this->adjacency[i]->adjacents;

            std::cout << bold << bright_yellow << i << reset;    

            volatile auto node_ptr_test = this->adjacency[i];
            if(node_ptr_test->state & nullary_cluster)
                std::cout << red << "N";
            if(node_ptr_test->state & unary_cluster)
                std::cout << yellow << "U";
            if(node_ptr_test->state & binary_cluster)
                std::cout << green << "B";
            if(node_ptr_test->state & update_eligible)
                std::cout << green << "E";
            
            if(this->adjacency[i] == this->first_contracted_node)
                std::cout << bright_blue << "C";
            if(this->adjacency[i]->state & affected)
                std::cout << bright_blue << "A";
            // if(is_update_eligible(this->adjacency[i]))
            //     std::cout << bright_blue << "E";
            
                std::cout << magenta << "[ ";
            
                for(const auto& ptr : node_ptr_arr)
                {
                    if(ptr != nullptr && ptr->cluster_ptr->index != -1)
                        std::cout << bold << blue << "(" << ptr->cluster_ptr->index <<") "<< reset << magenta;

                    // && ptr->state & (binary_cluster | base_edge)
                    if(ptr != nullptr )
                    {
                        const auto& nbr_nodes_list = ptr->adjacents;
                        for(const auto& nbr_node : nbr_nodes_list)
                            if(nbr_node != nullptr && nbr_node->cluster_ptr->index != this->index)
                                std::cout << nbr_node->cluster_ptr->index << " ";
                            // else if(nbr_node != nullptr)
                            //     std::cout << "x" << nbr_node->cluster_ptr->index << "";
                            
                    }
                    // else
                    //     std::cout << "nl ";
                    std::cout << " ";
                }
                
                std::cout << white << this->adjacency[i]->get_num_neighbours_live() << " " << magenta;
            std::cout << "]";
        }

        // std::cout << reset << bright_white << "F=" << this->finalize_time << " I=" << this->asked_to_increment;
        
        std::cout << reset << std::endl;
    }

};

template<typename T, typename D>
node<T,D>* get_other_side(node<T,D>* myself, node<T,D>* edge_ptr)
{

    if(edge_ptr != nullptr && edge_ptr->state & (base_edge | binary_cluster))
    {
        for(auto& other_ptr : edge_ptr->adjacents)
            if(other_ptr != nullptr && other_ptr != myself)
                return other_ptr;
        return nullptr;
    }
    else
        return nullptr;
}

template<typename T, typename D>
short get_parent_child_index(cluster<T,D>* cluster_ptr)
{
    if(cluster_ptr == nullptr || cluster_ptr->parent == nullptr)
        return -1;
    for(short k = 0; k < max_neighbours; k++)
    {
        auto& child = cluster_ptr->parent->children[k];
        if(child == cluster_ptr)
            return k;
    }
    return -1;
}


#endif