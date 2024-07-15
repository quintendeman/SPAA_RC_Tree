#ifndef CLUSTER_H
#define CLUSTER_H

// Declarations and definitions



#include <atomic>
#include <algorithm>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../examples/counting_sort.h"
#include "adjacency_linked_list.h"



const char neighbour_type = 1;
const char parent_type = 2;
const char child_type = 4;
const char edge_type = 8;
const char added_type = 16;
const char deleted_type = 32;



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
#include <atomic>
#include <cstring> // for std::memcpy if needed

template <typename T, typename D>
struct cluster
{
public:
    T index;
    adjacency_list adjacency;
    std::atomic<T> counter;
    T tiebreak;
    cluster<T,D>* parent = nullptr;
    D data;
    int state;

    cluster(void) : counter(0)
    {

    }

    cluster(const cluster& other) :
        counter(other.counter.load())
    {

    }
    
    short get_height(void)
    {
        if(this->adjacency.size() > 0)
        {
            return this->adjacency.get_tail()->contraction_level;
        }
        assert(true && "Shouldn't ever get the height of an uncontracted node");
    }

    short get_parent()
    {
        return this->parent;
    }

    unsigned long get_default_colour(void)
    {
        return reinterpret_cast<unsigned long>(this);
    }

    // Filled with null pointers and no state
    void add_empty_level(void)
    {
        std::array<node*, max_neighbours> arr = { nullptr };
        if(this->adjacency.size())
            this->adjacency.add_level(arr, empty_type, this->adjacency.get_tail()->contraction_level, this->index);
        else
            this->adjacency.add_level(arr, empty_type, 0, this->index);
    }

    void add_empty_level(int state, unsigned char level)
    {
        std::array<node*, max_neighbours> arr = { nullptr };
        this->adjacency.add_level(arr, state, level, this->index);
    }

    void add_ptr_to_highest_level(node* nodeptr)
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

    void print(void)
    {
        std::cout << bold << white << this->index << " " << reset;
        for(auto i = 0; i < this->adjacency.size(); i++)
        {
            std::cout << magenta << "[ ";


            const auto& node_ptr_arr = this->adjacency[i]->adjacents;
            for(const auto& ptr : node_ptr_arr)
            {
                if(ptr != nullptr && ptr->index != -1)
                    std::cout << bold << blue << "(" << ptr->index <<") "<< reset << magenta;

                if(ptr != nullptr)
                {
                    const auto& nbr_nodes_list = ptr->adjacents;
                    for(const auto& nbr_node : nbr_nodes_list)
                        if(nbr_node != nullptr && nbr_node->index != this->index)
                            std::cout << nbr_node->index << " ";
                }
                std::cout << " ";
            }
            std::cout << "] ";
        }
        std::cout << reset << std::endl;
    }

};



#endif