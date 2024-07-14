#ifndef ADJ_LL_H
#define ADJ_LL_H

#include "cluster.h"
#include <parlay/alloc.h>


const short empty_type = 0;
const short base_vertex = 1;
const short base_edge = 2;
const short unary_cluster = 4;
const short binary_cluster = 8;
const short nullary_cluster = 16;
const short IS_MIS_SET = 32;
const short live = 64;
const short affected = 128;
const short adjacency_changed = 256;
const short all_contracted_affected = 512;
const short debug_state = 1024;
const short contracts_this_round = 2048;
const short C1 = 4096;
const short C2 = 8192;
const short C3 = 8192 * 2;

const int max_neighbours = 3;


template <typename T>
struct node
{
    std::array<T, max_neighbours> adjacents; 
    node* next = nullptr;
    node* prev = nullptr;
    short state = 0; // example value
    unsigned char contraction_level = 0;
    
    node(std::array<T, max_neighbours> arr, short state_val, unsigned char contraction_level_val, node* next_ptr, node* prev_ptr)
        : adjacents(arr), state(state_val), contraction_level(contraction_level_val), next(next_ptr), prev(prev_ptr)
    {
        
    }
    T& operator[](const short index)
    {
        return this->adjacents[index];
    }

    inline short size(void)
    {
        return this->adjacents.size();
    }
};

template <typename T>
class adjacency_list
{
    using node_allocator = parlay::type_allocator<node<T>>;
    private:
        node<T>* head;
        node<T>* tail;
        unsigned char numel;

        

    public:
        adjacency_list(void)
        {
            this->head = nullptr;
            this->tail = nullptr;
            this->numel = 0;
        }
        adjacency_list(std::array<T, max_neighbours> arr)
        {
            this->head = node_allocator::create(arr, live, 0, nullptr, nullptr);
            // new node(arr, live, 0, nullptr, nullptr);/
            this->tail = this->head;
            this->numel = 1;
        }

        node<T>* add_tail(std::array<T, max_neighbours> arr, short state, unsigned char level)
        {
            if(this->numel == 0)
            {
                this->head = node_allocator::create(arr, live, 0, nullptr, nullptr);
                // new node(arr, live, 0, nullptr, nullptr);/
                this->tail = this->head;
                this->numel = 1;
                return this->head;
            }
            auto new_node = node_allocator::create(arr, state, level, nullptr, this->tail);
            this->tail->next = new_node;
            this->tail = new_node;
            this->numel++;
            return this->tail;
        }
        
        node<T>* add_level(std::array<T, max_neighbours> arr, short state, unsigned char contraction_time)
        {
            return this->add_tail(arr, state, contraction_time);
        }

        node<T>* add_level(node<T>* node_ptr)
        {
            return this->add_tail(node_ptr->adjacents, node_ptr->state, node_ptr->contraction_level);
        }

        unsigned char size(void)
        {
            return this->numel;
        }
        

        node<T>* delete_tail()
        {
            if(this->numel == 0)
            {
                return nullptr;
            }
            auto new_tail = this->tail->prev;
            node_allocator::destroy(this->tail);
            this->tail = new_tail;
            this->numel--;
            if(this->numel == 0)
                this->head = nullptr;
            else
                this->tail->next = nullptr;
            return this->tail;
        }

        node<T>* copy_till_level(adjacency_list* adj, const unsigned char& level)
        {
            this->clear_all(); // assume that there is nothing in this list

            node<T>* other_head = adj->get_head();
            while(other_head != nullptr && other_head->contraction_level <= level)
            {
                this->add_level(other_head);
                other_head = other_head->next;
            }
            
            return this->get_tail();
        }

        /*
            Best effort attempt
        */
        node<T>* operator[](const unsigned char level) const
        {
            if(this->numel == 0)
                return nullptr;
            else if (level == 0)
                return this->head;
            auto ret_ptr = this->head;
            auto prev_ptr = ret_ptr;
            while(ret_ptr != nullptr && ret_ptr->contraction_level != level)
            {
                prev_ptr = ret_ptr;
                ret_ptr = ret_ptr->next;
                if(ret_ptr != nullptr && ret_ptr->contraction_level == level)
                    return ret_ptr;
            }
            return prev_ptr;
        }

        node<T>* adopt(adjacency_list* alt_adj)
        {  
            this->clear_all();
            this->head = alt_adj->head;
            this->tail = alt_adj->tail;
            this->numel = alt_adj->size();
            alt_adj->clear_vars();
            return this->tail;
        }

        void clear_vars(void)
        {
            this->numel = 0;
            this->head = nullptr;
            this->tail = nullptr;
        }

        void clear_all(void)
        {
            while(this->delete_tail() != nullptr);
        }    

        node<T>* get_head()
        {
            return this->head;
        }
        
        node<T>* get_tail()
        {
            return this->tail;
        }

        ~adjacency_list(void)
        {
            this->clear_all();
        }




};

#endif