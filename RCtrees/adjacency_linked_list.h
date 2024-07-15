#ifndef ADJ_LL_H
#define ADJ_LL_H

#include "cluster.h"
#include <parlay/alloc.h>


const int empty_type = 0;
const int base_vertex = 1;
const int base_edge = 2;
const int unary_cluster = 4;
const int binary_cluster = 8;
const int nullary_cluster = 16;
const int IS_MIS_SET = 32;
const int live = 64;
const int affected = 128;
const int adjacency_changed = 256;
const int all_contracted_affected = 512;
const int debug_state = 1024;
const int contracts_this_round = 2048;

const int max_neighbours = 3;

struct node
{
    std::array<node*, max_neighbours> adjacents; 
    node* next = nullptr;
    node* prev = nullptr;
    long index = -1;
    int state = 0; // example value
    unsigned char contraction_level = 0;
    
    node(std::array<node*, max_neighbours> arr, int state_val, unsigned char contraction_level_val, node* next_ptr, node* prev_ptr)
        : adjacents(arr), state(state_val), contraction_level(contraction_level_val), next(next_ptr), prev(prev_ptr)
    {
        
    }
    node*& operator[](const int index)
    {
        return this->adjacents[index];
    }

    inline int size(void)
    {
        return this->adjacents.size();
    }
};

class adjacency_list
{
    using node_allocator = parlay::type_allocator<node>;
    private:
        node* head;
        node* tail;
        unsigned char numel;

    public:
        adjacency_list(void)
        {
            this->head = nullptr;
            this->tail = nullptr;
            this->numel = 0;
        }
        adjacency_list(std::array<node*, max_neighbours> arr)
        {
            this->head = node_allocator::create(arr, live, 0, nullptr, nullptr);
            // new node(arr, live, 0, nullptr, nullptr);/
            this->tail = this->head;
            this->numel = 1;
        }


        node* add_tail(std::array<node*, max_neighbours> arr, int state, unsigned char level, long index_in)
        {
            if(this->numel == 0)
            {
                this->head = node_allocator::create(arr, live, 0, nullptr, nullptr);
                // new node(arr, live, 0, nullptr, nullptr);/
                this->tail = this->head;
                this->numel = 1;
                this->tail->index = index_in;
                return this->head;
            }
            auto new_node = node_allocator::create(arr, state, level, nullptr, this->tail);
            this->tail->next = new_node;
            this->tail = new_node;
            this->numel++;
            this->tail->index = index_in;
            return this->tail;
        }
        
        node* add_level(std::array<node*, max_neighbours> arr, int state, unsigned char contraction_time, long index_in) // TODO change name from contraction time to level 
        {
            return this->add_tail(arr, state, contraction_time, index_in);
        }

        node* add_level(node* node_ptr, long index_in)
        {
            return this->add_tail(node_ptr->adjacents, node_ptr->state, node_ptr->contraction_level, index_in);
        }

        unsigned char size(void)
        {
            return this->numel;
        }
        

        node* delete_tail()
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

        /*
            Best effort attempt
        */
        node* operator[](const unsigned char level) const
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

        node* adopt(adjacency_list* alt_adj)
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

        node* get_head()
        {
            return this->head;
        }
        
        node* get_tail()
        {
            return this->tail;
        }

        ~adjacency_list(void)
        {
            this->clear_all();
        }
};

#endif