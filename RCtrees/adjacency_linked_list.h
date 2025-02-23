#ifndef ADL_H
#define ADL_H

#include "cluster.h"
template<typename T, typename D>
class cluster;


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
const int dont_finalize = 4096;
const int update_eligible = 8192;
const int one_sided = 8192 * 2;
const int is_marked = 8192 * 4;
const int is_marked_endpoint = 8192 * 8;
const int subtree_query_marked = 8192 * 16;


const int max_neighbours = 3;

template<typename T, typename D>
struct node
{
    std::array<node<T,D>*, max_neighbours> adjacents; 
    node<T,D>* next = nullptr;
    node<T,D>* prev = nullptr;
    cluster<T,D>* cluster_ptr = nullptr;
    int state = 0; // example value
    unsigned char contraction_level = 0;
    
    node(std::array<node<T,D>*, max_neighbours> arr, int state_val, unsigned char contraction_level_val, node<T,D>* next_ptr, node<T,D>* prev_ptr)
        : adjacents(arr), state(state_val), contraction_level(contraction_level_val), next(next_ptr), prev(prev_ptr)
    {
        
    }
    node*& operator[](const int index)
    {
        return this->adjacents[index];
    }

    T& tiebreak(void)
    {
        return this->cluster_ptr->tiebreak;
    }

    bool add_ptr(node<T,D>* new_node_ptr)
    {
        for(auto& ptr : this->adjacents)
            if(ptr == new_node_ptr)
            {
                return true;
            }

        for(auto& ptr : this->adjacents)
            if(ptr == nullptr)
            {
                ptr = new_node_ptr;
                return true;
            }

        return false;
    }

    T index(void)
    {
        return this->cluster_ptr->index;
    }

    inline int size(void)
    {
        return this->adjacents.size();
    }

    // Get number of neighbours assuming this node is live
    // To get boundary vertices, do something else
    short get_num_neighbours_live(void)
    {
        short ret_num = 0;
        for(const auto& adj_ptr : this->adjacents)
            if(adj_ptr != nullptr && (adj_ptr->state & (base_edge | binary_cluster)))
                ret_num++;
        return ret_num;
    }

    ~node(void)
    {
        this->adjacents.fill(nullptr);
        // this->state = 3;
    }

};


template<typename T, typename D>
class adjacency_list
{
    using node_allocator = parlay::type_allocator<node<T,D>>;
    private:
        node<T,D>* head;
        node<T,D>* tail;
        unsigned char numel;

    public:
        adjacency_list(void)
        {
            this->head = nullptr;
            this->tail = nullptr;
            this->numel = 0;
        }
        adjacency_list(std::array<node<T,D>*, max_neighbours> arr)
        {
            this->head = node_allocator::create(arr, live, 0, nullptr, nullptr);
            // new node(arr, live, 0, nullptr, nullptr);/
            this->tail = this->head;
            this->numel = 1;
        }


        node<T,D>* add_tail(std::array<node<T,D>*, max_neighbours> arr, int state, unsigned char level, cluster<T,D>* cluster_ptr_in)
        {
            if(this->numel == 0)
            {
                this->head = node_allocator::create(arr, state, 0, nullptr, nullptr);
                // new node(arr, live, 0, nullptr, nullptr);/
                this->tail = this->head;
                this->numel = 1;
                this->tail->cluster_ptr = cluster_ptr_in;
                return this->head;
            }
            auto new_node = node_allocator::create(arr, state, level, nullptr, this->tail);
            this->tail->next = new_node;
            this->tail = new_node;
            this->numel++;
            this->tail->cluster_ptr = cluster_ptr_in;
            return this->tail;
        }
        
        node<T,D>* add_level(std::array<node<T,D>*, max_neighbours> arr, int state, unsigned char contraction_time, cluster<T,D>* cluster_ptr_in) // TODO change name from contraction time to level 
        {
            return this->add_tail(arr, state, contraction_time, cluster_ptr_in);
        }

        node<T,D>* add_level(node<T,D>* node_ptr, long index_in)
        {
            return this->add_tail(node_ptr->adjacents, node_ptr->state, node_ptr->contraction_level, node_ptr->cluster_ptr_in);
        }

        unsigned char size(void)
        {
            return this->numel;
        }
        

        node<T,D>* delete_tail()
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
        node<T,D>* operator[](const unsigned char level) const
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

        node<T,D>* adopt(adjacency_list* alt_adj)
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

        node<T,D>* get_head() const
        {
            return this->head;
        }
        
        node<T,D>* get_tail()
        {
            return this->tail;
        }

        ~adjacency_list(void)
        {
            this->clear_all();
        }
};

#endif