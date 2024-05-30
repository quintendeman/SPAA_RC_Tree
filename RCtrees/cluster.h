#include <atomic>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../examples/helper/graph_utils.h"
#include "../examples/counting_sort.h"

const short empty_type = 0;
const short base_vertex = 1;
const short base_edge = 2;
const short unary_cluster = 4;
const short binary_cluster = 8;
const short nullary_cluster = 16;
const short live = 256;
const short carrying_weight = 512;
const short internal = 8192;



/*
    This represents a cluster in an RC tree.
    All these variables might be excessive but these flags and such are necessary since will be relying on pointer chasing
*/
template <typename T>
struct cluster
{
    public:
        parlay::sequence<std::atomic<cluster<T>*>> neighbours;
        parlay::sequence<std::atomic<cluster<T>*>> children;
        parlay::sequence<std::atomic<T>> counter = parlay::sequence<std::atomic<T>>(1);
        parlay::sequence<std::atomic<T>> height = parlay::sequence<std::atomic<T>>(1);
        cluster<T>* parent = nullptr;
        T index = -1;
        T colour = -1;
        T data = 0;
        bool is_MIS = false;
        short state = 0; // unaffected, affected or update eligible

        unsigned long get_colour(void)
        {
            return reinterpret_cast<T>(this);
        }

        T get_neighbour_count(void)
        {
            T ret_count = 0;
            for(T j = 0; j < this->neighbours.size(); j++)
            {
                if(this->neighbours[j].load() != nullptr)
                {
                   ret_count++;
                }
            }
            return ret_count;
        }

        /**
         * Should be async safe
         * sets "neighbour" to be anywhere in the neighbours array that was previously null
        */
        bool add_neighbour(cluster<T>* neighbour)
        {
            for(T j = 0; j < this->neighbours.size(); j++)
            {
                if(this->neighbours[j].load() == nullptr)
                {
                    cluster<T>* expected = nullptr;
                    bool swapped = this->neighbours[j].compare_exchange_strong(expected, neighbour);
                    if(swapped)
                    {
                        return true;
                    }
                }
            }
            return false;
        }
        /**
         * Not guaranteed to remove!
         * Race conditions can cause this to "skip" a count
         * Put it while(!remove_neighbour(neighbour))
         * to ensure it removes
        */
        bool remove_neighbour(cluster<T>* neighbour)
        {
            for(T j = 0; j < this->neighbours.size(); j++)
            {
                if(this->neighbours[j].load() == neighbour)
                {
                    cluster<T>* expected = neighbour;
                    bool swapped = this->neighbours[j].compare_exchange_strong(expected, nullptr);
                    if(swapped)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        T get_children_count(void)
        {
            T ret_count = 0;
            for(T j = 0; j < this->children.size(); j++)
            {
                if(this->children[j].load() != nullptr)
                {
                   ret_count++;
                }
            }
            return ret_count;
        }

        /**
         * Should be async safe
         * sets "child" to be anywhere in the child array that was previously null
        */
        bool add_child(cluster<T>* child)
        {
            for(T j = 0; j < this->children.size(); j++)
            {
                if(this->children[j].load() == nullptr)
                {
                    cluster<T>* expected = nullptr;
                    bool swapped = this->children[j].compare_exchange_strong(expected, child);
                    if(swapped)
                    {
                        return true;
                    }
                }
            }
            return false;
        }
        /**
         * Not guaranteed to remove!
         * Race conditions can cause this to "skip" a value
        */
        bool remove_child(cluster<T>* child)
        {
            for(T j = 0; j < this->children.size(); j++)
            {
                if(this->children[j].load() == child)
                {
                    cluster<T>* expected = child;
                    bool swapped = this->children[j].compare_exchange_strong(expected, nullptr);
                    if(swapped)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        /**
         * Get first edge in neighbours array that isn't equal to this_side and not null
         * Suffers from race conditions but MIS should keep it safe
        */
        cluster<T>* get_other_side(cluster<T>* this_side)
        {
            cluster<T>* ret_ptr = NULL;
            for(T j = 0; j < this->neighbours.size(); j++)
            {
                if(this->neighbours[j].load() != nullptr && this->neighbours[j].load() != this_side)
                {
                    return this->neighbours[j].load();
                }
            }
            return ret_ptr;
        }
        /**
         * Set all 2 hop neighbour's is_MIS as MIS
        */
        void set_neighbour_mis(bool MIS)
        {
            for(T i = 0; i < this->neighbours.size(); i++)
            {
                if(this->neighbours[i] == nullptr)
                    continue;
                auto edge = this->neighbours[i].load();
                for(T j = 0; j < edge->neighbours.size(); j++)
                {
                    if(edge->neighbours[j].load() != nullptr && edge->neighbours[j].load() != this)
                        edge->neighbours[j].load()->is_MIS = MIS;
                }
            }
        }

        /**
         * if any of 2 hop neighbours i.e. after edge are true, return true
        */
        bool get_neighbour_MIS(void)
        {
            for(T i = 0; i < this->neighbours.size(); i++)
            {
                if(this->neighbours[i] == nullptr)
                    continue;
                auto edge = this->neighbours[i].load();
                for(T j = 0; j < edge->neighbours.size(); j++)
                {
                    if(edge->neighbours[j].load() != nullptr && edge->neighbours[j].load() != this && edge->neighbours[j].load()->is_MIS == true)
                        return true;

                }
            }

            return false;
        }

        void get_two_neighbouring_edges(cluster<T>** left, cluster<T>** right)
        {
            bool left_found = false;

            for(T i = 0; i < this->neighbours.size(); i++)
            {
                if(this->neighbours[i].load() != nullptr)
                {
                    
                    if(left_found == false)
                    {
                        *left = this->neighbours[i].load();
                        left_found = true;
                    }
                    else
                    {
                        *right = this->neighbours[i].load();
                        break;
                    }
                }    
            }
        }

        bool is_neighbour(cluster<T>* neighbour)
        {
            for(T i = 0; i < this->neighbours.size(); i++)
            {
                if(this->neighbours[i].load() == neighbour)
                    return true;
            }
            return false;
        }
        
        bool is_child(cluster<T>* child)
        {
            for(T i = 0; i < this->children.size(); i++)
            {
                if(this->children[i].load() == child)
                    return true;
            }
            return false;
        }

        template <typename lambdafunc>
        T get_children_contribution(lambdafunc func, T def, cluster<T>* exclude = nullptr)
        {
            if(this->get_children_count() == 0)
            {
                return def;
            }
            bool atleast_one = false;
            T accumulated;
            for(uint i = 0; i < this->children.size(); i++)
            {
                auto child_ptr = this->children[i].load();
                if(child_ptr != nullptr && child_ptr != exclude
                    && (child_ptr->state & binary_cluster || child_ptr->state & base_edge))
                {
                    if(atleast_one == false)
                    {
                        accumulated = child_ptr->data;
                        atleast_one = true;
                    }
                    else
                    {
                        accumulated = func(child_ptr->data, accumulated);
                    }
                }
            }

            if(atleast_one == false)
                return def;

            return accumulated;
        }




        void print()
        {
            if(this == NULL)
            {
                std::cout << "ID: NULL" << std::endl;
            }
            std::cout << "ID:" << this->index; 
            if(this->state & binary_cluster)
                std::cout << " binary";
            else if (this->state & unary_cluster)
                std::cout << " unary";
            else if (this->state & nullary_cluster)
                std::cout << " nullary";
            std::cout << std::endl;
            if(this->parent)
                std::cout << "parent: " << parent->index << std::endl;
            std::cout << "neighbours: ";
            for(uint i = 0; i < this->neighbours.size(); i++)
            {
                if(this->neighbours[i].load() != nullptr)
                    std::cout << this->neighbours[i].load()->index << " ";
            }
            std::cout << std::endl;
            std::cout << "children: ";
            for(uint i = 0; i < this->children.size(); i++)
            {
                if(this->children[i].load() != nullptr)
                {
                    std::cout << this->children[i].load()->index;
                    this->children[i].load()->simple_print_neighbours();
                }
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }

        void simple_print_neighbours()
        {
            std::cout << "(";
            for(uint i = 0; i < this->neighbours.size(); i++)
                if(this->neighbours[i].load() != nullptr)
                    std::cout << this->neighbours[i].load()->index << " ";
            std::cout << ") ";
        }

        void print_ancestory()
        {
            if(this == nullptr)
            {
                std::cout << "NULL " << std::endl;
                return;
            }
            auto cluster_ptr = this;

            while(cluster_ptr)
            {
                std::cout << cluster_ptr->index << " ";
                cluster_ptr=cluster_ptr->parent;
            }
            std::cout << std::endl;
        }
};