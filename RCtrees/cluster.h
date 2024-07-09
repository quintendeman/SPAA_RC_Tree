#ifndef CLUSTER_H
#define CLUSTER_H

// Declarations and definitions



#include <atomic>
#include <algorithm>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../examples/counting_sort.h"

const short empty_type = 0;
const short base_vertex = 1;
const short base_edge = 2;
const short unary_cluster = 4;
const short binary_cluster = 8;
const short nullary_cluster = 16;
const short IS_MIS_SET = 32;
const short live = 256;
const short carrying_weight = 512;
const short internal = 1024;
const short affected = 2048;
const short adjacency_changed = 4096;


const char neighbour_type = 1;
const char parent_type = 2;
const char child_type = 4;
const char edge_type = 8;
const char added_type = 16;
const char deleted_type = 32;

const int max_neighbours = 3;


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
    static const short size = max_neighbours*2;
    std::array<cluster<T, D>*, size> ptrs;
    // cluster<T, D>* ptrs[max_neighbours * 2];
    parlay::sequence<T> adjacencies = parlay::sequence<T>(max_neighbours + 1); // representing a graph -- dynamic! [level node node node level node node node]
    parlay::sequence<T> alternate_adjacencies;
    T index = -1;
    T colour = -1; 
    std::atomic<T> counter; // a node will not have more than 255 edges
    D data = -1.0;
    short state = 0; // unaffected, affected or update eligible
    std::array<char, size> types;
    // char types[max_neighbours * 2];
    unsigned char contraction_time = 0;
    
    // Default constructor
    cluster() :  counter(0) {
        for(uint i = 0; i < this->size; i++)
        {
            this->ptrs[i] = nullptr;
            this->types[i] = 0;
        }
    }

    // Copy constructor only for initialization when creating a sequence
    cluster(const cluster& other)
        : index(other.index),
          colour(other.colour),
          data(other.data),
          counter(other.counter.load()), // Load the atomic value
          state(other.state),
          contraction_time(other.contraction_time)
    {
        for(uint i = 0; i < this->size; i++)
        {
            this->ptrs[i] = nullptr;
            this->types[i] = 0;
            this->adjacencies = other.adjacencies;
        }
        return;
    }

    // Method to get colour based off memory address
    unsigned long get_default_colour(void) const
    {
        return reinterpret_cast<unsigned long>(this);
    }

    // To be used only for connecting the base_edges to base_nodes
    void add_initial_neighbours(cluster<T, D>* V,  cluster<T, D>* W)
    {
        this->ptrs[0] = V;
        this->ptrs[1] = W;
        this->types[0] = neighbour_type;
        this->types[1] = neighbour_type;
    }


    short isPtr(T ngbr_index, short defretval = -1) const
    {
        for(short i = 0; i < this->size; i+=2)
            if(this->ptrs[i] != nullptr && this->ptrs[i]->index  == ngbr_index)
                return i;
        return defretval;
    }
    short isPtr(const cluster<T, D>* ngbr_ptr, short defretval = -1) const
    {
        if(ngbr_ptr == nullptr)
            return defretval;
        return this->isPtr(ngbr_ptr->index, defretval);
    }


    // index of neighbour if successful, -1 if not
    // Note, if that pointer already exists, set its type to neighbour
    // and returns that index
    // index i+1 will contain edge for easier management
    short add_neighbour(cluster<T, D>* neighbour, cluster<T, D>* edge)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] == nullptr || this->ptrs[i] == neighbour)
            {
                this->ptrs[i] = neighbour;
                this->ptrs[i+1] = edge;
                this->types[i] |= neighbour_type;
                this->types[i+1] |= edge_type;
                return i;
            }
        }
        return -1;
    }

    /**
     * Removes a node as a neighbour and returns the index
    */
    short remove_neighbour(cluster<T, D>* neighbour)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == neighbour)
            {
                this->types[i]&=(~neighbour_type);
                return i;
            }
        }
        return -1;
    }

    short change_to_neighbour(cluster<T, D>* neighbour)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == neighbour)
            {
                this->types[i]&=(~child_type);
                this->types[i]&=(~edge_type);
                this->types[i]|=neighbour_type;
                return i;
            }
        }
        return -1;
    }

    // changes a ptr from neighbour_type or edge_type to child_type
    short change_to_child(cluster<T, D>* child)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == child)
            {
                this->types[i]&=(~neighbour_type);
                this->types[i]&=(~edge_type);
                this->types[i]|=child_type;
                return i;
            }
        }
        return -1;
    }

    /**
     * designates a neighbour as a parent
    */
    short set_parent(cluster<T, D>* node)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == node)
            {
                this->types[i]=parent_type;
                this->types[i]&=(~neighbour_type);
                return i;
            }
        }
        return -1;
    }

    cluster<T, D>* get_parent(void) const
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i]&parent_type)
            {
                return this->ptrs[i];
            }
        }
        return nullptr;
    }

    cluster<T, D>* get_parent(short& parent_index) const
    {
        parent_index = -1;
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i] & parent_type)
            {
                parent_index = i;
                return this->ptrs[i];
            }
        }
        return nullptr;
    }

    // set the neighbours MIS to val
    // useful when colouring
    void set_neighbour_mis(bool val)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i])
                if(val)
                    this->ptrs[i]->state |= IS_MIS_SET;
                else
                    this->ptrs[i]->state &= (~IS_MIS_SET);
        }
    }

    void set_neighbour_colour(T input_colour)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->types[i] & neighbour_type)
                this->ptrs[i]->colour = input_colour;
        }
    }

    bool is_max_neighbour_colour(void)
    {
        T max_colour = this->colour;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->ptrs[i]->colour > max_colour && this->types[i] & neighbour_type)
                max_colour = this->ptrs[i]->colour;
        }
        if(max_colour == this->colour)
            return true;
        return false;
    }

    cluster<T, D>* get_one_neighbour(void) const
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i] & neighbour_type)
                return this->ptrs[i];
        }
        return nullptr;
    }

    /*
        return true if any neighbours are in MIS
        else return false
    */
    bool get_neighbour_MIS(void) const
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->ptrs[i]->state & IS_MIS_SET)
                return true;
        }
        return false;
    }

    short get_neighbour_count(void) const
    {
        short num_neighbours = 0;
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i] & neighbour_type)
                num_neighbours++;
        }

        return num_neighbours;
    }

    short get_children_count(void) const
    {
        short num_children = 0;
        for(const auto& t : this->types)
            if (t & child_type)
                num_children++;
        return num_children;
    }

    void get_two_neighbours_edges(cluster<T, D>*& neighbour1, cluster<T, D>*& edge1, cluster<T, D>*& neighbour2, cluster<T, D>*& edge2) const
    {
        bool first_neighbour_found = false;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->types[i]&neighbour_type)
            {
                if(first_neighbour_found == false)
                {
                    first_neighbour_found = true;
                    neighbour1 = this->ptrs[i];
                    edge1 = this->ptrs[i+1];
                }
                else
                {
                    neighbour2 = this->ptrs[i];
                    edge2 = this->ptrs[i+1];
                    return;
                }
            }
        }
    }

    // Find boundary vertices by observing own state and ajdacent edges
    void find_boundary_vertices(cluster<T, D>* &l, D& lval, cluster<T, D>* &r, D& rval, const D& defretval) const
    {
        lval = rval = defretval;
        if(this->state & nullary_cluster)
        {
            l = r = nullptr;
            return;
        }
        if(this->state & unary_cluster)
        {
            short parent_index = 0;
            l = r = this->get_parent(parent_index);
            /* l CAN NEVER BE NULLPTR*/
            lval = rval = this->ptrs[parent_index+1]->data;
            return;
        }
        if(this->state & binary_cluster)
        {
            short parent_index;
            r = this->get_parent(parent_index);
            rval = this->ptrs[parent_index+1]->data;
            for(uint i = 0; i < this->size; i+=2)
            {
                if(this->types[i] & neighbour_type && this->ptrs[i] != r)
                {
                    l = this->ptrs[i];
                    lval = this->ptrs[i+1]->data;
                    return;
                }
            }

        }
        return;
    }

    /*
        Mainly used when doing compress

    */
    short overwrite_neighbour(cluster<T, D>* old_neighbour, cluster<T, D>* new_neighbour, cluster<T, D>* new_edge)
    {

        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == old_neighbour)
            {
                this->ptrs[i] = new_neighbour;
                this->types[i]|=neighbour_type;

                this->ptrs[i+1] = new_edge;
                this->types[i+1]|= edge_type;
                this->types[i+1]&= (~neighbour_type);

                return i;
            }
        }

        return -1;
    }

    inline unsigned char get_height() const
    {
        return this->contraction_time;
    }

    void set_height()
    {
        this->height = -1;
        for(uint i = 0; i < this->size; i++)
        {
            if(this->types[i] & child_type && this->ptrs[i] != nullptr)
                this->height = std::max(this->height, this->ptrs[i]->height);
        }
        this->height = this->height + 1;
        return;
    }

    void print_as_edge(void)
    {
        if(this == nullptr)
            return;
        for(uint i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] != nullptr)
                std::cout << " " << this->ptrs[i]->index;
        }
        std::cout << ((this->state & base_edge) ? ""  : ("[" + std::to_string(this->index) + "]"));
        return;
    }

    void print(void)
    {
        if(this->state & binary_cluster)
            std::cout << red;
        else if (this->state & unary_cluster)
            std::cout << cyan;
        else if (this->state & nullary_cluster)
            std::cout << green;
        else if (this->state & live)
            std::cout << bright_yellow;
        else
            std::cout << reset;
        std::cout << "Index: " << this->index << reset;
        std::cout << blue << "[" << (int) this->get_height()  << "]" << reset;
        std::cout << yellow << "[" << this->data << "]" << reset;
        std::cout << magenta << "[" << (int) this->counter << "]" << reset;
        std::cout << "  ";
        if(this->state&live)
            std::cout << "live ";
        else
            std::cout << "dead ";
        for(uint i = 0; i < this->size; i+=2)
        {
            auto cluster_ptr = this->ptrs[i];
            if(cluster_ptr == nullptr)
                continue;
            std::cout << cluster_ptr->index;
            std::cout << " ";
            auto ptr_type = this->types[i];
            if(ptr_type & neighbour_type)
                std::cout << "N";
            if(ptr_type & edge_type)
                std::cout << "E";
            if(ptr_type & parent_type)
                std::cout << "P";
            if(ptr_type & child_type)
                std::cout << "C";
            std::cout << "(";
            this->ptrs[i+1]->print_as_edge();
            std::cout << " ";
            ptr_type = this->types[i+1];
            if(ptr_type & neighbour_type)
                std::cout << "N";
            if(ptr_type & edge_type)
                std::cout << "E";
            if(ptr_type & parent_type)
                std::cout << "P";
            if(ptr_type & child_type)
                std::cout << "C";
            std::cout << " " << green;
            std::cout << this->ptrs[i+1]->data;
            std::cout << reset << " ";
            std::cout << ")   ";
        }
        std::cout << std::endl;
        std::cout << bright_magenta;
        for(const auto& i : this->adjacencies)
            std::cout << i << " ";
        std::cout << reset << std::endl;

        std::cout << bright_green;
        for(const auto& i : this->alternate_adjacencies)
            std::cout << i << " ";
        std::cout << reset << std::endl;

        


        // std::cout << " colour(" << this->colour << ") ";
        std::cout << std::endl;
    }
};



#endif