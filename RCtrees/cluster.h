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


const std::string reset = "\033[0m";
const std::string black = "\033[30m";
const std::string red = "\033[31m";
const std::string green = "\033[32m";
const std::string yellow = "\033[33m";
const std::string blue = "\033[34m";
const std::string magenta = "\033[35m";
const std::string cyan = "\033[36m";
const std::string white = "\033[37m";


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
    cluster<T, D>* ptrs[max_neighbours * 2]; 
    T index = -1;
    T colour = -1;
    T initial_ngbrs[3];
    T height;
    D data = -1.0;
    static const short size = max_neighbours*2;
    short state = 0; // unaffected, affected or update eligible
    bool is_MIS = false;
    std::atomic<char> counter; // a node will not have more than 255 edges
    char types[max_neighbours * 2] = {};
    unsigned char contraction_time = 0;
    
    // Default constructor
    cluster() : height(0), counter(0) {
        for(uint i = 0; i < this->size; i++)
        {
            this->ptrs[i] = nullptr;
            this->types[i] = 0;
            if(i < max_neighbours)
                this->initial_ngbrs[i] = -1;
        }
    }

    // Copy constructor only for initialization when creating a sequence
    cluster(const cluster& other)
        : index(other.index),
          colour(other.colour),
          data(other.data),
          height(other.height), // Load the atomic value
          counter(other.counter.load()), // Load the atomic value
          state(other.state),
          contraction_time(other.contraction_time),
          is_MIS(other.is_MIS)
    {    
        for(uint i = 0; i < this->size; i++)
        {
            this->ptrs[i] = nullptr;
            this->types[i] = 0;
            if(i < max_neighbours)
                this->initial_ngbrs[i] = other.initial_ngbrs[i];
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

    cluster<T, D>* get_parent(void)
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

    // set the neighbours MIS to val
    // useful when colouring
    void set_neighbour_mis(bool val)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i])
                this->ptrs[i]->is_MIS = val;
        }
    }

    void set_neighbour_colour(T input_colour)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i])
                this->ptrs[i]->colour = input_colour;
        }   
    }

    bool is_max_neighbour_colour(void)
    {
        T max_colour = this->colour;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->ptrs[i]->colour > max_colour)
                max_colour = this->ptrs[i]->colour;
        }
        if(max_colour == this->colour)
            return true;
        return false;
    }

    /*
        return true if any neighbours are in MIS
        else return false
    */
    bool get_neighbour_MIS(void)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->ptrs[i]->is_MIS)
                return true;
        }
        return false;
    }

    short get_neighbour_count(void)
    {
        short num_neighbours = 0;
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i] & neighbour_type)
                num_neighbours++;
        }

        return num_neighbours;
    }

    short get_children_count(void)
    {
        short num_children;
        for(const auto& t : this->types)
            if (t & child_type)
                num_children++;
        return num_children;
    }

    void get_two_neighbours_edges(cluster<T, D>*& neighbour1, cluster<T, D>*& edge1, cluster<T, D>*& neighbour2, cluster<T, D>*& edge2)
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
    void find_boundary_vertices(cluster<T, D>* &l, D& lval, cluster<T, D>* &r, D& rval, const D& defretval)
    {
        if(this->state & nullary_cluster)
        {
            l = r = nullptr;
            return;
        }
        if(this->state & unary_cluster)
        {
            l = r = this->get_parent();
            /* l CAN NEVER BE NULLPTR*/
            lval = rval = l->data;
            return;
        }
        if(this->state & binary_cluster)
        {
            r = this->get_parent();
            rval = r->data;
            for(uint i = 0; i < this->size; i+=2)
            {
                if(this->types[i] & neighbour_type && this->ptrs[i] != r)
                {
                    l = this->ptrs[i];
                    lval = l->data;
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

    inline T get_height()
    {
        return this->height;
    }

    void set_height()
    {
        this->height = 0;
        for(uint i = 0; i < this->size; i+=2)
        {
            if(this->types[i] & child_type && this->ptrs[i] != nullptr)
                this->height = std::max(this->height, this->ptrs[i]->height);
        }
        this->height++;
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
        std::cout << "Index: " << this->index;
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
        // std::cout << " colour(" << this->colour << ") ";
        std::cout << std::endl;
    }
};
