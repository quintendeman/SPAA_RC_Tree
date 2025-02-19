#include <parlay/primitives.h>
#include <parlay/hash_table.h>
#include <iostream>
#include <sstream>

#include <cmath>


const double load_factor = 4.0f; // load factor for the hash table
const long extra_tern_node_factor = 3;

namespace { // we already have a tern_node in adjacency linked list but.. I couldn't come up with a better name
    

    template <typename T, typename D>
    struct tern_node{
        T dummy_index = -1; // the actual node that I am
        T v = -1; // The owner i.e. whose linked list I am in
        T w = -1; // the tern_node I am an extra edge to
        std::array<T, 3> outgoing_edgeindices = {-1};
        std::array<D, 3> outgoing_weights;
        double random_colour = 0.0;
        bool real = true; // TODO unnecessary?
        bool marked_for_deletion = false;
        bool shortlisted; // pun not intended

        void clear_data(D identity)
        {
            v = -1;
            w = -1;
            // assert(!real); // TODO remove
            marked_for_deletion = false;
            shortlisted = false;
            random_colour = 0.0f;
            outgoing_edgeindices[0] = outgoing_edgeindices[1] = outgoing_edgeindices[2] = -1;
            outgoing_weights[0] = outgoing_weights[1] = outgoing_weights[2] = identity;
        }

        friend std::ostream& operator<<(std::ostream& os, const tern_node& obj) {
            os << "[" << obj.v << "-" << obj.w << ": " << obj.dummy_index << " ";
            os << obj.outgoing_edgeindices[0] << "/" << obj.outgoing_weights[0] << " ";
            os << obj.outgoing_edgeindices[1] << "/" << obj.outgoing_weights[1] << " ";
            os << obj.outgoing_edgeindices[2] << "/" << obj.outgoing_weights[2] << " ";
            os << std::fixed << std::setprecision(2) << obj.random_colour << " ";
            if(obj.real)
                os << "REAL ";
            else
                os << "DMMY ";
            if(obj.marked_for_deletion)
                os << "M";
            if(obj.shortlisted)
                os << "S";
            os << "]";
            return os;
        }
    };

    template<typename T, typename D>
    struct hash_nodic {
        using kType = tern_node<T, D>; // Key type: the 3-element tuple
        using eType = void*; // Value type: pointer to the tuple

        // Define an empty value as nullptr
        eType empty() { return nullptr; }

        // Extract the key (tuple) from the value (pointer)
        kType getKey(eType v) {
            return *static_cast<kType*>(v); // Dereference pointer to get the tuple
        }

            // Hash function for the tuple
        size_t hash(kType V) {
            // Hash each element of the tuple and combine
            assert(V.v != -1 && V.w != -1); 
            T v = V.v;
            T w = V.w;
            if(v < w)
                std::swap(v, w);
            static const long random_prime1 = 1299059;
            static const long random_prime2 = 1149769;
            static const long random_prime3 = 602057;

            size_t h1 = v * random_prime1;
            size_t h2 = w * random_prime2;
            // std::cout << "[HASH] returning " << (h1 ^ h2) << " for " << v << " -- " << w << std::endl;

            return h1 ^ h2; // Combine hashes using XOR 
            // return v + w; // Combine hashes using XOR 
        }

        // Comparison function for tuples
        int cmp(kType v1, kType v2) {
            // Extract values
            T val1 = v1.v, wal1 = v1.w;
            T val2 = v2.v, wal2 = v2.w;



            // Ensure (val, wal) is always in sorted order
            if (val1 < wal1) std::swap(val1, wal1);
            if (val2 < wal2) std::swap(val2, wal2);

            // std::cout << "cmp" << std::endl;

            // Lexicographic comparison: (val1, wal1) vs (val2, wal2)
            if (val1 != val2) return (val1 < val2) ? -1 : 1;
            if (wal1 != wal2) return (wal1 < wal2) ? -1 : 1;
            return 0; // Equal
        }

        // Replace query: always return false (no replacement logic here)
        bool replaceQ(eType, eType) { return false; }

        // Update: return the current value
        eType update(eType v, eType) { return v; }

        // Atomic compare-and-swap (CAS) operation
        bool cas(eType* p, eType o, eType n) {
            return std::atomic_compare_exchange_strong_explicit(
            reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
        }
    };
}

template <typename T, typename D>
class ternarizer{
    using wedge = std::tuple<T,T,D>;
    private:
        T max_index;
        T num_vertices; // including dummy graphs
        D identity;
        parlay::sequence<T> free_list; // free dummy indices;
        parlay::sequence<tern_node<T,D>> simplified_tree;
        parlay::hashtable<hash_nodic<T,D>> HT;
        parlay::sequence<std::atomic<T>> tail_indices;
        long seed = 0;

        
        

        T cursor = 0; // location of free_list entries;
    public:

        ternarizer(T max_index, D identity) : HT(extra_tern_node_factor * max_index, hash_nodic<T,D>(), load_factor)
        {
            this->identity = identity;
            this->max_index = max_index;
            this->num_vertices = extra_tern_node_factor * max_index;
            this->free_list = parlay::tabulate(this->num_vertices - max_index, [max_index] (T i) {
                return i + max_index; // N to 2N-1
            });
            this->simplified_tree = parlay::tabulate(this->num_vertices, [max_index, identity] (T i) {
                auto ret_tern_node = tern_node<T,D>();
                ret_tern_node.clear_data(identity);
                ret_tern_node.dummy_index = i;
                ret_tern_node.real = i >= max_index ? false : true;
                ret_tern_node.v = i >= max_index ? -1 : i;
                return ret_tern_node;
            });
            // this->tail_indices = parlay::tabulate(this->max_index, [] (T i){
            //     return i;
            // });
            this->tail_indices = parlay::sequence<std::atomic<T>>(this->max_index);
            parlay::parallel_for(0, max_index, [&] (T i) {
                this->tail_indices[i] = i;
            });
            
            this->cursor = 0; // start at N  (this just is for the free list)
            // this->HT = parlay::hashtable<hash_nodic<T,D>> (num_vertices, hash_nodic<T,D>(), load_factor);
            return;
        }

        auto group_by_first(parlay::sequence<wedge> wedge_list)
        {
            auto keyed_sequence = parlay::delayed_tabulate(wedge_list.size(), [wedge_list] (T i) {
                return std::tuple<T,wedge>(std::get<0>(wedge_list[i]), wedge_list[i]);
            });

            return std::move(group_by_key(keyed_sequence));
        }

        auto group_by_second(parlay::sequence<wedge> wedge_list)
        {
            auto keyed_sequence = parlay::delayed_tabulate(wedge_list.size(), [wedge_list] (T i) {
                return std::tuple<T,wedge>(std::get<1>(wedge_list[i]), wedge_list[i]);
            });
            return std::move(group_by_key(keyed_sequence));
        }

        auto group_by_first(parlay::sequence<std::pair<T,T>> edge_list)
        {
            auto keyed_sequence = parlay::delayed_tabulate(edge_list.size(), [edge_list] (T i) {
                return std::tuple<T,std::pair<T,T>>(std::get<0>(edge_list[i]), edge_list[i]);
            });

            return std::move(group_by_key(keyed_sequence));
        }

        auto group_by_second(parlay::sequence<std::pair<T,T>> edge_list)
        {
            auto keyed_sequence = parlay::delayed_tabulate(edge_list.size(), [edge_list] (T i) {
                return std::tuple<T,std::pair<T,T>>(std::get<1>(edge_list[i]), edge_list[i]);
            });

            return std::move(group_by_key(keyed_sequence));
        }

        tern_node<T,D>& find_tern_node_ptr_in_HT(T v, T w) // smaller helper function
        {
            if(v < w)
                std::swap(v, w);
            tern_node<T,D> DN; // dummy tern_node
            DN.v = v;
            DN.w = w;
            void* ptr_to_tern_node = this->HT.find(DN);
            assert(ptr_to_tern_node != nullptr); // TODO remove
            return *static_cast<tern_node<T,D>*>(ptr_to_tern_node);
        }

        parlay::sequence<wedge> add_edges(parlay::sequence<wedge>& add_edges)
        {
            parlay::sequence<wedge> ret_modified_edges;
            parlay::parallel_for(0, add_edges.size(), [&add_edges] (T i) {
                auto& we = add_edges[i];
                T& v = std::get<0>(we);
                T& w = std::get<1>(we);
                if(v > w)
                    std::swap(v, w);
                return;
            });

            auto by_first = group_by_first(add_edges);
            auto by_second = group_by_second(add_edges);

            // if(add_edges.size() <= 500)
            // {
            //     std::cout << "Sorted by first: " << std::endl;;
            //     for(auto& og : by_first)
            //     {
            //         std::cout << og.first << ": ";
            //         for(wedge& we : og.second)
            //             std::cout << "[" << std::get<0>(we) << " " << std::get<1>(we) << " " << std::get<2>(we) << "] ";
            //         std::cout << std::endl;
            //     }

            //     std::cout << "Sorted by second: " << std::endl;
            //     for(auto& og : by_second)
            //     {
            //         std::cout << og.first << ": ";
            //         for(wedge& we : og.second)
            //             std::cout << "[" << std::get<0>(we) << " " << std::get<1>(we) << " " << std::get<2>(we) << "] ";
            //         std::cout << std::endl;
            //     }
            // } 

            parlay::sequence<wedge> add_pairs_gt;
            parlay::sequence<wedge> add_pairs_lt;

            {

                auto extra_counts = parlay::delayed_tabulate(by_first.size(), [&] (T i) {
                    auto& group = by_first[i];
                    T& index = group.first;
                    auto& edgelist = group.second;
                    T incoming_count = group.second.size();
                    T extra_count = incoming_count;
                    return extra_count;
                });

                auto count_pairs = parlay::scan(extra_counts);
                auto& entry_wise_counts = count_pairs.first;
                auto& total_entries = count_pairs.second;

                assert(this->cursor + total_entries <= num_vertices);

                add_pairs_gt = parlay::flatten(parlay::tabulate(by_first.size(), [&] (T I){
                    T my_extra_tern_node_index = this->cursor + entry_wise_counts[I];
                    // std::cout << "my_extra_tern_node_index " << my_extra_tern_node_index << std::endl;

                    auto& group = by_first[I];
                    T& index = group.first;
                    auto& edgelist = group.second;

                    auto& my_tern_node = this->simplified_tree[index];
                    
                    my_tern_node.real = true;

                    tern_node<T,D>& tail_tern_node = this->simplified_tree[this->tail_indices[index]];

                    assert(tail_tern_node.outgoing_edgeindices[2] == -1); // TODO

                    parlay::sequence<wedge> contributions;

                    contributions = parlay::flatten(parlay::delayed_tabulate(edgelist.size(), [&] (T i) {

                        parlay::short_sequence<wedge> local_contribution;
                        T free_tern_node_index = this->free_list[my_extra_tern_node_index + i];
                        
                        // std::cout << "Free tern_node index " << free_tern_node_index << std::endl;
                        tern_node<T,D>& dummy_tern_node = this->simplified_tree[free_tern_node_index];
                        dummy_tern_node.real = false;

                        T& v = std::get<0>(edgelist[i]);
                        T& w = std::get<1>(edgelist[i]);
                        D& weight = std::get<2>(edgelist[i]);
                        dummy_tern_node.v = index;
                        dummy_tern_node.w = w;
                        // std::cout << "free_tern_node_index: " << dummy_tern_node.dummy_index << " setting " << dummy_tern_node.v << "  " << dummy_tern_node.w << std::endl;
                        // std::cout << "Inserting " << v << " -- " << w  << " at dummy index " << dummy_tern_node.dummy_index << std::endl;
                        this->HT.insert(&dummy_tern_node);

                       

                        if(i == 0) // INSERTING AN EDGE
                        { // make sure it goes to the tail
                            dummy_tern_node.outgoing_edgeindices[0] = tail_tern_node.dummy_index;
                            dummy_tern_node.outgoing_weights[0] = identity;
                            tail_tern_node.outgoing_edgeindices[2] = free_tern_node_index;
                            tail_tern_node.outgoing_weights[2] = identity;
                            local_contribution.push_back(std::tuple<T,T,D>(tail_tern_node.dummy_index, free_tern_node_index, identity));
                        }
                        else // INSERTING AN EDGE
                        {
                            T prev_free_tern_node_index = this->free_list[my_extra_tern_node_index + i - 1];
                            dummy_tern_node.outgoing_edgeindices[0] = prev_free_tern_node_index;
                            dummy_tern_node.outgoing_weights[0] = identity;
                            this->simplified_tree[prev_free_tern_node_index].outgoing_edgeindices[2] = free_tern_node_index;
                            this->simplified_tree[prev_free_tern_node_index].outgoing_weights[2] = identity;
                            local_contribution.push_back(std::tuple<T,T,D>(prev_free_tern_node_index, free_tern_node_index, identity));
                        }
                        dummy_tern_node.outgoing_edgeindices[1] = -1; // unknown for now
                        dummy_tern_node.outgoing_weights[1] = weight;

                        if(i == (edgelist.size() - 1)) // last one
                        {
                            dummy_tern_node.outgoing_edgeindices[2] = -1;
                            this->tail_indices[index] = free_tern_node_index;
                        }

                        return local_contribution;
                    }));

                    return contributions;
                }));

                // std::cout << "cursor went from " << this->cursor << " to ";
                this->cursor += total_entries; 
                // std::cout << this->cursor << std::endl;
            }

            {
                auto extra_counts = parlay::delayed_tabulate(by_second.size(), [&] (T i) {
                    auto& group = by_second[i];
                    T& index = group.first;
                    auto& edgelist = group.second;
                    T incoming_count = group.second.size();
                    T extra_count = incoming_count;
                    return extra_count;
                });

                auto count_pairs = parlay::scan(extra_counts);
                auto& entry_wise_counts = count_pairs.first;
                auto& total_entries = count_pairs.second;

                assert(this->cursor + total_entries <= num_vertices);

                

                add_pairs_lt = parlay::flatten(parlay::tabulate(by_second.size(), [&] (T I){
                    T my_extra_tern_node_index = this->cursor + entry_wise_counts[I];
                    // std::cout << "my_extra_tern_node_index " << my_extra_tern_node_index << std::endl;

                    auto& group = by_second[I];
                    T& index = group.first;
                    auto& edgelist = group.second;

                    auto& my_tern_node = this->simplified_tree[index];

                    my_tern_node.real = true;

                    tern_node<T,D>& tail_tern_node = this->simplified_tree[this->tail_indices[index]];
                    

                    parlay::sequence<wedge> contributions;

                    contributions = parlay::flatten(parlay::delayed_tabulate(edgelist.size(), [&] (T i) {
                        parlay::short_sequence<wedge> local_contribution;
                        T free_tern_node_index = this->free_list[my_extra_tern_node_index + i];
                        
                        auto& dummy_tern_node = this->simplified_tree[free_tern_node_index];


                        dummy_tern_node.real = false;

                        T& v = std::get<1>(edgelist[i]);
                        T& w = std::get<0>(edgelist[i]);
                        D& weight = std::get<2>(edgelist[i]);
                        dummy_tern_node.v = index;
                        dummy_tern_node.w = w;
                        

                        if(i == 0) // INSERTING AN EDGE
                        {
                            dummy_tern_node.outgoing_edgeindices[0] = tail_tern_node.dummy_index;
                            dummy_tern_node.outgoing_weights[0] = identity;
                            tail_tern_node.outgoing_edgeindices[2] = free_tern_node_index;
                            tail_tern_node.outgoing_weights[2] = identity;
                            local_contribution.push_back(std::tuple<T,T,D>(tail_tern_node.dummy_index, free_tern_node_index, identity));
                        }
                        else // INSERTING AN EDGE
                        {
                            T prev_free_tern_node_index = this->free_list[my_extra_tern_node_index + i - 1];
                            dummy_tern_node.outgoing_edgeindices[0] = prev_free_tern_node_index;
                            dummy_tern_node.outgoing_weights[0] = identity;
                            this->simplified_tree[prev_free_tern_node_index].outgoing_edgeindices[2] = free_tern_node_index;
                            this->simplified_tree[prev_free_tern_node_index].outgoing_weights[2] = identity;
                            local_contribution.push_back(std::tuple<T,T,D>(prev_free_tern_node_index, free_tern_node_index, identity));
                        }
                        // dummy_tern_node.outgoing_edgeindices[1] = -1; // unknown for now

                        dummy_tern_node.outgoing_weights[1] = weight;
                        // std::cout << "attempting to look for " << v << " -- " << w << "... ";
                        auto& other_tern_node = this->find_tern_node_ptr_in_HT(v, w); // same v and w means this should match to the other side
                        // std::cout << " got " << other_tern_node.dummy_index << std::endl;
                        // std::cout << "Inserting " << v << " -- " << w  << " at dummy index " << dummy_tern_node.dummy_index << std::endl;
                        dummy_tern_node.outgoing_edgeindices[1] = other_tern_node.dummy_index; // INSERTING AN EDGE
                        other_tern_node.outgoing_edgeindices[1] = dummy_tern_node.dummy_index;

                        local_contribution.push_back(std::tuple<T,T,D>(other_tern_node.dummy_index, dummy_tern_node.dummy_index, weight));

                        if(i == (edgelist.size() - 1)) // last one
                        {
                            dummy_tern_node.outgoing_edgeindices[2] = -1;
                            this->tail_indices[index] = free_tern_node_index;
                        }

                        return local_contribution;
                    }));
                    return contributions;
                }));

                // std::cout << "cursor went from " << this->cursor << " to ";
                this->cursor += total_entries; 
                // std::cout << this->cursor << std::endl; 
            }

            

            auto& final_add_edges = ret_modified_edges;
            final_add_edges.append(add_pairs_lt);
            final_add_edges.append(add_pairs_gt);
            // parlay::append(std::move(final_add_edges), std::move(add_pairs_lt));
            // parlay::append(std::move(final_add_edges), std::move(add_pairs_gt));


            // parlay::append(std::move(del_pairs), std::move(del_pairs_lt));
            // parlay::append(std::move(del_pairs), std::move(del_pairs_gt));

            return std::move(ret_modified_edges);
        }

        std::pair<parlay::sequence<wedge>,parlay::sequence<std::pair<T,T>>> delete_edges(parlay::sequence<std::pair<T,T>> input_edges)
        {
            this->print_state();
            
            auto end_points = parlay::flatten(parlay::delayed_tabulate(input_edges.size(), [&] (T i) {
                parlay::sequence<tern_node<T,D>*> tern_node_ptrs;
                auto& edge = input_edges[i];
                tern_node<T,D>& first_end_point = this->find_tern_node_ptr_in_HT(edge.first, edge.second);
                tern_node<T,D>& second_end_point = this->simplified_tree[first_end_point.outgoing_edgeindices[1]];
                tern_node_ptrs.push_back(&first_end_point);
                tern_node_ptrs.push_back(&second_end_point);
                first_end_point.marked_for_deletion = true;
                first_end_point.shortlisted = false;
                second_end_point.marked_for_deletion = true;
                second_end_point.shortlisted = false;
                T free_list_entry_first_index = this->cursor - (2 * input_edges.size()) + (2 * i);
                T free_list_entry_second_index = this->cursor - (2 * input_edges.size()) + (2 * i + 1);
                this->free_list[free_list_entry_first_index] = first_end_point.dummy_index;
                this->free_list[free_list_entry_second_index] = second_end_point.dummy_index;
                return tern_node_ptrs;
            }));

            // std::cout << "cursor went from " << this->cursor << " to ";
            this->cursor-= (2 * input_edges.size());
            // std::cout << this->cursor << std::endl;
            assert(this->cursor >= 0);

            auto all_delete_edges = parlay::flatten(parlay::delayed_tabulate(end_points.size(), [&] (T i) {
                auto& my_node = *end_points[i];
                parlay::sequence<std::pair<T,T>> ret_seq;
                // left node always
                ret_seq.push_back(std::pair<T,T>(my_node.outgoing_edgeindices[0], my_node.dummy_index));

                // payload if tiebreak
                if(my_node.outgoing_edgeindices[1] < my_node.dummy_index)
                {
                    this->HT.deleteVal(my_node);
                    ret_seq.push_back(std::pair<T,T>(my_node.outgoing_edgeindices[1], my_node.dummy_index));
                }
                // right-most edge if it is not already marked for deletion and I am not tail
                if(my_node.outgoing_edgeindices[2] != -1 && !this->simplified_tree[my_node.outgoing_edgeindices[2]].marked_for_deletion) 
                    ret_seq.push_back(std::pair<T,T>(my_node.outgoing_edgeindices[2], my_node.dummy_index));
                
                return ret_seq;
            }));

            // start accumulating the add edges
            parlay::sequence<wedge> add_edges;

            // gather all the tails
            auto all_tails = parlay::filter(end_points, [&] (auto* tnp){
                auto& my_node = *tnp;
                if(this->tail_indices[my_node.v] == my_node.dummy_index)
                    return true;
                return false;
            });

            end_points = parlay::filter(end_points, [&] (auto* tnp){
                auto& my_node = *tnp;
                if(this->tail_indices[my_node.v] == my_node.dummy_index)
                    return false;
                return true;
            });

            // std::cout << "Gathered " << all_tails.size() << " tails and " << end_points.size() << " others " << std::endl;  
            


            int counter = 0;
            std::uniform_real_distribution<double> dis_ur(0, 1);
            while(end_points.size())
            {
                parlay::random_generator gen(15213 + seed++ + this->max_index);
                counter++;
                assert(counter < 10000); // To prevent, reasonably infinite loops
                // for(auto& ptr : end_points)
                //     std::cout << " " << *ptr;
                // std::cout << std::endl;
                // this->print_state();

                // std::cout << "end_points.size() = " << end_points.size() << std::endl;
                // randomly pick sqrt(n)
                const double max_prob = 0.5;
                const double min_prob = 0.4;
                double n = static_cast<double>(end_points.size());
                double probability =  (std::sqrt(n) * std::log(n) + 1)/n; // oversample
                if(probability > max_prob) // 
                    probability = max_prob;
                if(probability < min_prob)
                    probability = min_prob;
                
                


                

                parlay::parallel_for(0, end_points.size(), [&] (T i) {
                    auto r = gen[i];
                    double random_val = dis_ur(r);
                    end_points[i]->random_colour = random_val;
                    end_points[i]->shortlisted = false;
                });

                parlay::parallel_for(0, end_points.size(), [&] (T i) {
                    auto& my_node = *end_points[i];
                    const T& left = my_node.outgoing_edgeindices[0];
                    const T& right = my_node.outgoing_edgeindices[2];
                    double left_colour = 0.0f;
                    double right_colour = 0.0f;
                    if(left != -1)
                        left_colour = this->simplified_tree[left].random_colour;
                    if(right != -1)
                        right_colour = this->simplified_tree[right].random_colour;
                    const double& my_colour = my_node.random_colour;
                    if(my_colour > left_colour && my_colour > right_colour)
                        end_points[i]->shortlisted = true;

                });

                

                // parlay::sequence<double> colours = parlay::tabulate(end_points_shortlisted.size(), [] (T i) {
                //     auto r = gen[i];
                //     return dis_ur(r);
                // });

                // std::cout << "end_points_shortlisted.size() = " << end_points_shortlisted.size() << std::endl;

                auto add_edge_contributions = parlay::tabulate(end_points.size(), [&] (T i) {
                    auto ret_tuple = std::tuple<T,T,D>(0,0, this->identity); // empty by default
                    T& outv = std::get<0>(ret_tuple);
                    T& outw = std::get<1>(ret_tuple);

                    auto& my_node = *end_points[i];

                    assert(my_node.marked_for_deletion); // TODO remove
                    assert(!my_node.real);

                    if(!my_node.shortlisted)
                        return ret_tuple;

                    // // am I in the IS
                    T left = my_node.outgoing_edgeindices[0];
                    T right = my_node.outgoing_edgeindices[2];

                    assert(left != -1); // TODO remove

                    // double left_colour = this->simplified_tree[left].random_colour;
                    // double right_colour = right == -1 ? 0.0f : this->simplified_tree[right].random_colour;
                    // if(left_colour >= my_node.random_colour || right_colour >= my_node.random_colour)
                    //     return ret_tuple; // I wasn't picked :(
                    
                    auto& left_node = this->simplified_tree[left];

                    if(!left_node.marked_for_deletion) // left_node has finished
                    {
                        if(right == -1) // right is the end
                        {
                            assert(this->tail_indices[my_node.v] == my_node.dummy_index); // TODO remove
                            return ret_tuple;
                            this->tail_indices[my_node.v] = left;
                            left_node.outgoing_edgeindices[2] = -1;
                            // left_node.v = my_node.v;
                            std::ostringstream ss;
                            // ss << "deleting " << my_node << " left=" << left << "\n";
                            // std::cout << ss.str();
                            my_node.clear_data(this->identity);
                            return ret_tuple; // I am done, I contracted into nothingness
                        }
                        auto& right_node = this->simplified_tree[right];
                        if(!right_node.marked_for_deletion) // I am in the middle, connect the two edges and leave
                        {
                            right_node.outgoing_edgeindices[0] = left;
                            left_node.outgoing_edgeindices[2] = right;
                            outv = left;
                            outw = right;
                            // std::cout << ("splicing " + std::to_string(left) + "-" + std::to_string(my_node.dummy_index) + "-" + std::to_string(right) + "\n");
                            // std::cout << "deleting " << my_node << std::endl;
                            my_node.clear_data(this->identity);
                            return ret_tuple;
                        }
                    }

                    if(right == -1)
                        return ret_tuple;
                    auto& right_node = this->simplified_tree[right];
                    right_node.outgoing_edgeindices[0] = left;
                    left_node.outgoing_edgeindices[2] = right;
                    // std::cout << ("splicing " + std::to_string(left) + "-" + std::to_string(my_node.dummy_index) + "-" + std::to_string(right) + "\n");
                    // std::cout << "deleting " << my_node << std::endl;
                    my_node.clear_data(this->identity);
                    return ret_tuple;

                });

                    end_points = parlay::filter(end_points, [] (auto* tnp){
                        return tnp->marked_for_deletion;
                    });

                    add_edge_contributions = parlay::filter(add_edge_contributions, [] (auto ent){
                        return std::get<0>(ent) != std::get<1>(ent);
    
                });
                add_edges.append(add_edge_contributions);
            }

            parlay::parallel_for(0, all_tails.size(), [&] (T i){
                auto& my_node = *all_tails[i];
                assert(this->tail_indices[my_node.v] == my_node.dummy_index); // TODO remove both
                assert(my_node.outgoing_edgeindices[2] == -1);
                assert(simplified_tree[my_node.outgoing_edgeindices[0]].outgoing_edgeindices[2] = my_node.dummy_index);

                this->tail_indices[my_node.v] = my_node.outgoing_edgeindices[0];
                this->simplified_tree[my_node.outgoing_edgeindices[0]].outgoing_edgeindices[2] = -1;
                // std::cout << ("Setting tail of " + std::to_string(my_node.v) + " to " + std::to_string(my_node.outgoing_edgeindices[0]) + " from " + std::to_string(my_node.dummy_index) + " \n ");
            });

            parlay::parallel_for(0, all_tails.size(), [&] (T i){
                all_tails[i]->clear_data(this->identity);
            });

            this->print_state();
            return std::make_pair(add_edges, all_delete_edges);
        }

        void verify_simple_tree(void)
        {
            // do my neighbours have me too?
            parlay::parallel_for(0, this->num_vertices, [&] (T i) {
                auto& my_tern_node = this->simplified_tree[i];
                for(short I = 0; I < 3; I ++)
                {
                    T& neighbour_index = my_tern_node.outgoing_edgeindices[I];
                    if(neighbour_index < 0)
                        continue;
                    auto& neighbour_tern_node = this->simplified_tree[neighbour_index];
                    for(short j = 0; j < 3; j++)
                    {
                        if(neighbour_tern_node.outgoing_edgeindices[j] == my_tern_node.dummy_index)
                            return;
                    }
                    assert("Neighbour doesn't have me?" && false);
                }
            });
            parlay::parallel_for(0, this->num_vertices, [&] (T i) {
                auto& my_node = this->simplified_tree[i];
                if(my_node.v != -1)
                {
                    T v = my_node.v;
                    assert(this->simplified_tree[this->tail_indices[v]].outgoing_edgeindices[2] == -1);
                }
            });
            parlay::parallel_for(0, this->max_index, [&] (T i) {
                long tail_index = this->tail_indices[i];
                auto& tail_node = this->simplified_tree[tail_index];
                if(tail_node.outgoing_edgeindices[2] != -1)
                {
                    std::cout << tail_node << std::endl;
                    std::cout << this->simplified_tree[tail_node.outgoing_edgeindices[2]] << std::endl;
                }
                assert(tail_node.outgoing_edgeindices[2] == -1);
                if(tail_node.v != i)
                    std::cout << tail_node << " " << i << std::endl;
                assert(tail_node.v == i);
            });
        }

        void print_state(void)
        {
            return;
            for(long i = 0; i < this->max_index; i++)
            {
                std::cout << i << " " << this->tail_indices[i] << std::endl;
            }

            for(auto& n : this->simplified_tree)
            { 
                std::cout << n << std::endl;
            }
            return;
        }

};