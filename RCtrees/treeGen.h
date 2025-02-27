#include "../include/parlay/sequence.h"
#include "../include/parlay/primitives.h"
#include <cmath>



enum distribution {
        constant,
        uniform,
        geometric,
        exponential
    };

template<typename T, typename D> // T is index type, D is data type
class TreeGen
{

    using subgraph = std::pair<T,T>; // subgraph
    using edge = std::tuple<T,T,D>;
    using edgelist = parlay::sequence<edge>;

    // parent subgraph, child subgraph, connecting edge, boolean
    using interGraphConnect = std::tuple<subgraph, subgraph, edge, bool>; 

    private:
    public:
        distribution dist;

        T num_vertices;

        double lg;
        double ln;
        double mean;
        double min_weight;
        double max_weight;
        long seed = 15213;

        parlay::sequence<subgraph> subgraphs;
        // parlay::sequence<std::atomic<short>> counts;
        parlay::sequence<T> parents;
        parlay::sequence<interGraphConnect> interconnects;
        parlay::sequence<D> weights;
        parlay::sequence<T> random_perm_map;
        // parlay::sequence<bool> exists;

        
        /*
            
        */
        TreeGen(T num_vertices, D min_weight = 0.0, D max_weight = 1.0, double l = 0.95f, double mean = 1000.0f, distribution dist = uniform, bool randomize_map = true, long seed = 15213)
        {
            assert(mean <= num_vertices); //otherwise subgraph assertion fails later
            assert(num_vertices > 0);
            assert(l > 0.0f && l < 1.0f); 
            assert(mean > 1.0f);
            double ln = 1.0; // probability of picking the immediate left node;
            double lg = l; // probability of picking the immediate left lg;
            this->min_weight = min_weight;
            this->max_weight = max_weight;

            this->dist = dist;
            this->ln = 1.0;
            this->lg = lg;
            this->mean = mean;
            this->num_vertices = num_vertices;
            this->seed = seed;
            if(randomize_map)
                this->random_perm_map = parlay::random_permutation(num_vertices);
            else
                this->random_perm_map = parlay::tabulate(num_vertices, [] (T i) {return i;});
        }


        void populateSubGraph(subgraph sg)
        {
            if(sg.first == sg.second)
            {
                parents[sg.first] = sg.first;                
                return;
            }
            
            parlay::random_generator gen(seed++);
            std::uniform_real_distribution<double> dis_ur(0, 1);

            parlay::parallel_for(sg.first, sg.second, [&] (T index) {
                std::uniform_int_distribution<T> dis(sg.first, index);
                auto r = gen[index];
                double random_val = dis_ur(r);
                T random_parent; 

                if(random_val <= ln)
                    random_parent = (index - 1) >= sg.first ? (index-1) : sg.first;
                else
                    random_parent = dis(r);
    

                if(random_parent == index)
                {
                    parents[index] = index;
                    return;
                }
                // try to increment myself
                // std::atomic<short>& my_count = counts[index];
                // std::atomic<short>& parent_count = counts[random_parent];

                // std::cout << "sgtrying " << index << " -- " << random_parent << std::endl;

                // short old_count = my_count.fetch_add(1);
                // if(old_count >= max_edges) // failure
                // {
                    // my_count.fetch_add(-1);
                    // parents[index] = index;
                    // return;
                // }
                // try to increment parent
                // short old_parent_count = parent_count.fetch_add(1);
                // if(old_parent_count >= max_edges)
                // {
                    // parent_count.fetch_add(-1);
                    // my_count.fetch_add(-1);
                    // parents[index] = index;
                    // return;
                // }
                // can add intergraph edge!

                parents[index] = random_parent;
                std::uniform_real_distribution<D> disint(this->min_weight, this->max_weight);
                D random_weight = disint(r);
                weights[index] = random_weight;
            });
        }

        void verifyParents()
        {
            // parlay::sequence<std::atomic<short>> actual_counts = parlay::sequence<std::atomic<short>>(this->num_vertices);

            // parlay::parallel_for(0, this->parents.size(), [&] (T i) {
            //     if(parents[i] != i)
            //     {
            //         actual_counts[i].fetch_add(1);
            //         actual_counts[parents[i]].fetch_add(1);
            //     }
            // });

            // parlay::parallel_for(0, this->parents.size(), [&] (T i) {
            //     if(this->counts[i] != actual_counts[i])
            //         std::cout << "[" << i << "] does not match " << this->counts[i] << " vs " << actual_counts[i] << std::endl;
            //     assert(this->counts[i] == actual_counts[i]);
            //     assert(this->counts[i] <= max_edges);
            // });
        }

        T get_num_dynamic_edges()
        {
            return this->subgraphs.size() > 0 ? this->subgraphs.size() - 1 : 0;
        }

        void generateInterconnects()
        {
            parlay::random_generator gen(seed++);
            std::uniform_real_distribution<double> dis_ur(0, 1);
            
            interconnects = parlay::tabulate(this->subgraphs.size(), [&] (T i) {
                if(i == 0)
                {
                    edge default_retedge = edge(i, i, (D) i); // if i == i, invalid edge which will get filtered out
                    interGraphConnect defretIGC = interGraphConnect(this->subgraphs[i], this->subgraphs[i], default_retedge, false); // invalid
                    return defretIGC;
                }
                    
                std::uniform_int_distribution<T> dis(0, std::numeric_limits<T>::max());
                auto r = gen[i];
                T parent_subgraph_index = dis(r);
                parent_subgraph_index = parent_subgraph_index % i; // 0 to i - 1 now hopefully
                
                double random_val = dis_ur(r);
                
                edge default_retedge = edge(i, i, (D) i); // if i == i, invalid edge which will get filtered out
                interGraphConnect defretIGC = interGraphConnect(this->subgraphs[i], this->subgraphs[i], default_retedge, false); // invalid
                
                

                if(random_val <= this->lg)
                {
                    parent_subgraph_index = i - 1;
                    if(parent_subgraph_index < 0)
                        return defretIGC;
                }
                
                if(parent_subgraph_index == i || parent_subgraph_index < 0)
                {
                    return defretIGC;
                }
                
                subgraph& my_subgraph = this->subgraphs[i];
                subgraph& parent_subgraph = this->subgraphs[parent_subgraph_index];
                
                T my_index = my_subgraph.first;

                assert(parent_subgraph.second - 1 >= parent_subgraph.first);
                assert(parent_subgraph.first >= 0);

                std::uniform_int_distribution<T> dis_parent(parent_subgraph.first, parent_subgraph.second - 1);
                
                T parent_index = dis_parent(r);

                assert(parents[my_index] == my_index);
                assert(parent_index != my_index);

                // std::cout << "Trying " << my_index << "/" << counts[my_index] << " -- " << parent_index << "/" << counts[parent_index] << std::endl;

                // attempt to connect child and parent
                // short old_value = counts[my_index].fetch_add(1);
                // if(old_value >= max_edges)
                // {
                    // counts[my_index].fetch_add(-1);
                    // return defretIGC;
                // }

                // short old_parent_value = counts[parent_index].fetch_add(1);
                // if(old_parent_value >= max_edges)
                // {
                    // counts[parent_index].fetch_add(-1);
                    // counts[my_index].fetch_add(-1);
                    // return defretIGC;
                // }
                
                parents[my_index] = parent_index;

                // generate a random weight
                std::uniform_real_distribution<D> disint(this->min_weight, this->max_weight);
                D random_weight = disint(r);
                weights[my_index] = random_weight;
                
                return interGraphConnect(my_subgraph, parent_subgraph, edge(my_index, parent_index, random_weight), true);

            });  

            interconnects = parlay::filter(interconnects, [] (interGraphConnect I) {
                return std::get<3>(I);
            });          
            

        }

        parlay::sequence<std::pair<T,T>> generateDeleteEdges(double prob = 1.0) // probability of
        {
            parlay::random_generator gen(seed++);
            std::uniform_real_distribution<double> dis_ur(0, 1);

            assert(prob <= 1.0 && prob >= 0.0);

            parlay::sequence<std::pair<T,T>> retpairs;
            retpairs = parlay::tabulate(this->interconnects.size(), [&] (T i) {
                auto r = gen[i];
                double random_roll = dis_ur(r);
                interGraphConnect& ICG = this->interconnects[i];
                bool& edge_live = std::get<3>(ICG);
                
                if(edge_live == false || random_roll > prob)
                {
                    return std::pair<T,T>(i, i);
                }
                else
                {
                    edge_live = false;
                    // counts[std::get<0>(std::get<2>(ICG))].fetch_add(-1);
                    // counts[std::get<1>(std::get<2>(ICG))].fetch_add(-1);
                    parents[std::get<0>(std::get<2>(ICG))] = std::get<0>(std::get<2>(ICG));
                    // return std::pair<T,T>(std::get<0>(std::get<2>(ICG)), std::get<1>(std::get<2>(ICG)));
                    return std::pair<T,T>(random_perm_map[std::get<0>(std::get<2>(ICG))], random_perm_map[std::get<1>(std::get<2>(ICG))]);
                }
            });
            retpairs = parlay::filter(retpairs, [&] (std::pair<T,T> pr) {
                return pr.first != pr.second;
            }); 

            return retpairs;
        }

        edgelist generateAddEdges(double prob = 1.0)
        {

            edgelist retedgelist;
            parlay::random_generator gen(seed++);
            std::uniform_real_distribution<double> dis_ur(0, 1);

            assert(prob <= 1.0 && prob >= 0.0);

            retedgelist = parlay::tabulate(this->interconnects.size(), [&] (T i) {
                interGraphConnect& ICG = this->interconnects[i];
                bool& exists = std::get<3>(ICG);
                auto r = gen[i];
                double random_roll = dis_ur(r);

                if(exists || random_roll > prob) // can't add it again
                {
                    return edge(i, i, (D) i); // invalid edge to itself
                }

                const T& my_index = std::get<0>(std::get<2>(ICG));
                T& parent_index = std::get<1>(std::get<2>(ICG));
                D& old_weight = std::get<2>(std::get<2>(ICG));

                // std::cout << "Trying " << my_index << " -- " << parent_index << std::endl;
                assert(parent_index != my_index);

                // decrement both myself and my parent count
                // counts[my_index].fetch_add(-1);
                // counts[parent_index].fetch_add(-1);
                // parents[my_index] = my_index;
                // exists = false;

                // Pick a random parent
                subgraph& parent_graph = std::get<1>(ICG);
                std::uniform_int_distribution<T> dis(parent_graph.first, parent_graph.second-1);
                
                T random_parent_index = dis(r);

                assert(random_parent_index != my_index);

                // Can I accept new counts
                // short old_my_val = counts[my_index].fetch_add(1);
                // std::cout << "Just incremented " << my_index << " counts from " << old_my_val << std::endl;
                // if(old_my_val >= max_edges)
                // {
                //     counts[my_index].fetch_add(-1);
                //     // std::cout << "Just decremented " << my_index << " counts" << std::endl;
                
                //     return edge(i, i, (D) i); // invalid edge to yourself
                // }

                // short old_parent_val = counts[random_parent_index].fetch_add(1);
                // std::cout << "Justt incremented " << random_parent_index << " counts" << std::endl;
                
                // if(old_parent_val >= max_edges)
                // {
                //     counts[random_parent_index].fetch_add(-1);
                //     // std::cout << "Justt decremented " << random_parent_index << " counts" << std::endl;
                
                //     counts[my_index].fetch_add(-1);
                //     // std::cout << "Justt decremented " << my_index << " counts" << std::endl;
                
                //     return edge(i, i, (D) i);
                // }

                parent_index = random_parent_index;

                std::uniform_real_distribution<D> disint(this->min_weight, this->max_weight);
                D random_weight = disint(r);
                old_weight = random_weight;

                parents[my_index] = parent_index;
                weights[my_index] = random_weight;

                exists = true;

                // return std::get<2>(ICG);
                edge edge_to_return = std::get<2>(ICG);
                std::get<0>(edge_to_return) = random_perm_map[std::get<0>(edge_to_return)];
                std::get<1>(edge_to_return) = random_perm_map[std::get<1>(edge_to_return)];
                // return std::get<2>(ICG);
                return edge_to_return;
            });

            retedgelist = parlay::filter(retedgelist, [] (edge E) {
                return std::get<0>(E) != std::get<1>(E);
            });


            return retedgelist;
            
        }

        void generateInitialEdges()
        {
            parlay::sequence<T> sizes;
            T total_size;
            if(dist == constant)
            {
                sizes = parlay::tabulate((T)((num_vertices)/(mean) + 2), [&] (T i) {
                    return (T) mean;
                });
                total_size = parlay::reduce(sizes);
            }
            else if(dist == uniform)
            {
                parlay::random_generator gen(seed++);
                std::uniform_int_distribution<T> dis(1, (T)(this->mean * 2));


                sizes = parlay::tabulate((T)(num_vertices/(mean)), [&] (T i) {
                    auto r = gen[i];
                    return dis(r);
                    // return (T) (i % 3);
                });
                total_size = parlay::reduce(sizes);
            }
            else if(dist == geometric)
            {
                parlay::random_generator gen(seed++); 
                std::geometric_distribution<T> dis(1/mean);
                // std::geometric_distribution<T> dis(0.01);

                sizes = parlay::tabulate((T)(num_vertices/(mean)), [&] (T i) {
                    auto r = gen[i];
                    T retval = dis(r);
                    return retval > 0 ? retval : 1;
                    // return (T) (i % 3);
                });
                total_size = parlay::reduce(sizes);
            }
            else if(dist == exponential)
            {
                parlay::random_generator gen(seed++); // CMU's postal code
                std::exponential_distribution<double> dis(1/mean);
                // std::geometric_distribution<T> dis(0.01);

                sizes = parlay::tabulate((T)(num_vertices/(mean)), [&] (T i) {
                    auto r = gen[i];
                    double retval = dis(r);
                    T retT = ceil(retval);
                    return retT;
                    // return (T) (i % 3);
                });
                total_size = parlay::reduce(sizes);
            }
            

            auto retpair = parlay::scan(sizes);
            auto& cumul_sizes = retpair.first;
            auto& actual_total_size = retpair.second;

            // for(auto& sz : sizes)
            //     std::cout << sz << " ";
            // std::cout << std::endl;
            // for(auto& sz : cumul_sizes)
            //     std::cout << sz << " ";
            // std::cout << std::endl;
            
            if(total_size < this->num_vertices)
            {
                // std::cout << "Undersampled, last subgraph is  " <<  (this->num_vertices) - total_size << " in size" << std::endl;
                sizes.push_back((this->num_vertices) - total_size); 
                actual_total_size+=(this->num_vertices) - total_size;
            }
            else
            {
                // std::cout << "Oversampled, last subgraph goes till " <<  total_size << " in size" << std::endl;
                if(sizes.size() == 1) // there is only one subgraph
                {
                    sizes[0]-= total_size-this->num_vertices;
                }
                else
                {

                    auto indices = parlay::tabulate(cumul_sizes.size(), [] (T i) {return i;});
                    // std::cout << "Num vertices " << this->num_vertices << std::endl;
                    // std::cout << "Last cumul " << cumul_sizes[cumul_sizes.size()-1] << std::endl;

                    T& last_cumul = cumul_sizes[cumul_sizes.size()-1];
                    if(last_cumul < this->num_vertices)
                    {
                        last_cumul+= (this->num_vertices - last_cumul);
                    }

                    auto location = parlay::find_if(indices, [&] (T i) {return cumul_sizes[i] >= this->num_vertices;});

                    assert(location != indices.end());
                    T index_which_exceeds = *location;
                    
                    T last_valid_index = index_which_exceeds;
                    // std::cout << "last valid index: " <<  *location << std::endl;
                    T overshoot = cumul_sizes[last_valid_index-1] - this->num_vertices;
                    // std::cout << "overshoot in last bin " << overshoot << std::endl;
                    if(overshoot == 0)
                        last_valid_index--;
                    
                    cumul_sizes = parlay::tabulate(last_valid_index, [&] (T i) {
                        return cumul_sizes[i];
                    });
                }
            }

            // std::cout << "Total subgraphs: " << cumul_sizes.size() << std::endl;

            this->subgraphs= parlay::tabulate(cumul_sizes.size(), [&] (T i) {
                if(i == (cumul_sizes.size() - 1))
                    return subgraph(cumul_sizes[i], this->num_vertices);
                return subgraph(cumul_sizes[i], cumul_sizes[i+1]);
            });

            // for(auto& pr : this->subgraphs)
            //     std::cout << pr.first << " " << pr.second << std::endl;
            
            parlay::parallel_for(0, subgraphs.size(), [&] (T i) { // I am fine with leaving asserts in code not being measured
                auto& sg = this->subgraphs[i];
                assert(sg.first != sg.second);
            });
            if(this->subgraphs.size() == 0)
            {
                this->subgraphs.push_back(subgraph(0,0));
            }

            if(subgraphs.size() == 1) // only one entry
            {
                auto& sg = this->subgraphs[0];
                sg.first = 0;
                sg.second = this->num_vertices;
            }
            assert(subgraphs[subgraphs.size()-1].second == this->num_vertices);
            
            // this->counts = parlay::sequence<std::atomic<short>>(this->num_vertices);
            this->parents = parlay::sequence<T>(this->num_vertices);
            this->weights = parlay::sequence<D>(this->num_vertices);

            parlay::parallel_for(0, this->subgraphs.size(), [&] (T i) {
                this->populateSubGraph(this->subgraphs[i]);
            });


            // for(unsigned int i = 0; i < parents.size(); i++)
            // {
            //     std::cout << i <<" -> " << parents[i] << std::endl; 
            // }

            verifyParents(); // TODO
            // std::cout << "initial fine, interconnect?" << std::endl;
            generateInterconnects();
        }

        void clearWeights(D value = 0.0)
        {
            parlay::parallel_for(0, this->weights.size(), [&] (T i) {
                this->weights[i] = value;
            });
        }

        parlay::sequence<edge> getAllEdges()
        {
            auto retedges = parlay::tabulate(num_vertices, [&] (T i) {
                // return edge(i, parents[i], weights[i]);
                return edge(random_perm_map[i], random_perm_map[parents[i]], weights[i]);
            });
            retedges = parlay::filter(retedges, [] (edge E) {
                return std::get<0>(E) != std::get<1>(E);
            });
            return retedges;
        }

};