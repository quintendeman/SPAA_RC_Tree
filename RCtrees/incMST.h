#ifndef INCMST_H
#define INCMST_H


#include "RC.h"
#include <queue>
#include <set>
#include "../examples/filter_kruskal.h"
#include "RCdynamic.h"


template<typename T, typename D> 
void verifyCompressPathTree(parlay::sequence<cluster<T,D>>& clusters, const parlay::sequence<std::tuple<T,T,D>>& newEdges, const parlay::sequence<T>& vertices);


// Takes a RC tree and generates a sequence of truples (without duplicates) of <Index, Index, Weight>
template <typename T, typename D>
parlay::sequence<std::tuple<T,T,D>> convertClustersToTruples(const parlay::sequence<cluster<T,D>>& clusters)
{

    parlay::sequence<std::tuple<T,T,D>> Truples;

    Truples = parlay::flatten(parlay::tabulate(clusters.size(), [&] (T index) {
        parlay::sequence<std::tuple<T,T,D>> Truple_contribution;

        const auto& relevant_cluster = clusters[index];
        for(const auto& edge : relevant_cluster.adjacency.get_head()->adjacents)
        {
            if(edge == nullptr)
                continue;
            const auto& other_side = get_other_side(relevant_cluster.adjacency.get_head(), edge);
            const auto& other_base_cluster_ptr = other_side->cluster_ptr;
            const auto& base_edge_cluster_ptr = edge->cluster_ptr;

            // Only consider the other_edge if the other side has an index greater than mine
            if(other_base_cluster_ptr->index > relevant_cluster.index)
                Truple_contribution.push_back(std::tuple<T,T,D>(other_base_cluster_ptr->index, relevant_cluster.index, base_edge_cluster_ptr->data));
        }
        return Truple_contribution;
    }));

    return Truples;
}

template<typename T>
void set_insert(std::vector<T>& vertices, T& new_vertex)
{
    assert(false && "Should no longer be called!");
    for(auto& ver : vertices)
        if(ver == new_vertex)
            return;
    vertices.push_back(new_vertex);
}



template<typename T, typename D>
void set_insert(std::vector<std::tuple<T,T,D>>& edges, std::tuple<T,T,D>& new_edge)
{
    assert(false && "Should no longer be called!");
    for(auto& edge: edges)
    {
        if(std::get<0>(edge) == std::get<0>(new_edge))
            if(std::get<1>(edge) == std::get<1>(new_edge))
                if(std::get<2>(edge) == std::get<2>(new_edge))
                    return;
    }
    edges.push_back(new_edge);
}

template<typename T>
void set_remove(std::vector<T>& vertices, const T& target_vertex)
{
    assert(false && "Should no longer be called!");
    for (auto it = vertices.begin(); it != vertices.end(); ++it)
    {
        if (*it == target_vertex)
        {
            vertices.erase(it);
            return;
        }
    }
}

template<typename T, typename D>
void set_remove(std::vector<std::tuple<T, T, D>>& edges, const std::tuple<T, T, D>& target_edge)
{
    assert(false && "Should no longer be called!");
    for (auto it = edges.begin(); it != edges.end(); ++it)
    {
        if (std::get<0>(*it) == std::get<0>(target_edge) &&
            std::get<1>(*it) == std::get<1>(target_edge) &&
            std::get<2>(*it) == std::get<2>(target_edge))
        {
            edges.erase(it);
            return;
        }
    }
}

// TODO: determine if a particular cluster is contributing a vertex or not and also create edges
template<typename T, typename D>
void  top_down_ct(cluster<T,D>*& cluster_ptr, parlay::sequence<cluster<T,D>>& clusters) // I call this on roots
{
    if(cluster_ptr == nullptr)
        return;
    cluster_ptr->vertex_count = 0;

    if(cluster_ptr->state & is_marked_endpoint)
    {
        if(!(cluster_ptr->state & is_marked))
        {
            cluster_ptr->print();
        }
        assert(cluster_ptr->state & is_marked);
    }

    if(!(cluster_ptr->state & is_marked))
    {
        if(isBinary(cluster_ptr)) // line 7 
        {
            cluster<T,D>* lptr = nullptr;
            cluster<T,D>* rptr = nullptr;
            D lval;
            D rval;
            D defval;
            cluster_ptr->find_boundary_vertices(lptr, lval, rptr, rval, defval);
            assert(lptr && rptr && lptr != rptr);
            T lindex = lptr->index;
            T rindex = rptr->index;
            D retval;
            if(lval > rval) // let dce handle it
                retval = lval;
            else
                retval = rval;
            retval = cluster_ptr->data;
            lptr->insertCGentry(rindex, retval, cluster_ptr->max_weight_edge);
            rptr->insertCGentry(lindex, retval, cluster_ptr->max_weight_edge);
            if(lindex > rindex) // assert some ordering
                std::swap(lindex, rindex);
            cluster_ptr->vertex_count = 0;
        }
        else if(isUnary(cluster_ptr) || isNullary(cluster_ptr)) // line 9
        {
            cluster<T,D>* lptr = nullptr;
            cluster<T,D>* rptr = nullptr;
            D lval;
            D rval;
            D defval;

            cluster_ptr->find_boundary_vertices(lptr, lval, rptr, rval, defval);
            if(!(lptr && rptr && lptr == rptr))
            {
                std::cout << cluster_ptr->index << std::endl;
                std::cout << isUnary(cluster_ptr) << " " << isNullary(cluster_ptr) << std::endl;
                cluster_ptr->print();
            }
            assert(lptr && rptr && lptr == rptr);
            T lindex = lptr->index;
            T rindex = rptr->index;
            // ret_vertices.push_back(lindex);
            cluster_ptr->vertex_count = 0; // you can't actually have a nullary child so this doesn't matter

        }
        else if(cluster_ptr->parent && (cluster_ptr->parent->state & is_marked)) // part of line 13 -- must be base edge by elimination
        { // Do we even need is marked here?? I suspect not since parent must have been marked for recursion to occur
            assert(cluster_ptr->state & base_edge); // TODO remove asserts
            cluster<T,D>* lptr = nullptr;
            cluster<T,D>* rptr = nullptr;
            cluster_ptr->find_endpoints(lptr, rptr);
            D retval = cluster_ptr->data;
            assert(lptr && rptr && lptr != rptr);
            T lindex = lptr->index;
            T rindex = rptr->index;
            lptr->insertCGentry(rindex, retval, cluster_ptr); // should just be cluster_ptr as base edge
            rptr->insertCGentry(lindex, retval, cluster_ptr);
            if(lindex > rindex) // assert some ordering
                std::swap(lindex, rindex); 
            cluster_ptr->vertex_count = 0;
        }
        return;
    }
    
    cluster_ptr->vertex_count = 1; // my  representative vertex
    
    parlay::parallel_for(0, max_neighbours, [&] (T i) { 
        auto& child_ptr = cluster_ptr->children[i];
        top_down_ct(child_ptr, clusters);
        return;
    });
    

    T index1 = -1;
    T index2 = -2;
    T degree = 0;
    T curindex = 0;
    T myIndex = cluster_ptr->index;
    
    for(auto& valid : cluster_ptr->CGEntryValid)
    {
        if(valid)
        {
            if(index1 == -1)
                index1 = curindex;
            else if(index2 == -1)
                index2 = curindex;
            degree++;
        }
        curindex++;
    }

    T to_decrease_edgecount = 0;

    T toRemove = myIndex;
    
    if(degree == 1 && !(cluster_ptr->state & is_marked_endpoint))
    {
        auto EntryToRemove = cluster_ptr->CGGraphEntry[index1];
        auto vertex_to_remove = myIndex;
        auto lptr = &clusters[EntryToRemove.first];
        lptr->clearCGgraphEntry(myIndex);
        auto rptr = &clusters[myIndex];
        rptr->clearCGgraphEntry(EntryToRemove.first);
        
        // set_remove(ret_vertices, vertex_to_remove);
        
        cluster_ptr->vertex_count = 0;

        toRemove = EntryToRemove.first;

    }

    // check if degree is 2 (again)
    degree = 0;
    curindex = 0;
    index1 = -1;
    index2 = -1;
    myIndex = toRemove;

    auto old_cluster_ptr = cluster_ptr;

    cluster_ptr = &clusters[myIndex];

    for(auto& valid : cluster_ptr->CGEntryValid)
    {
        if(valid)
        {
            if(index1 == -1)
                index1 = curindex;
            else if(index2 == -1)
                index2 = curindex;
            degree++;
        }
        curindex++;
    }

    if(degree == 2 && !(cluster_ptr->state & is_marked_endpoint)) //
    {

        auto relevant_edge1 = cluster_ptr->CGGraphEntry[index1];
        auto relevant_edge2 = cluster_ptr->CGGraphEntry[index2];
        D new_weight = relevant_edge1.second > relevant_edge2.second ? relevant_edge1.second : relevant_edge2.second;
        cluster<T,D>* new_heaviest_edge = relevant_edge1.second > relevant_edge2.second ? cluster_ptr->CGGraphEdgePtr[index1] : cluster_ptr->CGGraphEdgePtr[index2];
        T newlindex = relevant_edge1.first;
        T newrindex = relevant_edge2.first;
        auto lptr = &clusters[newlindex];
        lptr->clearCGgraphEntry(myIndex);
        auto rptr = &clusters[newrindex];
        rptr->clearCGgraphEntry(myIndex);
        cluster_ptr->clearCGgraphEntry(newlindex);
        cluster_ptr->clearCGgraphEntry(newrindex);
        // cluster<T,D>* new_heaviest_edge = cluster_ptr->max_weight_edge;

        lptr->insertCGentry(newrindex, new_weight, new_heaviest_edge);
        rptr->insertCGentry(newlindex, new_weight, new_heaviest_edge);
        
        cluster_ptr->vertex_count = 0;
        
        // old_cluster_ptr->vertex_count-=1; // TODO: only if the one being spliced out is myself.. is that not necessary? 
        
    }

    cluster_ptr = old_cluster_ptr;

    return;
}

template<typename T, typename D>
void sum_up_child_vertex_counts(cluster<T,D>* root)
{
    if(root == nullptr)
        return;
    auto& cluster_ptr = root;

    parlay::parallel_for(0, max_neighbours, [&] (T i) { 
        auto& child_ptr = cluster_ptr->children[i];
        sum_up_child_vertex_counts(child_ptr);
        return;
    });

    for(short i = 0; i < max_neighbours; i++)
    {
        auto& child_ptr = cluster_ptr->children[i];
        if(child_ptr)
            cluster_ptr->vertex_count+=child_ptr->vertex_count;
    }

}

template<typename T, typename D> // calculate vertex sizes
void write_vertices(cluster<T,D>* root, const T& offset, parlay::sequence<T>& all_vertices, parlay::sequence<cluster<T,D>>& clusters)
{
    auto& cluster_ptr = root; // so I can copy paste code from above
    if(cluster_ptr == nullptr)
        return;

    if(!(cluster_ptr->state & is_marked))
    {
        return;
    }

    std::array<T, max_neighbours> children_offsets;
    std::array<T, max_neighbours> children_sizes;
    T sum_till_now = 0;
    for(short i  = 0; i < max_neighbours; i++)
    {
        auto&child = cluster_ptr->children[i];
        children_sizes[i] = 0;
        children_offsets[i] = 0;
        if(!child)
            continue;
        children_sizes[i] = child->vertex_count;
    }
    children_offsets[0] = offset;

    for(short i = 1; i < max_neighbours; i++)
    {
        children_offsets[i] = children_offsets[i-1] + children_sizes[i-1];
    }
    T local_offset = children_offsets[max_neighbours-1] + children_sizes[max_neighbours-1];
    
    parlay::parallel_for(0, max_neighbours, [&] (T i) { 
        auto& child_ptr = cluster_ptr->children[i];
        write_vertices(child_ptr, children_offsets[i], all_vertices, clusters);
        return;
    });

    // std::string output;
    // output += "I, " + std::to_string(cluster_ptr->index) + " child of " + std::to_string((cluster_ptr->parent ? cluster_ptr->parent->index : -1)) + " have offset " + std::to_string(offset) + " and children:\t";
    // for (short i = 0; i < max_neighbours; i++) {
    //     auto& child = cluster_ptr->children[i];
    //     if (!child)
    //         output += "nl/" + std::to_string(children_sizes[i]) + "/" + std::to_string(children_offsets[i]) + "\t";
    //     else
    //         output += std::to_string(child->index) + "/" + std::to_string(children_sizes[i]) + "/" + std::to_string(children_offsets[i]) + "\t";
    // }
    // output += "\tMaking my write location be " + std::to_string(local_offset) + "\t";
    
    // sum of children sizes

    T sum_child_sizes = 0;
    for(auto& cs : children_sizes)
        sum_child_sizes+=cs;
    if(sum_child_sizes < cluster_ptr->vertex_count)
    {
        all_vertices[local_offset] = cluster_ptr->index;
        // output+= " and I wrote there\n";
    }
    // else
    //     output += "\n";
    // std::cout << output;

}

template<typename T, typename D>
void checkCGEdgesValid(parlay::sequence<T> cluster_indices, parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, cluster_indices.size(), [&] (T I){
        auto cluster_ptr = &clusters[cluster_indices[I]];
        for(short i = 0; i < max_neighbours; i++)
        {
            if(!cluster_ptr->CGEntryValid[i])
                continue;
            if(!cluster_ptr->CGGraphEdgePtr[i] || !(cluster_ptr->CGGraphEdgePtr[i]->data == cluster_ptr->CGGraphEntry[i].second))
            {
                cluster_ptr->print();
                if(cluster_ptr->CGGraphEdgePtr[i])
                    std::cout << "Entry is " << cluster_ptr->CGGraphEdgePtr[i]->index << "/" << cluster_ptr->CGGraphEdgePtr[i]->data << " vs " << cluster_ptr->CGGraphEntry[i].second << std::endl;
                assert(false && "CGEntry not valid");
            }
            
        }
    });

}


template<typename T, typename D>
std::pair<parlay::sequence<std::tuple<T,T,D>>, parlay::sequence<std::pair<T,T>>> createCompressedPathTree(parlay::sequence<cluster<T,D>>& clusters, const parlay::sequence<std::tuple<T,T,D>>& newEdges)
{
    parlay::sequence<std::tuple<T,T,D>> return_truples;

    parlay::sequence<cluster<T,D>*> relevant_null_clusters = parlay::sequence<cluster<T,D>*>(2*newEdges.size(), nullptr);
    
    static const int remove_mark_mask = ~(is_marked | is_marked_endpoint);

    // mark the endpoints and their parents as marked
    parlay::parallel_for(0, newEdges.size(), [&] (T i) {
        const auto& v = std::get<0>(newEdges[i]);
        const auto& w = std::get<1>(newEdges[i]);

        auto cluster_v_ptr = &clusters[v];
        auto cluster_w_ptr = &clusters[w];
        cluster_v_ptr->state &= remove_mark_mask;
        cluster_w_ptr->state &= remove_mark_mask;

        while(cluster_v_ptr != nullptr)
        {
            if(cluster_v_ptr->counter.fetch_add(1))
                break; // 
            cluster_v_ptr->state &= remove_mark_mask;
            for(auto& child : cluster_v_ptr->children)
                if(child)
                    child->state &= remove_mark_mask;

            cluster_v_ptr = cluster_v_ptr->parent;
        }

        while(cluster_w_ptr != nullptr)
        {
            if(cluster_w_ptr->counter.fetch_add(1))
                break; // 
            cluster_w_ptr->state &= remove_mark_mask;
            for(auto& child : cluster_w_ptr->children)
                if(child)
                    child->state &= remove_mark_mask;
            cluster_w_ptr = cluster_w_ptr->parent;
        }
    });

    // Mark them
    parlay::parallel_for(0, newEdges.size(), [&] (T i) {
        const auto& v = std::get<0>(newEdges[i]);
        const auto& w = std::get<1>(newEdges[i]);

        auto cluster_v_ptr = &clusters[v];
        auto cluster_w_ptr = &clusters[w];

        cluster_v_ptr->state |= is_marked_endpoint;
        cluster_w_ptr->state |= is_marked_endpoint;

        while(cluster_v_ptr != nullptr)
        {
            if(cluster_v_ptr->counter.fetch_add(-1) > 1)
                break; // we are not the last
            // std::cout << "Relevant cluster " << cluster_v_ptr->index << std::endl;
            cluster_v_ptr->state |= is_marked;
            if(cluster_v_ptr->parent == nullptr)
            {
                top_down_ct(cluster_v_ptr, clusters);
                relevant_null_clusters[i*2] = cluster_v_ptr;
            }
            cluster_v_ptr = cluster_v_ptr->parent;
        }

        while(cluster_w_ptr != nullptr)
        {
            if(cluster_w_ptr->counter.fetch_add(-1) > 1)
                break; // we are not the last
            cluster_w_ptr->state |= is_marked;
            if(cluster_w_ptr->parent == nullptr)
            {
                top_down_ct(cluster_w_ptr, clusters);
                relevant_null_clusters[i*2 + 1] = cluster_w_ptr;
            }
            cluster_w_ptr = cluster_w_ptr->parent;
        }
    });



    relevant_null_clusters = parlay::filter(relevant_null_clusters, [] (cluster<T,D>* ptr) {
        return ptr != nullptr;
    });

    parlay::parallel_for(0, relevant_null_clusters.size(), [&] (T i) {
        sum_up_child_vertex_counts(relevant_null_clusters[i]);
    });

    auto countpair = parlay::scan(
        parlay::tabulate(relevant_null_clusters.size(), [&] (T i) {
            return relevant_null_clusters[i]->vertex_count;
        })
    );

    parlay::sequence<T> per_forest_vertices = countpair.first;
    T totalcount = countpair.second;

    parlay::sequence<T> vertices = parlay::sequence<T>(totalcount, -1);
    parlay::parallel_for(0, relevant_null_clusters.size(), [&] (T i) {
        write_vertices(relevant_null_clusters[i], per_forest_vertices[i], vertices, clusters);
    });

    // verifyCompressPathTree(clusters, newEdges, vertices); // TODO REMOVE!!
    // checkCGEdgesValid(vertices, clusters); // TODO also remove!!

    parlay::sequence<T> edge_counts = parlay::tabulate(vertices.size(), [&] (T i) {
        auto cluster_ptr = &clusters[vertices[i]];
        assert(cluster_ptr->state & is_marked);
        return (T) cluster_ptr->getCGNumEdges();
    });

    auto edgeoffsetpair = parlay::scan(edge_counts);
    T& total_edge_count = edgeoffsetpair.second;
    parlay::sequence<T>& edge_offsets = edgeoffsetpair.first;

    
    return_truples = parlay::sequence<std::tuple<T,T,D>>(total_edge_count);
    auto return_endpoints = parlay::sequence<std::pair<T,T>>(total_edge_count);


    // write edges
    parlay::parallel_for(0, vertices.size(), [&] (T I) {
        auto cluster_ptr = &clusters[vertices[I]];
        auto local_offset = edge_offsets[I];
        for(short i = 0; i < max_neighbours; i++)
        {
            if(cluster_ptr->CGEntryValid[i] && (cluster_ptr->CGGraphEntry[i].first > cluster_ptr->index))
            {
                return_truples[local_offset] = std::tuple<T,T,D>(cluster_ptr->index, cluster_ptr->CGGraphEntry[i].first, cluster_ptr->CGGraphEntry[i].second);
                
                cluster<T,D>* heaviestEdgePtr = cluster_ptr->CGGraphEdgePtr[i];
                assert(heaviestEdgePtr);
                assert(heaviestEdgePtr->state & base_edge);
                cluster<T,D>* l;
                cluster<T,D>* r;
                heaviestEdgePtr->find_endpoints(l, r);
                assert(l && r);
                return_endpoints[local_offset] = {l->index, r->index};
                
                cluster_ptr->CGEntryValid[i] = false; // Now that we have this tuple, delete this edge
                local_offset++;
            }
        }
    });
    

    /*
        Just printing stuff below 
    */

    std::cout << "There are " << relevant_null_clusters.size() << " relevant forests" << std::endl;



    if(clusters.size() <= 100)
    {
        std::cout << "Counts [V/E]:" << std::endl;
        for(unsigned int i = 0; i < clusters.size(); i++)
        {
            std::cout << yellow << i << ": " << reset; 
            std::cout << yellow << clusters[i].vertex_count << " ";
            std::cout << red << "[ ";
            for(unsigned int cgi = 0; cgi < max_neighbours; cgi++)
            {
                if(clusters[i].CGEntryValid[cgi])
                    std::cout << "{ " << clusters[i].CGGraphEntry[cgi].first << " " << clusters[i].CGGraphEntry[cgi].second << "} ";
            }
            std::cout  << "] " << reset;
            std::cout << std::endl;
        }
        std::cout << "relevant null clusters: ";
        for(auto& ptr : relevant_null_clusters)
            std::cout << ptr->index << " ";
        std::cout << std::endl;
        std::cout << "relevant counts: ";
        for(auto& val : per_forest_vertices)
            std::cout << val << " ";
        std::cout << totalcount;
        std::cout << std::endl;
        std::cout << "relevant vertices final: ";
        for(auto& val : vertices)
            if(val != -1)
                std::cout << val << " ";
        std::cout << std::endl;
        std::cout << "edge counts: ";
        for(auto& edg : edge_counts)
            std::cout << edg << " ";
        std::cout << std::endl;

        std::cout << std::endl;
    }

    return {return_truples, return_endpoints};
}


template<typename T, typename D> //TODO make iterative  
std::pair<bool, D> vanillaQueryHelper(cluster<T,D>* ignore_ptr, cluster<T,D>* came_from_ptr, cluster<T,D>* destination)
{
    assert(destination);
    if(came_from_ptr == nullptr)
        return std::pair<bool, D>(false, 0);
    else if(came_from_ptr == destination)
        return std::pair<bool, D>(true, 0);
    D retval = 0;
    bool retbool = false;
    auto l1node = came_from_ptr->adjacency.get_head();
    assert(l1node);
    for(auto& l1ngbr : l1node->adjacents)
    {
        if(l1ngbr == nullptr)
            continue;
        cluster<T,D>* other_side = get_other_side(l1node, l1ngbr)->cluster_ptr;
        assert(other_side);
        D edgedata = l1ngbr->cluster_ptr->data;
        if(other_side == ignore_ptr)
            continue;
        // std::cout << "[" << came_from_ptr->index << "] Recursively checking out " << other_side->index << std::endl;
        auto retpair = vanillaQueryHelper(came_from_ptr, other_side, destination);
        if(retpair.first)
        {
            retbool = true;
            retval = retval > retpair.second ? retval : retpair.second;
            retval = retval > edgedata ? retval : edgedata;
        }
    }
    return std::pair<bool,D>(retbool, retval);
}


template<typename T, typename D>
D getVanillaQuery(parlay::sequence<cluster<T,D>>& clusters, T A, T B)
{
    auto retpair = vanillaQueryHelper(&clusters[A], &clusters[A], &clusters[B]);
    if(retpair.first)
        return retpair.second;
    return (D) 0;
}

template<typename T, typename D> //TODO make iterative
std::pair<bool, D> compressQueryHelper(parlay::sequence<cluster<T,D>>& clusters, cluster<T,D>* ignore_ptr, cluster<T,D>* came_from_ptr, cluster<T,D>* destination, const std::set<T>& all_vertices)
{
    if(came_from_ptr == destination)
        return std::pair<bool,D>(true, 0);
    D retval = 0;
    bool retbool = 0;
    for(unsigned int i = 0; i < max_neighbours; i++)
    {
        if(came_from_ptr->CGEntryValid[i])
        {
            T other_index = came_from_ptr->CGGraphEntry[i].first;
            if(!all_vertices.count(other_index))
            {
                std::cout << other_index << " not in set " << std::endl;
            }
            assert(all_vertices.count(other_index) && " vertex not in final set?");
            T data = came_from_ptr->CGGraphEntry[i].second;
            cluster<T,D>* other_side = &clusters[other_index];
            if(other_side == ignore_ptr)
                continue;
            // std::cout << "[" << came_from_ptr->index << "] Recursively checking out " << other_side->index << std::endl;
            auto retpair = compressQueryHelper(clusters, came_from_ptr, other_side, destination, all_vertices);
            if(retpair.first)
            {
                // std::cout << "[" << came_from_ptr->index << "] -> " << other_side->index << " was the correct path and had data " << data << std::endl;
                retbool = true;
                retval = retval > retpair.second ? retval : retpair.second;
                retval = retval > data ? retval : data;
            }
        }
    }
    return std::pair<bool,D>(retbool, retval);
}


template<typename T, typename D>
D getCompressQuery(parlay::sequence<cluster<T,D>>& clusters, T A, T B, const std::set<T>& all_vertices)
{
    
    auto retpair = compressQueryHelper(clusters, &clusters[A], &clusters[A], &clusters[B], all_vertices);

    if(retpair.first)
        return retpair.second;
    return (D) 0;
}

template<typename T, typename D> // TODO Haha this is O(n^3) I think but who cares, It is pretty trivial to replace these with par fors and they will not conflict since queries don't change state
// Plus, the helpers need to made iterative and parlay based (std allocator == death)
void verifyCompressPathTree(parlay::sequence<cluster<T,D>>& clusters, const parlay::sequence<std::tuple<T,T,D>>& newEdges, const parlay::sequence<T>& vertices)
{

    std::cout << "Checking to see if a random sample of endpoints match in compressed and actual tree" << std::endl;
    std::cout << "And to see if there are any vertices encountered that aren't in the vertex set. This might take a while." << std::endl; 

    if(!newEdges.size())
        return;
    std::set<T> all_vertices;
    for(auto& vert: vertices)
        all_vertices.insert(vert);


    static const int max_check_num = 10; // checks max_check_num times max_check_num
    
    parlay::parallel_for(0, max_check_num, [&] (T o) {
        o = rand() % newEdges.size();

        auto edge = newEdges[o];

        for(T j = 0; j < max_check_num; j++)
        {
            // T v = rand() % newEdges.size();
            // auto edge_ = newEdges[v];
            T i = j;
            if(i % 2) // mix of nearby and random
                i = rand() % newEdges.size();
            if((i+o) >= newEdges.size()) 
                break;
            auto edge_ = newEdges[i+o];
            T A = std::get<0>(edge);
            T B = std::get<1>(edge_);
            D actualPathQuery = getVanillaQuery(clusters, A, B);
            D compressPathQuery = getCompressQuery(clusters, A, B, all_vertices);
            // std::cout << "Checked out " << A << " -- " << B << blue << actualPathQuery << reset << " vs " << red << compressPathQuery << reset << std::endl;
            // assert(actualPathQuery == compressPathQuery && "Compress tree is invalid for queries");
            if(actualPathQuery != compressPathQuery) 
            {
                std::cout << "MISMATCH!" << std::endl;
                std::cout << "Actual path query from " << A << " -- " << B << ": " << red << actualPathQuery << reset << std::endl;
                std::cout << "Compress path query from " << A << " -- " << B << ": " << blue << compressPathQuery << reset << std::endl;
                if(clusters.size() <= 100)
                {
                    printTree(clusters);
                }
                assert("Actual path query must match compressed path query" && false);
            }
        }
    });

}



template<typename T, typename D>
void checkEdgePtrValid(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T I) {
        auto cluster_ptr = &clusters[I];
        auto heaviestEdgePtr = cluster_ptr->max_weight_edge;
        if(isBinary(cluster_ptr))
        {
            if(!heaviestEdgePtr || !(heaviestEdgePtr->data == cluster_ptr->data) || heaviestEdgePtr == cluster_ptr)
            {
                cluster_ptr->print();
                assert(false && "heaviest edge not maintained correctly");
            }
            
        }
    });
}

template<typename T, typename D>
void verifyInsertionObeysDegreeCap(parlay::sequence<cluster<T,D>>& clusters, parlay::sequence<std::tuple<T,T,D>>& add_edges)
{
    std::vector<std::atomic<short>> atomic_edgecounts = std::vector<std::atomic<short>>(clusters.size());

    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        short count = 0;
        for(auto& nd : clusters[i].adjacency.get_head()->adjacents)
            if(nd)
                count++;
        atomic_edgecounts[i]+=1;
    });

    parlay::parallel_for(0, add_edges.size(), [&] (T i) {
        T a = std::get<0>(add_edges[i]);
        T b = std::get<0>(add_edges[i]);
        short old_val_a = atomic_edgecounts[a].fetch_add(1);
        short old_val_b = atomic_edgecounts[b].fetch_add(1);
        assert(old_val_a <= (max_neighbours-1) && old_val_b <= (max_neighbours-1));

    });

    return;
}

/*
    Warning, changes cluster and MST in place
*/
template<typename T, typename D>
void incrementMST(parlay::sequence<cluster<T,D>>& clusters, const parlay::sequence<std::tuple<T,T,D>>& newEdges)
{
    auto compressed_tree_pair = createCompressedPathTree(clusters, newEdges);

    auto compressed_tree = compressed_tree_pair.first;
    auto compressed_tree_endpoints = compressed_tree_pair.second;

    if(clusters.size() <= 100)
    {
        std::cout << green << "Compressed tree: " << std::endl;
        for(const auto& edge : compressed_tree)
            std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
        std::cout << reset << std::endl;    
    }

    auto AllNewEdges = parlay::append(compressed_tree, newEdges);
    auto finalEdgeIndices = min_spanning_forest(AllNewEdges, clusters.size());

    parlay::sequence<bool> wasDeleted = parlay::sequence<bool>(AllNewEdges.size(), true);

    parlay::parallel_for(0, finalEdgeIndices.size(), [&] (T i) {
        wasDeleted[finalEdgeIndices[i]] = false;
    });


    if(clusters.size() <= 100)
    {
        std::cout << reset << "final edges: " << std::endl;
        for(long i = 0; i < AllNewEdges.size(); i++)
        {
            auto& edge = AllNewEdges[i];
            if(wasDeleted[i])
                std::cout << red;
            else
                std::cout << green;

            std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << " ";
            if(i < compressed_tree.size())
                std::cout << blue << compressed_tree_endpoints[i].first << " -- " << compressed_tree_endpoints[i].second << reset;
            std::cout << std::endl;
        }
        std::cout << reset << std::endl;    
    }

    // gather add edges via a filter which is O(n) work
    // TODO replace with a pack by breaking the boolean array into two
    auto add_edges = parlay::filter(parlay::tabulate(newEdges.size(), [&] (T i) {
        if(wasDeleted[i+compressed_tree.size()])
            return std::tuple<T,T,D>(0,0, (T) 0);
        return newEdges[i];
    }), [] (auto tup) {
        return std::get<0>(tup) != std::get<1>(tup);
    });

    auto delete_edges = parlay::filter(parlay::tabulate(compressed_tree.size(), [&] (T i) {
        if(!wasDeleted[i])
            return std::pair<T,T>(0,0);
        return compressed_tree_endpoints[i];
    }), [] (auto pr) {
        return pr.first != pr.second;
    });

    if(clusters.size() <= 100)
    {
        std::cout << "Final add edges " << std::endl;
        for(auto& edge : add_edges)
            std::cout << std::get<0>(edge)  << " -- " << std::get<1>(edge) << ": " << std::get<2>(edge) << std::endl;
        std::cout << std::endl;

        std::cout << "Final delete edges " << std::endl;
        for(auto& edge : delete_edges)
            std::cout << edge.first  << " -- " << edge.second << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Deleting " << delete_edges.size() << " edges " << std::endl;
    std::cout << "Adding " << add_edges.size() << " edges " << std::endl;

    parlay::sequence<std::tuple<T,T,D>> empty_add_edges;    
    parlay::sequence<std::pair<T,T>> empty_delete_edges;    

    batchInsertEdge(delete_edges, add_edges, clusters, (D) 0, [] (D a, D b) {return a > b ? a : b;});

    // std::cout << "Done deleting " << std::endl;
    // if(clusters.size() <= 100)
    //     printTree(clusters);

    // verifyInsertionObeysDegreeCap(clusters, add_edges); // TODO remove

    // batchInsertEdge(empty_delete_edges, add_edges, clusters, (D) 0, [] (D a, D b) {return a > b ? a : b;});


    return;
}

#endif