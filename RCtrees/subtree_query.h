#ifndef SUBTREE_QUERY_H
#define SUBTREE_QUERY_H

#include "RC.h"
#include "ternarizer.h"


template<typename T, typename D, typename assocfunc>
D subtree_query_helper(cluster<T, D>* root, cluster<T, D>* dir_giver, D defretval, assocfunc func)
{
    D ret_val = defretval;
    if(root == nullptr || dir_giver == nullptr || root == dir_giver)
        return defretval;
    if(PRINT_QUERY)
        std::cout << bold << bright_yellow << "Root: " << root->index << " dir: " << dir_giver->index << reset << std::endl;

    D lval = defretval;
    D rval = defretval;
    cluster<T,D>* dir_l = nullptr;
    cluster<T,D>* dir_r = nullptr;

    dir_giver->find_boundary_vertices(dir_l, lval, dir_r, rval, defretval);

    cluster<T,D>* root_l = nullptr;
    cluster<T,D>* root_r = nullptr;

    root->find_boundary_vertices(root_l, lval, root_r, lval, defretval);

    if(PRINT_QUERY)
    {
        std::cout << bold << yellow << "Root's boundary vertices ";
        if(root_l == nullptr)
            std::cout << "nl ";
        else
            std::cout << root_l->index << " ";
        if(root_r == nullptr)
            std::cout << "nl ";
        else
            std::cout << root_r->index << " ";
        std::cout << reset << std::endl;

        std::cout << bold << yellow << "dir_giver's boundary vertices ";
        if(dir_l == nullptr)
            std::cout << "nl ";
        else
            std::cout << dir_l->index << " ";
        if(dir_r == nullptr)
            std::cout << "nl ";
        else
            std::cout << dir_r->index << " ";
        std::cout << reset << std::endl;
    }

    if(root == dir_l || root == dir_r)
    {
        if(PRINT_QUERY)
            std::cout << bold << bright_cyan << "root already boundary vertex of dir giver" << reset << std::endl;

        auto child_with_dir = dir_giver;
        while(child_with_dir->parent != root)
            child_with_dir = child_with_dir->parent;

        auto& child_that_contains_u = child_with_dir; // rename

        cluster<T,D>* ign_l = nullptr;
        cluster<T,D>* ign_r = nullptr;
        child_that_contains_u->find_boundary_vertices(ign_l, lval, ign_r, rval, defretval);

        if(PRINT_QUERY)
            std::cout << bold << bright_cyan << "Child to ignore is " << child_with_dir->index << reset << std::endl;

        bool first = true;

        for(const auto& child : root->children)
        {
            if(child == nullptr)
                continue;

            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, defretval);
            if(PRINT_QUERY)
            {
                std::cout << bold << bright_cyan << "[" << root->index << "] child " << child->index << " boundary vertices ";
                if(child_l == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_l->index << " ";
                if(child_r == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_r->index << " ";
                std::cout << reset << std::endl;
            }

            if(child == child_with_dir)
                continue;
            
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "Adding value of " << child->index << "/" << child->data << reset << std::endl;
            
            if(first)
            {
                first = false;
                ret_val = child->data;
            }
            else
            {
                ret_val = func(ret_val, child->data);
            }
        }

        if(root_l != nullptr && root_l != ign_l && root_l != ign_r)
        {
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "recursive call to  " << root_l->index << " wrt " << root->index << reset << std::endl;
            auto child_ret_val = subtree_query_helper(root_l, root, defretval, func);
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "returning from recursive call to  " << root_l->index << " wrt " << root->index << ", adding " << child_ret_val << reset << std::endl;
            ret_val = func(ret_val, child_ret_val);
        }
        if(root_r != root_l && root_r != nullptr && root_r != ign_l && root_r != ign_r)
        {
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "recursive call to  " << root_r->index << " wrt " << root->index << reset << std::endl;
            auto child_ret_val = subtree_query_helper(root_r, root, defretval, func);
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "returning from recursive call to  " << root_r->index << " wrt " << root->index << ", adding " << child_ret_val << reset << std::endl;
            ret_val = func(ret_val, child_ret_val);    
        }
    }
    else
    {
        if(PRINT_QUERY)
            std::cout << bold << bright_green << "Should only go here once" << reset << std::endl;

        bool first = true;
        for(auto& child : root->children)
        {
            if(child == nullptr)
                continue;
            
            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, defretval);
            if(PRINT_QUERY)
            {
                std::cout << bold << bright_green << "[" << root->index << "] child " << child->index << " boundary vertices ";
                if(child_l == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_l->index << " ";
                if(child_r == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_r->index << " ";
                std::cout << reset << std::endl;
            }

            if((child->adjacency.get_head()->state & base_edge) || (child->first_contracted_node->next->state & binary_cluster))
            {
                if(child_l == dir_giver || child_r == dir_giver)
                    continue;
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                }
                else
                {
                    ret_val = func(ret_val, child->data);
                }
                if(child_l != root)
                {               
                    auto child_ret_val = subtree_query_helper(child_l, root, defretval, func);
                    ret_val = func(ret_val, child_ret_val);
                }
                if(child_l != child_r && child_r != root)
                {
                    auto child_ret_val = subtree_query_helper(child_r, root, defretval, func);
                    ret_val = func(ret_val, child_ret_val);
                }
            }
            else
            {
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                }
                else
                {
                    ret_val = func(ret_val, child->data);
                }
            }
            
        }
    }
    return ret_val;
}

template<typename T, typename D, typename assocfunc>
D subtree_query(T root_index, T dir_giver_index, parlay::sequence<cluster<T,D>>& clusters, ternarizer<T,D>& TR, D identity, assocfunc func)
{
    auto ternerized_edge = TR.translate_edge(root_index, dir_giver_index);
    return subtree_query_helper(&clusters[ternerized_edge.first], &clusters[ternerized_edge.second], identity, func);
}

template<typename T, typename D, typename assocfunc>
D manual_subtree_sum_helper(cluster<T,D>* root, cluster<T,D>* dir_giver, parlay::sequence<cluster<T,D>>& clusters, D identity, assocfunc func)
{
    D ret_val = identity;
    if (root == nullptr || dir_giver == nullptr || root == dir_giver)
        return identity;
    std::stack<std::pair<cluster<T, D>*, cluster<T, D>*>> stack;
    stack.push({root, dir_giver});

    while (!stack.empty()) {
        auto [current_root, parent] = stack.top();
        stack.pop();

        T r = current_root->index;

        for (auto& edge_ptr : current_root->adjacency.get_head()->adjacents) {
            if (edge_ptr == nullptr)
                continue;

            auto other_cluster = get_other_side(current_root->adjacency.get_head(), edge_ptr)->cluster_ptr;
            if (other_cluster == parent)
                continue;

            T other_index = other_cluster->index;
            ret_val = func(ret_val,get_edge(r, other_index, clusters)->cluster_ptr->data);

            // Push the next cluster to the stack with the current cluster as its parent
            stack.push({other_cluster, current_root});
        }
    }

    return ret_val;

}

template<typename T, typename D, typename assocfunc>
D manual_subtree_sum(T root_index, T dir_giver, parlay::sequence<cluster<T,D>>& clusters, ternarizer<T,D>& TR, D identity, assocfunc func)
{
    auto ternerized_edge = TR.translate_edge(root_index, dir_giver);
    return manual_subtree_sum_helper(&clusters[ternerized_edge.first], &clusters[ternerized_edge.second], clusters, identity, func);
}

template<typename T, typename D>
short find_child_index(cluster<T,D>* child, cluster<T,D>* ancestor) // find the child index of ancestor which contains the child cluster
{
    auto curr = child;
    while(curr->parent != ancestor)
    {
        curr = curr->parent;
        if(curr == nullptr)
            return -1;
    }
    for(short i = 0; i < max_neighbours; i++)
    {
        if(ancestor->children[i] == curr)
            return i;
    }
    return -1;
}

template<typename T, typename D>
D find_child_partial_sum(cluster<T,D>* child, cluster<T,D>* ancestor) // find the partial sum corresponding to this child i.e. the entire tree without this child
{
    // std::cout << "Got asked for " << ancestor->index << " wrt " << child->index << std::endl;
    return ancestor->partial_sums[find_child_index(child, ancestor)];
}


template<typename T, typename D>
void mark_entryway(cluster<T,D>* root, cluster<T,D>* dir_giver)
{

    if(root == nullptr || dir_giver == nullptr || root == dir_giver)
        return;
    // std::cout << root->index << " w.r.t. " << dir_giver->index << " attempted" << std::endl;

    D lval;
    D rval;
    cluster<T,D>* dir_l = nullptr;
    cluster<T,D>* dir_r = nullptr;
    dir_giver->find_boundary_vertices(dir_l, lval, dir_r, rval, lval);

    cluster<T,D>* root_l = nullptr;
    cluster<T,D>* root_r = nullptr;
    root->find_boundary_vertices(root_l, lval, root_r, lval, lval);

    if(root == dir_l || root == dir_r) // root is already a boundary vertex of the direction giver
    { // this means that the root is an ancestor of the direction giver
        short root_as_child_index = find_child_index(dir_giver, root);
        // was I the first one to reach this
        auto& relevant_atomic_char = root->partial_sum_complete[root_as_child_index];
        char expected = 0;
        if(relevant_atomic_char.load() != expected)
            return;
        char desired = 1;
        if (!relevant_atomic_char.compare_exchange_strong(expected, desired)) // If not zero, return
        {
            // std::cout << "Stopped" << std::endl;
            return;
        }
        auto child_that_contains_u = root->children[root_as_child_index];
        
        cluster<T,D>* ign_l = nullptr;
        cluster<T,D>* ign_r = nullptr;
        child_that_contains_u->find_boundary_vertices(ign_l, lval, ign_r, rval, lval);

        // if(root_l != nullptr && root_l != ign_l && root_l != ign_r) 
        // {
            mark_entryway(root_l, root);
        // }
        // if(root_r != root_l && root_r != nullptr && root_r != ign_l && root_r != ign_r)
        // {
            mark_entryway(root_r, root);
        // }
    }
    else
    {
        for(auto& child : root->children)
        {
            if(child == nullptr)
                continue;
            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, lval);

            if((child->adjacency.get_head()->state & base_edge) || (child->first_contracted_node->next->state & binary_cluster))
            {
                if(child_l == dir_giver || child_r == dir_giver)
                    continue;
                
                if(child_l != root)
                {               
                    mark_entryway(child_l, root);
                }
                if(child_l != child_r && child_r != root)
                {
                    mark_entryway(child_r, root);
                }
            }
        }
    }
}

template<typename T, typename D, typename ds>
parlay::sequence<cluster<T,D>*> get_roots(ds& clusters)
{
    auto potential_roots = parlay::filter(parlay::tabulate(clusters.size(), [&] (T i) {
        cluster<T,D>* curr = clusters[i];
        cluster<T,D>* ret_null = nullptr;
        do
        {
        T old_val = curr->counter.fetch_add(1);
        if(old_val)
            return ret_null;
        // for(auto& val : curr->partial_sum_complete)
        //     val = 0;
        if(curr->parent == nullptr)
            return curr;
        curr = curr->parent;
        }while(curr != nullptr);
        return ret_null;
    }), [] (auto ptr) {return ptr != nullptr;});
    return std::move(potential_roots);
}

template<typename T, typename D, typename ds>
void reset_roots(ds& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        cluster<T,D>* curr = clusters[i];
        do
        {
        T old_val = curr->counter.fetch_add(-1);
        if(old_val != 1) 
            return;
        for(auto& val : curr->partial_sum_complete)
            val = 0;
        if(curr->parent == nullptr)
            return;
        curr = curr->parent;
        }while(curr != nullptr);
    });
    return;
}

template<typename T, typename D>
D get_cached_query(cluster<T,D>* root, cluster<T,D>* dir_giver)
{
    short child_index = find_child_index(dir_giver, root);
    // assert(child_index != -1); // dir_giver should be a child of root
    // assert(root->partial_sum_complete[child_index] == 2); // should be cached already
    return root->partial_sums[child_index];
}

template<typename T, typename D, typename assocfunc>
D use_cached_query(cluster<T,D>* root, cluster<T,D>* dir_giver, D identity, assocfunc func)
{
    if(root == nullptr || dir_giver == nullptr || root == dir_giver)
        return identity;
    D ret_val = identity;
    D lval = identity;
    D rval = identity;
    cluster<T,D>* dir_l = nullptr;
    cluster<T,D>* dir_r = nullptr;
    dir_giver->find_boundary_vertices(dir_l, lval, dir_r, rval, identity);

    cluster<T,D>* root_l = nullptr;
    cluster<T,D>* root_r = nullptr;
    root->find_boundary_vertices(root_l, lval, root_r, lval, identity);

    if(root == dir_l || root == dir_r)
    {
        auto child_with_dir = dir_giver;
        while(child_with_dir->parent != root)
            child_with_dir = child_with_dir->parent;

        auto& child_that_contains_u = child_with_dir; // rename

        cluster<T,D>* ign_l = nullptr;
        cluster<T,D>* ign_r = nullptr;
        child_that_contains_u->find_boundary_vertices(ign_l, lval, ign_r, rval, identity);

        bool first = true;

        for(const auto& child : root->children)
        {
            if(child == nullptr)
                continue;

            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, identity);

            if(child == child_with_dir)
                continue;
                            
            if(first)
            {
                first = false;
                ret_val = child->data;
            }
            else
            {
                ret_val = func(ret_val, child->data);
            }
        }

        if(root_l != nullptr && root_l != ign_l && root_l != ign_r)
        {
            auto child_ret_val = get_cached_query(root_l, root);
            ret_val = func(ret_val, child_ret_val);
        }
        if(root_r != root_l && root_r != nullptr && root_r != ign_l && root_r != ign_r)
        {
            auto child_ret_val = get_cached_query(root_r, root);
            ret_val = func(ret_val, child_ret_val);    
        }
    }
    else
    {
        // std::cout << "Went here " << root->index << " -- " << dir_giver->index << std::endl;
        bool first = true;
        for(auto& child : root->children)
        {
            if(child == nullptr)
                continue;
            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, identity);

            if((child->adjacency.get_head()->state & base_edge) || (child->first_contracted_node->next->state & binary_cluster))
            {
                if(child_l == dir_giver || child_r == dir_giver)
                    continue;
                // if(child_l)
                //     std::cout << "Went here for " << child_l->index;
                // if(child_r && child_r != child_l)
                //     std::cout << " and " << child_r->index;
                // std::cout << std::endl;
                // std::cout << "Note that the dir_giver was " << dir_giver->index  << std::endl;
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                }
                else
                {
                    ret_val = func(ret_val, child->data);
                }
                if(child_l != root)
                {               
                    // auto child_ret_val = subtree_query_helper(child_l, root, identity, func);
                    auto child_ret_val = get_cached_query(child_l, root);
                    // std::cout << "I asked for a partial sum of " << root->index << " from " << child_l->index << " which was " << child_ret_val << std::endl;
                    ret_val = func(ret_val, child_ret_val);
                }
                if(child_l != child_r && child_r != root)
                {
                    // auto child_ret_val = subtree_query_helper(child_r, root, identity, func);
                    auto child_ret_val = get_cached_query(child_r, root);
                    // std::cout << "And I asked for a partial sum of " << root->index << " from " << child_r->index << " which was " << child_ret_val << std::endl;
                    ret_val = func(ret_val, child_ret_val);
                }
            }
            else
            {
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                    // std::cout << "Adding value of child/value " << child->index << "/" << child->data << std::endl;
                }
                else
                {
                    // std::cout << "Adding value of child/value " << child->index << "/" << child->data << std::endl;
                    ret_val = func(ret_val, child->data);
                }
            }
            
        }
    }
    return ret_val;
}

template<typename T, typename D, typename assocfunc>
void build_partial_sums(cluster<T,D>* root, D identity, assocfunc func) // top-down
{
    root->counter = 0; // reset counter
    parlay::parallel_for(0, max_neighbours, [&] (T i) {
        if(root->children[i] == nullptr)
            return;
        if(root->children[i]->counter == 0)
            return;
        // if(root->partial_sum_complete[i] != 1) // it has either been done before it isn't necessary
        //     return; // not functionally needed but speedup I guess
        // char expected = 1;
        // char desired = 2;
        // if (!root->partial_sum_complete[i].compare_exchange_strong(expected, desired)) // If not zero, return
        // {
        //     // someone did this before me
        //     return;
        // }

        root->partial_sums[i] = use_cached_query(root, root->children[i], identity, func);

        build_partial_sums(root->children[i], identity, func);
        return;
    });
} //TODO* Q* when does this code ever finalize the root's value from the partial sums? 

template<typename T, typename D, typename assocfunc, typename ds>
parlay::sequence<D> do_constant_time_batch_query(ds& indices, parlay::sequence<cluster<T,D>>& clusters, D identity, assocfunc func)
{
    auto ret_vals = parlay::tabulate(indices.size(), [&] (T II) {
        auto root = &clusters[indices[II].first];
        auto dir_giver = &clusters[indices[II].second];
        D ret_val = identity;
        
        D lval = identity;
        D rval = identity;
        cluster<T,D>* dir_l = nullptr;
        cluster<T,D>* dir_r = nullptr;
        dir_giver->find_boundary_vertices(dir_l, lval, dir_r, rval, identity);

        cluster<T,D>* root_l = nullptr;
        cluster<T,D>* root_r = nullptr;
        root->find_boundary_vertices(root_l, lval, root_r, lval, identity);

        if(root == dir_l || root == dir_r)
        {
            auto child_with_dir = dir_giver;
            while(child_with_dir->parent != root)
                child_with_dir = child_with_dir->parent;

            auto& child_that_contains_u = child_with_dir; // rename

            cluster<T,D>* ign_l = nullptr;
            cluster<T,D>* ign_r = nullptr;
            child_that_contains_u->find_boundary_vertices(ign_l, lval, ign_r, rval, identity);

            bool first = true;

            for(const auto& child : root->children)
            {
                if(child == nullptr)
                    continue;

                cluster<T,D>* child_l = nullptr;
                cluster<T,D>* child_r = nullptr;
                if(child->adjacency.get_head()->state & base_edge)
                    child->find_endpoints(child_l, child_r);
                else
                    child->find_boundary_vertices(child_l, lval, child_r, rval, identity);

                if(child == child_with_dir)
                    continue;
                                
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                }
                else
                {
                    ret_val = func(ret_val, child->data);
                }
            }

            if(root_l != nullptr && root_l != ign_l && root_l != ign_r)
            {
                auto child_ret_val = get_cached_query(root_l, root);
                ret_val = func(ret_val, child_ret_val);
            }
            if(root_r != root_l && root_r != nullptr && root_r != ign_l && root_r != ign_r)
            {
                auto child_ret_val = get_cached_query(root_r, root);
                ret_val = func(ret_val, child_ret_val);    
            }
        }
        else
        {
            // std::cout << "Went here " << root->index << " -- " << dir_giver->index << std::endl;
            bool first = true;
            for(auto& child : root->children)
            {
                if(child == nullptr)
                    continue;
                cluster<T,D>* child_l = nullptr;
                cluster<T,D>* child_r = nullptr;
                if(child->adjacency.get_head()->state & base_edge)
                    child->find_endpoints(child_l, child_r);
                else
                    child->find_boundary_vertices(child_l, lval, child_r, rval, identity);

                if((child->adjacency.get_head()->state & base_edge) || (child->first_contracted_node->next->state & binary_cluster))
                {
                    if(child_l == dir_giver || child_r == dir_giver)
                        continue;
                    // if(child_l)
                    //     std::cout << "Went here for " << child_l->index;
                    // if(child_r && child_r != child_l)
                    //     std::cout << " and " << child_r->index;
                    // std::cout << std::endl;
                    // std::cout << "Note that the dir_giver was " << dir_giver->index  << std::endl;
                    if(first)
                    {
                        first = false;
                        ret_val = child->data;
                    }
                    else
                    {
                        ret_val = func(ret_val, child->data);
                    }
                    if(child_l != root)
                    {               
                        // auto child_ret_val = subtree_query_helper(child_l, root, identity, func);
                        auto child_ret_val = get_cached_query(child_l, root);
                        // std::cout << "I asked for a partial sum of " << root->index << " from " << child_l->index << " which was " << child_ret_val << std::endl;
                        ret_val = func(ret_val, child_ret_val);
                    }
                    if(child_l != child_r && child_r != root)
                    {
                        // auto child_ret_val = subtree_query_helper(child_r, root, identity, func);
                        auto child_ret_val = get_cached_query(child_r, root);
                        // std::cout << "And I asked for a partial sum of " << root->index << " from " << child_r->index << " which was " << child_ret_val << std::endl;
                        ret_val = func(ret_val, child_ret_val);
                    }
                }
                else
                {
                    if(first)
                    {
                        first = false;
                        ret_val = child->data;
                        // std::cout << "Adding value of child/value " << child->index << "/" << child->data << std::endl;
                    }
                    else
                    {
                        // std::cout << "Adding value of child/value " << child->index << "/" << child->data << std::endl;
                        ret_val = func(ret_val, child->data);
                    }
                }
                
            }
        }
        return ret_val;
    });

    return std::move(ret_vals);
}


template<typename T, typename D, typename assocfunc, typename ds>
parlay::sequence<D> batched_subtree_query_helper(ds& indices, parlay::sequence<cluster<T,D>>& clusters, D identity, assocfunc func)
{
    parlay::sequence<D> ret_seq;

    auto cluster_ptrs = parlay::delayed_tabulate(indices.size(), [&] (T i) {
        return &clusters[indices[i].first];
    });

    auto roots = get_roots<T,D>(cluster_ptrs);

    // parlay::parallel_for(0, indices.size(), [&] (T i) {
    //     const auto& index_pair = indices[i];
    //     cluster<T,D>* root = &clusters[index_pair.first];
    //     cluster<T,D>* dir_giver = &clusters[index_pair.second];
    //     mark_entryway(root, dir_giver);
    //     return;
    // });

    // std::cout << "There are " << roots.size() << " unique roots " << std::endl;

    // for(auto& root : roots)
    //     std::cout << root->index << " ";
    // std::cout << std::endl;

    parlay::parallel_for(0, roots.size(), [&] (T i){
        build_partial_sums(roots[i], identity, func);
    });

    ret_seq = do_constant_time_batch_query(indices, clusters, identity, func);

    // printTree(clusters);

    // reset_roots<T,D>(cluster_ptrs);

    return std::move(ret_seq);   
}

template<typename T, typename D, typename assocfunc>
parlay::sequence<D> batched_subtree_query(parlay::sequence<std::pair<T,T>>& indices, parlay::sequence<cluster<T,D>>& clusters, ternarizer<T,D>& TR, D identity, assocfunc func)
{
    auto ternerized_indices = parlay::delayed_tabulate(indices.size(), [&] (T i){
        T v = indices[i].first;
        T w = indices[i].second;
        return TR.translate_edge(v,w);
    }); 
    return batched_subtree_query_helper(ternerized_indices, clusters, identity, func);
}


#endif