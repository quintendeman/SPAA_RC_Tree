#include "cluster.h"
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"

// Use chatgpt to convert from recursive to iterative
//CITE*? (TOD2* do we cite chatgpt here?)
template<typename T, typename D>
D manual_subtree_sum(cluster<T, D>* root, cluster<T, D>* dir_giver, parlay::sequence<cluster<T,D>>& clusters) {
    D ret_val = 0.0;
    if (root == nullptr || dir_giver == nullptr || root == dir_giver)
        return 0.0;

    // Stack to store pairs of clusters and their corresponding parent
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
            ret_val += get_edge(r, other_index, clusters)->cluster_ptr->data;

            // Push the next cluster to the stack with the current cluster as its parent
            stack.push({other_cluster, current_root});
        }
    }

    return ret_val;
}



template<typename T, typename D, typename assocfunc>
D subtree_query(cluster<T, D>* root, cluster<T, D>* dir_giver, D defretval, assocfunc func)
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
            auto child_ret_val = subtree_query(root_l, root, defretval, func);
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "returning from recursive call to  " << root_l->index << " wrt " << root->index << ", adding " << child_ret_val << reset << std::endl;
            ret_val = func(ret_val, child_ret_val);
        }
        if(root_r != root_l && root_r != nullptr && root_r != ign_l && root_r != ign_r)
        {
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "recursive call to  " << root_r->index << " wrt " << root->index << reset << std::endl;
            auto child_ret_val = subtree_query(root_r, root, defretval, func);
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
                    auto child_ret_val = subtree_query(child_l, root, defretval, func);
                    ret_val = func(ret_val, child_ret_val);
                }
                if(child_l != child_r && child_r != root)
                {
                    auto child_ret_val = subtree_query(child_r, root, defretval, func);
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

//TOD2* Add batch subtree query
