//
// Created by Bojie Shen on 19/8/20.
//

#ifndef ROUTINGKIT_MODIFIED_DIJKSTRA_H
#define ROUTINGKIT_MODIFIED_DIJKSTRA_H
#pragma once

#include <searchnode.h>
#include <vector>
#include <queue>
#include <set>
#include <cpool.h>
namespace RoutingKit{
    template<typename T, typename Compare = std::greater<T> >
    struct PointerComp
    {
        bool operator()(const T* x,
                        const T* y) const
        {
            return Compare()(*x, *y);
        }
    };


    class Modified_dijkstra {
        typedef std::priority_queue <SearchNodePtr, std::vector<SearchNodePtr>,
        PointerComp<SearchNode>> pq;

    public:
        void gen_initial_nodes(unsigned source_node)
        {
            dist[source_node] = 0;
            unsigned fm_index = 0;
            auto source = new (node_pool->allocate()) SearchNode {nullptr, source_node, 0, true, 0,(*ranking_mapper)[source_node]};
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a], fm_index, true,
                                                                  cur_weight,(*ranking_mapper)[(*head)[a]]};
                nptr->parent = source;
                nptr->is_forward = (*ranking_mapper)[source_node] < (*ranking_mapper)[(*head)[a]];
                open_list.push(nptr);
                fm_index++;
            }
        }

        void run(unsigned source_node){
            init_search();
            gen_initial_nodes(source_node);
            while(!open_list.empty()){
                SearchNodePtr cur_node = open_list.top(); open_list.pop();
                unsigned cur_root = cur_node->root;
                unsigned cur_distance = cur_node->g;
                // record node final nodes.

                if(cur_distance <= dist[cur_root] ) {
                    bool first_move_exist = false;
                    for (SearchNodePtr final:final_nodes[cur_root]) {
                        if (final->first_move == cur_node->first_move) {
                            first_move_exist = true;
                            break;
                        }
                    }
                    if (!first_move_exist) {
                        final_nodes[cur_root].push_back(cur_node);
                        dist[cur_root] = cur_distance;
                        for (unsigned a = (*first_out)[cur_root]; a < (*first_out)[cur_root + 1]; ++a) {
                            unsigned cur_weight = (*weight)[a];
                            unsigned distance_to_next = cur_distance + cur_weight;
                            if (distance_to_next <= dist[(*head)[a]]) {
                                if (cur_node->is_forward) {
                                    // generate all possible nodes;
                                    auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a],
                                                                                      cur_node->first_move, true,
                                                                                      distance_to_next,
                                                                                      (*ranking_mapper)[(*head)[a]]};
                                    nptr->parent = cur_node;
                                    nptr->is_forward = (*ranking_mapper)[cur_root] < (*ranking_mapper)[(*head)[a]];
                                    open_list.push(nptr);
                                } else {
                                    if ((*ranking_mapper)[cur_root] > (*ranking_mapper)[(*head)[a]]) {
                                        auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a],
                                                                                          cur_node->first_move, false,
                                                                                          distance_to_next,
                                                                                          (*ranking_mapper)[(*head)[a]]};
                                        nptr->parent = cur_node;
                                        open_list.push(nptr);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        void run_dijsktra(unsigned source_node){
            init_search();
            dist[source_node] = 0;
            unsigned fm_index = 0;
            auto source = new (node_pool->allocate()) SearchNode {nullptr, 0, 0, true, 0,(*ranking_mapper)[source_node]};
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight =  (*weight)[a];

                auto nptr = new (node_pool->allocate()) SearchNode {nullptr, (*head)[a]+1,
                                                                    fm_index, true,cur_weight,(*ranking_mapper)[(*head)[a]]};
                nptr->parent = source;
                open_list.push(nptr);
                fm_index++;
            }
            while(!open_list.empty()){
                SearchNodePtr cur_node = open_list.top(); open_list.pop();
                unsigned cur_root = cur_node->root-1;
                unsigned cur_distance = cur_node->g;
                // record node final nodes.
                if(cur_distance <= dist[cur_root] ){
                    bool first_move_exist = false;
                    for(SearchNodePtr final:final_nodes[cur_root]){
                        if(final->first_move == cur_node->first_move){
                            first_move_exist = true;
                            break;
                        }
                    }
                    if(!first_move_exist) {
                        final_nodes[cur_root].push_back(cur_node);
                        dist[cur_root] = cur_distance;
                        for (unsigned a = (*first_out)[cur_root]; a < (*first_out)[cur_root + 1]; ++a) {

                            unsigned cur_weight = (*weight)[a];

                            unsigned distance_to_next = cur_distance + cur_weight;

                            if (distance_to_next <= dist[(*head)[a]]) {
                                // generate all possible nodes;
                                auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a] + 1,
                                                                                  cur_node->first_move, true,
                                                                                  distance_to_next,(*ranking_mapper)[(*head)[a]]};
//                                if((*head)[a] == cur_root || cur_weight == 0){
//                                    std::cout<<"graph error"<< std::endl;
//                                }
                                nptr->parent = cur_node;
                                open_list.push(nptr);
                            }


                        }
                    }
                }

//                if(dist[cur_node->root-1] != std::numeric_limits<int>::max()){
//                    if(cur_distance < dist[cur_node->root-1] ){
//                        std::cout<<"This should not happed"<<std::endl;
//                    }
//                }

                // expand nodes


            }
        }

        vector<unsigned> run_dijsktra_to_get_distance(unsigned source_node){
            vector<unsigned> distance_list  =  vector<unsigned>(dist.size());
            init_search();
            dist[source_node] = 0;
            unsigned fm_index = 0;
            auto source = new (node_pool->allocate()) SearchNode {nullptr, 0, 0, true, 0,(*ranking_mapper)[source_node]};
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight =  (*weight)[a];

                auto nptr = new (node_pool->allocate()) SearchNode {nullptr, (*head)[a]+1,
                                                                    fm_index, true,cur_weight,(*ranking_mapper)[(*head)[a]]};
                nptr->parent = source;
                open_list.push(nptr);
                fm_index++;
            }
            while(!open_list.empty()){
                SearchNodePtr cur_node = open_list.top(); open_list.pop();
                unsigned cur_root = cur_node->root-1;
                unsigned cur_distance = cur_node->g;
                // record node final nodes.
                if(cur_distance <= dist[cur_root] ){
                    bool first_move_exist = false;
                    for(SearchNodePtr final:final_nodes[cur_root]){
                        if(final->first_move == cur_node->first_move){
                            first_move_exist = true;
                            break;
                        }
                    }
                    if(!first_move_exist) {
                        final_nodes[cur_root].push_back(cur_node);
                        dist[cur_root] = cur_distance;
                        for (unsigned a = (*first_out)[cur_root]; a < (*first_out)[cur_root + 1]; ++a) {

                            unsigned cur_weight = (*weight)[a];

                            unsigned distance_to_next = cur_distance + cur_weight;

                            if (distance_to_next <= dist[(*head)[a]]) {
                                // generate all possible nodes;
                                auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a] + 1,
                                                                                  cur_node->first_move, true,
                                                                                  distance_to_next,(*ranking_mapper)[(*head)[a]]};
//                                if((*head)[a] == cur_root || cur_weight == 0){
//                                    std::cout<<"graph error"<< std::endl;
//                                }
                                nptr->parent = cur_node;
                                open_list.push(nptr);
                            }


                        }
                    }
                }

//                if(dist[cur_node->root-1] != std::numeric_limits<int>::max()){
//                    if(cur_distance < dist[cur_node->root-1] ){
//                        std::cout<<"This should not happed"<<std::endl;
//                    }
//                }

                // expand nodes


            }
            for(int i = 0; i < final_nodes.size(); i ++){
                if(i==source_node){
                    distance_list[i] = 0;
                }else{
                    distance_list[i] = final_nodes[i][0]->g;
                }

            }

            return distance_list;
        }
        vector<unsigned> run_down_dijsktra(unsigned source_node, int& number_of_reached_nodes,vector<unsigned>distance_list){
            init_search();
            dist[source_node] = 0;
            unsigned fm_index = 0;
            auto source = new (node_pool->allocate()) SearchNode {nullptr, 0, 0, true, 0,(*ranking_mapper)[source_node]};
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){

                if((*ranking_mapper)[source_node]>(*ranking_mapper)[(*head)[a]]) {
                    unsigned cur_weight = (*weight)[a];

                    auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a] + 1,
                                                                      fm_index, true, cur_weight,
                                                                      (*ranking_mapper)[(*head)[a]]};
                    nptr->parent = source;
                    open_list.push(nptr);
                }
                fm_index++;
            }
            while(!open_list.empty()){
                SearchNodePtr cur_node = open_list.top(); open_list.pop();
                unsigned cur_root = cur_node->root-1;
                unsigned cur_distance = cur_node->g;
                // record node final nodes.
                if(cur_distance <= dist[cur_root] ){
//                    bool first_move_exist = false;
//                    for(SearchNodePtr final:final_nodes[cur_root]){
//                        if(final->first_move == cur_node->first_move){
//                            first_move_exist = true;
//                            break;
//                        }
//                    }
//                    if(!first_move_exist) {
                        final_nodes[cur_root].push_back(cur_node);
                        dist[cur_root] = cur_distance;
                        for (unsigned a = (*first_out)[cur_root]; a < (*first_out)[cur_root + 1]; ++a) {

                            unsigned cur_weight = (*weight)[a];

                            unsigned distance_to_next = cur_distance + cur_weight;

                            if (distance_to_next <= dist[(*head)[a]]) {
                                // generate all possible nodes;
                                if((*ranking_mapper)[cur_root]>(*ranking_mapper)[(*head)[a]]) {
                                    auto nptr = new(node_pool->allocate()) SearchNode{nullptr, (*head)[a] + 1,
                                                                                      cur_node->first_move, true,
                                                                                      distance_to_next,
                                                                                      (*ranking_mapper)[(*head)[a]]};
//                                if((*head)[a] == cur_root || cur_weight == 0){
//                                    std::cout<<"graph error"<< std::endl;
//                                }
                                    nptr->parent = cur_node;
                                    open_list.push(nptr);
                                }
                            }


                        }
//                    }
                }

//                if(dist[cur_node->root-1] != std::numeric_limits<int>::max()){
//                    if(cur_distance < dist[cur_node->root-1] ){
//                        std::cout<<"This should not happed"<<std::endl;
//                    }
//                }

                // expand nodes


            }
            vector<unsigned> reachable_nodes;
            unsigned index =0;
            for(auto nodes :final_nodes){
                if(index != source_node) {
                    if (!nodes.empty()) {
                        for(SearchNodePtr node : nodes){
                            if(node->g == distance_list[index]){
                                reachable_nodes.push_back(index);
                                number_of_reached_nodes++;
                                break;
                            }
                        }
                    }
                }
                index++;
            }
            return reachable_nodes;
        }
        void init_search()
        {
            assert(node_pool);
            node_pool->reclaim();
            open_list = pq();
            std::fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
            final_nodes.clear();
            final_nodes.resize(first_out->size()-1);
        }


        const std::vector<std::set<unsigned short>>& get_first_move_table(unsigned source, int& number_of_down_symbols){
            first_move_table.clear();
            first_move_table.resize(first_out->size()-1);
            unsigned a=(*first_out)[source];
            for(int i = 0; i < final_nodes.size(); i++){
                if(source == i ){
                    continue;
                }
                for(SearchNodePtr nptr: final_nodes[i]){
                    if((*ranking_mapper)[source] < (*ranking_mapper)[(*head)[a+nptr->first_move]]){
                        //only add upward moves;
                        first_move_table[i].insert(nptr->first_move);
                    }else{
                        if( first_move_table[i].find(0XFF) ==first_move_table[i].end() ){
                            first_move_table[i].insert(0xFF);
                            number_of_down_symbols++;
                        }
                    }

//                    first_move_table[i].insert(nptr->first_move+2);
                }
                if( first_move_table[i].empty() ){
                    first_move_table[i].insert(0xFF);
                }
            }
            return first_move_table;
        }

        std::vector<unsigned > get_down_symbol_only(unsigned source, int& number_of_down_symbols){
            first_move_table.clear();
            first_move_table.resize(first_out->size()-1);
            std::vector<unsigned> reachable;
            unsigned a=(*first_out)[source];
            for(int i = 0; i < final_nodes.size(); i++){
                if(source == i ){
                    continue;
                }
                for(SearchNodePtr nptr: final_nodes[i]){
                    if((*ranking_mapper)[source] < (*ranking_mapper)[(*head)[a+nptr->first_move]]){
                        //only add upward moves;
                        first_move_table[i].insert(nptr->first_move);
                    }else{
                        if( first_move_table[i].find(0XFF) ==first_move_table[i].end() ){
                            first_move_table[i].insert(0xFF);
                            number_of_down_symbols++;
                            reachable.push_back(i);
                        }
                    }
                }
            }
            return reachable;
        }


        std::vector<std::vector<unsigned>> get_peek_nodes_only(unsigned source){
            first_move_table.clear();
            first_move_table.resize(first_out->size()-1);
            std::vector<std::vector<unsigned>> peek_nodes = std::vector<std::vector<unsigned>>(first_out->size()-1) ;
            unsigned a=(*first_out)[source];
            for(int i = 0; i < final_nodes.size(); i++){
                if(source == i ){
                    continue;
                }
                for(SearchNodePtr nptr: final_nodes[i]){
                    unsigned best_node = nptr->root;
                    unsigned best_peek = nptr->rank;
                    while(nptr->parent != nullptr){
                        if(best_peek<nptr->rank){
                            best_node = nptr->root;
                            best_peek = nptr->rank;
                        }
                        nptr = nptr->parent;
                    }
                    peek_nodes[i].push_back(best_node-1);
                }
            }
            return  peek_nodes;
        }

        const std::vector<std::set<unsigned short>>& get_first_move_table_with_down_symbol(unsigned source){
            first_move_table.clear();
            first_move_table.resize(first_out->size()-1);
            unsigned a=(*first_out)[source];
            for(int i = 0; i < final_nodes.size(); i++){
                if(source == i ){
                    continue;
                }
                for(SearchNodePtr nptr: final_nodes[i]){
                    first_move_table[i].insert(nptr->first_move);
                }
                if( first_move_table[i].empty() ){
                    first_move_table[i].insert(0xFF);
                }
            }
            return first_move_table;
        }


        const std::vector<std::set<unsigned short>>& get_first_direction_table(unsigned source){
            first_move_table.clear();
            first_move_table.resize(first_out->size()-1);
            unsigned a=(*first_out)[source];
            for(int i = 0; i < final_nodes.size(); i++){
                if(source == i ){
                    continue;
                }
                for(SearchNodePtr nptr: final_nodes[i]){
                    if((*ranking_mapper)[source] < (*ranking_mapper)[(*head)[a+nptr->first_move]]){
                        //only add upward moves;
                        first_move_table[i].insert(2);
                    }else{
                        first_move_table[i].clear();
                        first_move_table[i].insert(3);
                        break;
                    }

//                    first_move_table[i].insert(nptr->first_move+2);
                }
                if( first_move_table[i].empty() ){
                    first_move_table[i].insert(1);
                }
            }
            return first_move_table;
        }


        Modified_dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail,
                const std::vector<unsigned>&head,
                const std::vector<unsigned>&weight,const std::vector<unsigned>&rank ):
//                tentative_distance(first_out.size()-1),
//                predecessor_arc(first_out.size()-1),
//                was_popped(first_out.size()-1),
//                queue(first_out.size()-1),
            dist(first_out.size()-1),
            final_nodes(first_out.size()-1),
            first_move_table(first_out.size()-1),
            first_out(&first_out),
            tail(&tail),
            head(&head),
            weight(&weight),
            ranking_mapper(&rank){
            node_pool = new warthog::mem::cpool(sizeof(SearchNode));
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());


        }

    private:

        pq open_list;

        std::vector<unsigned>dist;
        std::vector<std::vector<SearchNodePtr>>final_nodes;
        std::vector<std::set<unsigned short>>first_move_table;
        warthog::mem::cpool* node_pool;
        const std::vector<unsigned>*first_out;
        const std::vector<unsigned>*tail;
        const std::vector<unsigned>*head;
        const std::vector<unsigned>*weight;
        const std::vector<unsigned>*ranking_mapper;
    };
}










#endif //ROUTINGKIT_MODIFIED_DIJKSTRA_H
