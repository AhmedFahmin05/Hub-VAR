#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <routingkit/id_queue.h>
#include <routingkit/constants.h>
#include <routingkit/timestamp_flag.h>
#include <vector>
#include <iostream>
#include <set>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sstream>
#include <memory.h>
#include <unordered_map>
#include "../src/tree_node.h"

namespace RoutingKit{
    using namespace std;
    class ScalarGetWeight{
    public:
        explicit ScalarGetWeight(const std::vector<unsigned>&weight):weight(&weight){}

        unsigned operator()(unsigned arc, unsigned departure_time)const{
            (void)departure_time;
            return (*weight)[arc];
        }

    private:
        const std::vector<unsigned>*weight;
    };

    class Dijkstra{
    public:
        unsigned number_of_nodes_expanded;
        unsigned number_of_nodes_generated;
        unsigned number_of_nodes;
        unsigned number_of_landmark;
        unsigned number_of_landmark_hit;
        unsigned landmark_start;


        Dijkstra():first_out(nullptr){}


        Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail, const std::vector<unsigned>&head,const std::vector<unsigned>&weight):
                tentative_distance(first_out.size()-1),
                predecessor_arc(first_out.size()-1),
                was_popped(first_out.size()-1),
                queue(first_out.size()-1),
                first_out(&first_out),
                tail(&tail),
                head(&head),
                was_distance_set(first_out.size()-1),weight(&weight){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());
            allowed.resize(first_out.size()-1);
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            allowed_unsigned.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
        }

        Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail, const std::vector<unsigned>&head):
                tentative_distance(first_out.size()-1),
                predecessor_arc(first_out.size()-1),
                was_popped(first_out.size()-1),
                queue(first_out.size()-1),
                first_out(&first_out),
                tail(&tail),
                head(&head),was_distance_set(first_out.size()-1){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());
            allowed.resize(first_out.size()-1);
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            allowed_unsigned.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
        }

        Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail, const std::vector<unsigned>&head,const std::vector<unsigned>&weight,const std::vector<unsigned>&rank):
                tentative_distance(first_out.size()-1),
                predecessor_arc(first_out.size()-1),
                was_popped(first_out.size()-1),
                queue(first_out.size()-1),
                first_out(&first_out),
                tail(&tail),
                head(&head),
                weight(&weight),
                rank_mapper(&rank),was_distance_set(first_out.size()-1){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());
            allowed.resize(first_out.size()-1);
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            allowed_unsigned.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
        }


        Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>&rank,const std::vector<unsigned>&landmark_list):
                tentative_distance(first_out.size()-1),
                predecessor_arc(first_out.size()-1),
                was_popped(first_out.size()-1),
                queue(first_out.size()-1),
                first_out(&first_out),
                tail(&tail),
                head(&head),
                weight(&weight),
                rank_mapper(&rank),was_distance_set(first_out.size()-1){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());
            allowed.resize(first_out.size()-1);
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            allowed_unsigned.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
        }

        Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail, const std::vector<unsigned>&head, const std::vector<unsigned>&weight, const std::vector<unsigned>&rank, const std::vector<unsigned*>&landmark_pointer,
                 unsigned landmark_number):
                tentative_distance(first_out.size()-1),
                predecessor_arc(first_out.size()-1),
                was_popped(first_out.size()-1),
                queue(first_out.size()-1),
                first_out(&first_out),
                tail(&tail),
                head(&head),
                weight(&weight),
                rank_mapper(&rank),landmark_pointer(&landmark_pointer),was_distance_set(first_out.size()-1){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());
            allowed.resize(first_out.size()-1);
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            allowed_unsigned.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
            number_of_landmark  = landmark_number;
        }






        Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper,
                 unsigned landmark_number,unsigned l_start,const std::vector<unsigned*>& landmark_pointer):
                tentative_distance(first_out.size()-1),
                predecessor_arc(first_out.size()-1),
                was_popped(first_out.size()-1),
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper),landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            allowed.resize(first_out.size()-1);
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            allowed_unsigned.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
//            row_size = number_of_nodes*4;
//            file = open(file_name.c_str(), O_RDONLY);
//            memmap = map;
            //preserve row_size memory + 4096;
//
//            preserved_memory= (mmap(nullptr, row_size+4096, PROT_NONE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
//            preserved_memory = (mmap(nullptr, row_size+4096, PROT_READ, MAP_PRIVATE,file, 0));
        }




//    unsigned get_max_landmark_distance(unsigned start, unsigned end){
//        unsigned max = 0 ;
//        for(int i = 0 ; i < 16; i++){
//
//            unsigned start_distance = (*landmark_list)[start+i*number_of_nodes];
//            unsigned end_distance =(*landmark_list)[end+i*number_of_nodes];
//            unsigned landmard_distance = abs((long long int)start_distance - end_distance);
//            if(max < landmard_distance){
//                max = landmard_distance;
//            }
//        }
//        return max;
//    }


//    unsigned get_min_landmark_upperbound(unsigned start, unsigned end){
//        unsigned min = inf_weight;
//        unsigned start_index = start * 16;
//        unsigned end_index =  end * 16;
//        for(int i = 0 ; i < 16; i++){
//            unsigned landmard_distance = (*landmark_list)[start_index+i] +  (*landmark_list)[end_index+i];
//            if(min > landmard_distance){
//                min = landmard_distance;
//            }
//        }
//        return min;
//    }

        bool pruned_by_landmark (unsigned start, unsigned end, unsigned current_cost){
            unsigned* startPtr = (*landmark_pointer)[start];
            unsigned* endPtr = (*landmark_pointer)[end];
            for(int i = 0 ; i < number_of_landmark; i++){
                if( *(startPtr+i*2) +  *(endPtr+i*2) <= current_cost){
                    return true;
                }
            }
            return false;
        }


        Dijkstra&reset(){
            queue.clear();
            was_popped.reset_all();
            return *this;
        }

        Dijkstra&reset(const std::vector<unsigned>&first_out, const std::vector<unsigned>&tail, const std::vector<unsigned>&head){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == tail.size());
            assert(first_out.back() == head.size());

            if(this->first_out != nullptr && first_out.size() == this->first_out->size()){
                this->first_out = &first_out;
                this->head = &head;
                this->tail = &tail;
                queue.clear();
                was_popped.reset_all();
                return *this;
            }else{
                this->first_out = &first_out;
                this->head = &head;
                this->tail = &tail;
                tentative_distance.resize(first_out.size()-1);
                predecessor_arc.resize(first_out.size()-1);
                was_popped = TimestampFlags(first_out.size()-1);
                queue = MinIDQueue(first_out.size()-1);
                return *this;
            }
        }

        Dijkstra&add_source(unsigned id, unsigned departure_time = 0){
            assert(id < first_out->size()-1);
            tentative_distance[id] = departure_time;
            predecessor_arc[id] = invalid_id;
            queue.push({id, departure_time});
            return *this;
        }

        bool is_finished()const{
            return queue.empty();
        }

        bool was_node_reached(unsigned x)const{
            assert(x < first_out->size()-1);
            return was_popped.is_set(x);
        }

        struct SettleResult{
            unsigned node;
            unsigned distance;
        };

        template<class GetWeightFunc>
        SettleResult settle(const GetWeightFunc&get_weight){
            assert(!is_finished());

            auto p = queue.pop();
            tentative_distance[p.id] = p.key;
            was_popped.set(p.id);

            for(unsigned a=(*first_out)[p.id]; a<(*first_out)[p.id+1]; ++a){
                if(!was_popped.is_set((*head)[a])){
                    unsigned w = get_weight(a, p.key);
                    if(w < inf_weight){
                        if(queue.contains_id((*head)[a])){
                            if(queue.decrease_key({(*head)[a], p.key + w})){
                                predecessor_arc[(*head)[a]] = a;
                            }
                        } else {
                            queue.push({(*head)[a], p.key + w});
                            predecessor_arc[(*head)[a]] = a;
                        }
                    }
                }
            }
            return SettleResult{p.id, p.key};
        }


        template<class GetWeightFunc>
        SettleResult create_spt(const GetWeightFunc&get_weight, vector<RoutingKit::tree_node>& spt){
            assert(!is_finished());

            auto p = queue.pop();
            tentative_distance[p.id] = p.key;
            was_popped.set(p.id);
            for(unsigned a=(*first_out)[p.id]; a<(*first_out)[p.id+1]; ++a){
                if(!was_popped.is_set((*head)[a])){
                    unsigned w = get_weight(a, p.key);
                    if(w < inf_weight){
                        if(queue.contains_id((*head)[a])){
                            if(queue.decrease_key({(*head)[a], p.key + w})){
                                predecessor_arc[(*head)[a]] = a;
                                                           }
                        } else {
                            queue.push({(*head)[a], p.key + w});
                            predecessor_arc[(*head)[a]] = a;
                        }
                    }
                }
            }
            spt[p.id].set_id(p.id);
            // spt[p.id].set_parent((*head)[predecessor_arc[p.id]]);
            spt[p.id].set_g(p.key);
            unsigned prev = predecessor_arc[p.id];
            if(prev != invalid_id){
                spt[p.id].set_parent((*tail)[predecessor_arc[p.id]]);
            }
            return SettleResult{p.id, p.key};
        }

        template<class GetWeightFunc>
        SettleResult create_spt(const GetWeightFunc&get_weight, vector<RoutingKit::tree_node>& spt, std::unordered_map<unsigned, unsigned >& idx){
            assert(!is_finished());

            auto p = queue.pop();
            tentative_distance[p.id] = p.key;
            was_popped.set(p.id);
            for(unsigned a=(*first_out)[p.id]; a<(*first_out)[p.id+1]; ++a){
                if(!was_popped.is_set((*head)[a])){
                    unsigned w = get_weight(a, p.key);
                    if(w < inf_weight ){
                        if(queue.contains_id((*head)[a])){
                            if(queue.decrease_key({(*head)[a], p.key + w})){
                                predecessor_arc[(*head)[a]] = a;
                            }
                        } else {
                            queue.push({(*head)[a], p.key + w});
                            predecessor_arc[(*head)[a]] = a;
                        }
                    }
                }
            }
            //spt[p.id].set_id(p.id);
            //spt[p.id].set_g(p.key);
            unsigned prev = predecessor_arc[p.id];
            if(prev != invalid_id){
                //spt[p.id].set_parent((*tail)[predecessor_arc[p.id]]);
                RoutingKit::tree_node current{p.id, (*tail)[predecessor_arc[p.id]], p.key};
                current.set_expanded(0);
                spt.push_back(current);
                idx[p.id] = spt.size();
            } else{
                RoutingKit::tree_node current{p.id, INT32_MAX, p.key};
                current.set_expanded(0);
                spt.push_back(current);
                idx[p.id] = spt.size();
            }
            return SettleResult{p.id, p.key};
        }

//    void reach(const unsigned & next_vertex, unsigned d, const std::set<unsigned>& fmove, const unsigned &number_of_vertices){
//        if(d < dist[next_vertex]){
//            if(queue.contains_id(next_vertex)){
//                queue.decrease_key({next_vertex, d});
//            } else {
//                queue.push({next_vertex,d});
//            }
//            parent[next_vertex] = number_of_vertices;
//            dist[next_vertex] = d;
//            allowed_unsigned[next_vertex]= fmove;
//        }else if(d == dist[next_vertex]){
//            if(parent[next_vertex] > number_of_vertices){
//                allowed_unsigned[next_vertex]= fmove;
//                parent[next_vertex] = number_of_vertices;
//            }else if(parent[next_vertex] == number_of_vertices){
//                for(unsigned fm: fmove){
//                    allowed_unsigned[next_vertex].insert(fm);
//                }
//            }
////            for(unsigned fm: fmove){
////                allowed[next_vertex].insert(fm);
////            }
//        }
//    }
//    template<class GetWeightFunc>
//    const std::vector<std::set<unsigned>>& run(int source_node,const GetWeightFunc&get_weight){
//        std::fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
//        std::fill( allowed_unsigned.begin(), allowed_unsigned.end(),set<unsigned>{ 0xFF});
//        std::fill( parent.begin(), parent.end(),std::numeric_limits<int>::max());
//
//        queue.clear();
//
//        dist[source_node] = 0;
//        unsigned fm_index = 0;
//        for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
//            unsigned weight =  get_weight(a,0);
//            reach((*head)[a], weight, set<unsigned>{fm_index},1);
//            fm_index++;
//        }
//
//        while(!queue.empty()){
//            auto x = queue.pop();
//            for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                unsigned weight = get_weight(a,x.key);
//                reach((*head)[a], x.key + weight, allowed_unsigned[x.id],parent[x.id]+1);
//            }
//        }
//        allowed_unsigned[source_node] ={};
//        return allowed_unsigned;
//    }


        void reach(const unsigned & next_vertex, unsigned d, const std::set<unsigned short>& fmove, const unsigned &number_of_vertices){
            if(d < dist[next_vertex]){
                if(queue.contains_id(next_vertex)){
                    queue.decrease_key({next_vertex, d});
                } else {
                    queue.push({next_vertex,d});
                }
                parent[next_vertex] = number_of_vertices;
                dist[next_vertex] = d;
                allowed[next_vertex] = fmove;
            }else if(d == dist[next_vertex]){
                if(parent[next_vertex] > number_of_vertices){
                    allowed[next_vertex] = fmove;
                    parent[next_vertex] = number_of_vertices;
                }else if(parent[next_vertex] == number_of_vertices){
                    for(const unsigned short& fm: fmove){
                        allowed[next_vertex].insert(fm);
                    }
                }
            }

        }
        template<class GetWeightFunc>
        const std::vector<std::set<unsigned short>>& run(int source_node,const GetWeightFunc&get_weight){
            std::fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
            std::fill( allowed.begin(), allowed.end(),set<unsigned short>{ 0xFF});
            std::fill( parent.begin(), parent.end(),std::numeric_limits<int>::max());

            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned weight =  get_weight(a,0);
                reach((*head)[a], weight, set<unsigned short>{fm_index},1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned weight = get_weight(a,x.key);
                    reach((*head)[a], x.key + weight, allowed[x.id],parent[x.id]+1);
                }
            }
            allowed[source_node] ={};
            return allowed;
        }


        void reach_shortest(const unsigned & next_vertex, unsigned d, const std::set<unsigned short>& fmove, const unsigned &number_of_vertices){
            if(was_distance_set.is_set(next_vertex)){
                if(d <  dist[next_vertex]){
                    queue.decrease_key({next_vertex, d});
                    parent[next_vertex] = number_of_vertices;
                    dist[next_vertex] = d;
                    allowed[next_vertex] = fmove;
                }else if(d ==  dist[next_vertex]){
                    if(parent[next_vertex] > number_of_vertices){
                        allowed[next_vertex] = fmove;
                        parent[next_vertex] = number_of_vertices;
                    }else if(parent[next_vertex] == number_of_vertices){
                        for(const unsigned short& fm: fmove){
                            allowed[next_vertex].insert(fm);
                        }
                    }
                }

            }else{
                queue.push({next_vertex,d});
                parent[next_vertex] = number_of_vertices;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                allowed[next_vertex] = fmove;
            }
        }

        const std::vector<std::set<unsigned short>>& run_shortest(int source_node){
            std::fill( allowed.begin(), allowed.end(),set<unsigned short>{ 0xFF});
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);

            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_shortest((*head)[a], cur_weight, set<unsigned short>{fm_index},1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_shortest((*head)[a], x.key + cur_weight, allowed[x.id],parent[x.id]+1);
                }
            }
            allowed[source_node] ={};
            return allowed;
        }


        void reach_greedy_shortest(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned &number_of_vertices){
            if(was_distance_set.is_set(next_vertex)){
                if(d <  dist[next_vertex]){
                    number_of_nodes_generated++;
                    queue.decrease_key({next_vertex, d});
                    parent[next_vertex] = number_of_vertices;
                    dist[next_vertex] = d;
                    d_allowed[next_vertex] = fmove;
                }else if(d == dist[next_vertex]){
                    if(parent[next_vertex] > number_of_vertices){
                        d_allowed[next_vertex] = fmove;
                        parent[next_vertex] = number_of_vertices;
                    }
                }
            }else{
                queue.push({next_vertex,d});
                parent[next_vertex] = number_of_vertices;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                d_allowed[next_vertex] = fmove;
                number_of_nodes_generated ++;
            }

        }


        const std::vector<unsigned short>& run_greedy_shortest_with_landmark(int source_node){
//        if((*rank_mapper)[source_node]>=264330){
//            return run_greedy_shortest(source_node);
//        }
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( d_allowed.begin(), d_allowed.end(),0xFF);

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_greedy_shortest((*head)[a], cur_weight, fm_index,1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
//            if(x.key< get_max_landmark_distance(source_node,x.id)){
//               std::cout<<"error"<<std::endl;
//            }
//            if(x.key>= get_min_landmark_upperbound(source_node,x.id)){
//                continue;
//            }
                if(pruned_by_landmark(source_node,x.id,x.key)){
                    d_allowed[x.id] =0xFF;
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_greedy_shortest((*head)[a], x.key + cur_weight, d_allowed[x.id],parent[x.id]+1);
                }
            }
            d_allowed[source_node] =0xFE;

            return d_allowed;
        }

        const std::vector<unsigned short>& run_greedy_CH_dijsktra_with_landmark(int source_node){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( d_allowed.begin(), d_allowed.end(), 0xFF);
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                if((*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]]){
                    // mask 10 indicate moving up
                    reach_greedy_CH((*head)[a], cur_weight, fm_index,2);
                }else{
                    // mask 01 indicate moving down
                    reach_greedy_CH((*head)[a], cur_weight, fm_index,1);
                }
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                if(pruned_by_landmark(source_node,x.id,x.key)){
                    d_allowed[x.id] =0xFF;
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_greedy_CH((*head)[a], x.key + cur_weight, d_allowed[x.id],1);
                        }
                    }else{
                        if((*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]]){
                            // mask 10 indicate moving up
                            reach_greedy_CH((*head)[a],  x.key + cur_weight, d_allowed[x.id],2);
                        }else{
                            // mask 01 indicate moving down
                            reach_greedy_CH((*head)[a], x.key + cur_weight, d_allowed[x.id],1);
                        }
                    }

                }
            }
            d_allowed[source_node] =0xFE;
            return d_allowed;
        }
        const std::vector<unsigned short>& run_greedy_shortest_with_landmark_check(int source_node,const vector<unsigned>& distance){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( d_allowed.begin(), d_allowed.end(),0xFF);

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_greedy_shortest((*head)[a], cur_weight, fm_index,1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
//            if(x.key< get_max_landmark_distance(source_node,x.id)){
//               std::cout<<"error"<<std::endl;
//            }
//            if(x.key>= get_min_landmark_upperbound(source_node,x.id)){
//                continue;
//            }
                unsigned expect = distance[x.id];

                if(pruned_by_landmark(source_node,x.id,x.key)){
                    d_allowed[x.id] =0xFF;
                    continue;
                }
                if(expect != x.key){
                    bool a =0;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_greedy_shortest((*head)[a], x.key + cur_weight, d_allowed[x.id],parent[x.id]+1);
                }
            }
            d_allowed[source_node] =0xFE;

            return d_allowed;
        }


        const std::vector<unsigned>& run_to_get_fm_distance(int source_node, const vector<unsigned short>& predecessor,vector<unsigned short>& backward_first_move,vector<unsigned short>&first_move){
            std::fill( backward_first_move.begin(), backward_first_move.end(),0xFF);
            std::fill( first_move.begin(), first_move.end(),0xFF);
            std::fill( dist.begin(), dist.end(),inf_weight);

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& predecessor ,const unsigned short& fmove, const unsigned &number_of_vertices) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d <  dist[next_vertex]){
                        queue.decrease_key({next_vertex, d});
                        parent[next_vertex] = number_of_vertices;
                        dist[next_vertex] = d;
                        backward_first_move[next_vertex] =  predecessor;
                        first_move[next_vertex] = fmove;
                    }else if(d == dist[next_vertex]){
                        if(parent[next_vertex] > number_of_vertices){
                            backward_first_move[next_vertex] =  predecessor;
                            first_move[next_vertex] = fmove;
                            parent[next_vertex] = number_of_vertices;
                        }
                    }
                }else{
                    queue.push({next_vertex,d});
                    parent[next_vertex] = number_of_vertices;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    backward_first_move[next_vertex] =  predecessor;
                    first_move[next_vertex] = fmove;
                }
            };
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_nodes((*head)[a], cur_weight, predecessor[a],fm,1);
                fm++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_nodes((*head)[a], x.key + cur_weight, predecessor[a],first_move[x.id],parent[x.id]+1);
                }
            }
            backward_first_move[source_node] =0xFE;
            first_move[source_node] =0xFE;
            return dist;
        }






        const std::vector<unsigned >& get_distance_table(){
            return dist;
        }

        const std::vector<unsigned short>& get_fm_and_distance(unsigned source_node){
            std::fill( d_allowed.begin(), d_allowed.end(),0xFF);
//            std::fill( dist.begin(), dist.end(),std::numeric_limits<unsigned>::max());
            memset(&dist[0],-1,dist.size()*4);

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove,  const unsigned &number_of_vertices) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d <  dist[next_vertex]){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        parent[next_vertex] = number_of_vertices;
                        dist[next_vertex] = d;
                        d_allowed[next_vertex] = fmove;
                    }else if(d == dist[next_vertex]){
                        if(parent[next_vertex] > number_of_vertices){
                            d_allowed[next_vertex] = fmove;
                            parent[next_vertex] = number_of_vertices;
                        }
                    }
                }else{
                    queue.push({next_vertex,d});
                    parent[next_vertex] = number_of_vertices;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    d_allowed[next_vertex] = fmove;
                    number_of_nodes_generated ++;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_nodes((*head)[a], cur_weight, fm_index,1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_nodes((*head)[a], x.key + cur_weight, d_allowed[x.id],parent[x.id]+1);
                }
            }
            d_allowed[source_node] =0xFE;
            return d_allowed;
        }


        const std::vector<unsigned short>& run_greedy_shortest(int source_node){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( d_allowed.begin(), d_allowed.end(),0xFF);

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_greedy_shortest((*head)[a], cur_weight, fm_index,1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_greedy_shortest((*head)[a], x.key + cur_weight, d_allowed[x.id],parent[x.id]+1);
                }
            }
            d_allowed[source_node] =0xFE;

            return d_allowed;
        }













        void reach_ori_CH(const unsigned & next_vertex, unsigned d, const std::set<unsigned short>& fmove, const unsigned &up_down_signal){

            if(d < dist[next_vertex]){
                if(queue.contains_id(next_vertex)){
                    queue.decrease_key({next_vertex, d});
                } else {
                    queue.push({next_vertex,d});
                }
                up_or_down[next_vertex] = up_down_signal;
                dist[next_vertex] = d;
                allowed[next_vertex] = fmove;
            }else if(d == dist[next_vertex]){
                up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                for(const unsigned short & fm: fmove){
                    allowed[next_vertex].insert(fm);
                }
            }

        }

        const std::vector<std::set<unsigned short>>& run_ori_CH_dijsktra(int source_node){
            std::fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
            std::fill( allowed.begin(), allowed.end(),set<unsigned short>{ 0xFF});
            std::fill( up_or_down.begin(), up_or_down.end(),0);

            queue.clear();

            dist[source_node] = 0;

            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                if((*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]]){
                    // mask 10 indicate moving up
                    reach_ori_CH((*head)[a], cur_weight, set<unsigned short>{fm_index},2);
                }else{
                    // mask 01 indicate moving down
                    reach_ori_CH((*head)[a], cur_weight, set<unsigned short>{fm_index},1);
                }
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_ori_CH((*head)[a], x.key + cur_weight, allowed[x.id],1);
                        }
                    }else{
                        if((*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]]){
                            // mask 10 indicate moving up
                            reach_ori_CH((*head)[a],  x.key + cur_weight, allowed[x.id],2);
                        }else{
                            // mask 01 indicate moving down
                            reach_ori_CH((*head)[a], x.key + cur_weight, allowed[x.id],1);
                        }
                    }

                }
            }
            allowed[source_node] ={};
            return allowed;
        }

        void reach_CH(const unsigned & next_vertex, unsigned d, const std::set<unsigned short>& fmove, const unsigned &up_down_signal){
            if(was_distance_set.is_set(next_vertex)){
                if(d < dist[next_vertex] ){
                    queue.decrease_key({next_vertex, d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    allowed[next_vertex] = fmove;
                }else if(d == dist[next_vertex] ){
                    up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    for(const unsigned short & fm: fmove){
                        allowed[next_vertex].insert(fm);
                    }
                }
            }else{
                queue.push({next_vertex,d});
                up_or_down[next_vertex] = up_down_signal;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                allowed[next_vertex] = fmove;
            }
        }

        const std::vector<std::set<unsigned short>>& run_CH_dijsktra(int source_node){
            std::fill( allowed.begin(), allowed.end(),set<unsigned short>{ 0xFF});
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                if((*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]]){
                    // mask 10 indicate moving up
                    reach_CH((*head)[a], cur_weight, set<unsigned short>{fm_index},2);
                }else{
                    // mask 01 indicate moving down
                    reach_CH((*head)[a], cur_weight, set<unsigned short>{fm_index},1);
                }
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_CH((*head)[a], x.key + cur_weight, allowed[x.id],1);
                        }
                    }else{
                        if((*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]]){
                            // mask 10 indicate moving up
                            reach_CH((*head)[a],  x.key + cur_weight, allowed[x.id],2);
                        }else{
                            // mask 01 indicate moving down
                            reach_CH((*head)[a], x.key + cur_weight, allowed[x.id],1);
                        }
                    }

                }
            }
            allowed[source_node] ={};
            return allowed;
        }




        void reach_greedy_CH(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned &up_down_signal){
            if(was_distance_set.is_set(next_vertex)){
                if(d < dist[next_vertex]){
                    number_of_nodes_generated++;
                    queue.decrease_key({next_vertex, d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    d_allowed[next_vertex] = fmove;
                }else if(d == dist[next_vertex]){
                    up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                }
            }else{
                number_of_nodes_generated++;
                queue.push({next_vertex,d});
                up_or_down[next_vertex] = up_down_signal;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                d_allowed[next_vertex] = fmove;
            }

        }

        const std::vector<unsigned short>& run_greedy_CH_dijsktra(int source_node){
            std::fill( d_allowed.begin(), d_allowed.end(), 0xFF);
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                if((*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]]){
                    // mask 10 indicate moving up
                    reach_greedy_CH((*head)[a], cur_weight, fm_index,2);
                }else{
                    // mask 01 indicate moving down
                    reach_greedy_CH((*head)[a], cur_weight, fm_index,1);
                }
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_greedy_CH((*head)[a], x.key + cur_weight, d_allowed[x.id],1);
                        }
                    }else{
                        if((*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]]){
                            // mask 10 indicate moving up
                            reach_greedy_CH((*head)[a],  x.key + cur_weight, d_allowed[x.id],2);
                        }else{
                            // mask 01 indicate moving down
                            reach_greedy_CH((*head)[a], x.key + cur_weight, d_allowed[x.id],1);
                        }
                    }

                }
            }
            d_allowed[source_node] =0xFE;
            return d_allowed;
        }

        void reach_all_move(const unsigned & next_vertex, unsigned d, const std::set<unsigned short>& fmove){
            if(d < dist[next_vertex]){
                if(queue.contains_id(next_vertex)){
                    queue.decrease_key({next_vertex, d});
                } else {
                    queue.push({next_vertex,d});
                }
                dist[next_vertex] = d;
                allowed[next_vertex] = fmove;
            }else if(d == dist[next_vertex]){
                for(unsigned short fm: fmove){
                    allowed[next_vertex].insert(fm);
                }
            }
        }



        const std::vector<std::set<unsigned short>>& run_to_get_all_move(int source_node){
            std::fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
            std::fill( allowed.begin(), allowed.end(),set<unsigned short>{ 0xFF});

            queue.clear();

            dist[source_node] = 0;
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_all_move((*head)[a], cur_weight, set<unsigned short>{fm_index});
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_all_move((*head)[a], x.key + cur_weight, allowed[x.id]);
                }
            }
            allowed[source_node] ={};
            return allowed;
        }



        const std::vector<unsigned short>&run_daniel_dijkstra(int source_node){
            std::fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
            std::fill(d_allowed.begin(), d_allowed.end(), 0);
            queue.clear();
            dist[source_node] = 0;
            d_allowed[source_node] = 0;
            auto d_reach = [&](unsigned v, unsigned d, unsigned short first_move){
                if(d < dist[v]){
                    if(queue.contains_id(v)){
                        queue.decrease_key({v, d});
                    } else {
                        queue.push({v,d});
                    }
                    dist[v] = d;
                    d_allowed[v] = first_move;
                }else if(d == dist[v]){
                    d_allowed[v] |= first_move;
                }
            };

            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                d_reach((*head)[a], cur_weight, 1<<fm_index);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    d_reach((*head)[a], x.key + cur_weight, d_allowed[x.id]);
                }
            }
            return d_allowed;
        }



        vector<unsigned> get_distance_list( ) {
            return dist;
        }

        unsigned get_distance_to(unsigned x) const {
            assert(x < first_out->size()-1);
            if(was_popped.is_set(x))
                return tentative_distance[x];
            else
                return inf_weight;
        }

        std::vector<unsigned>get_node_path_to(unsigned x) const {
            assert(x < first_out->size()-1);
            std::vector<unsigned>path;
            if(was_node_reached(x)){
                assert(was_node_reached(x));

                unsigned p;
                while(p = predecessor_arc[x], p != invalid_id){
                    path.push_back(x);
                    x = (*tail)[p];
                }
                path.push_back(x);
                std::reverse(path.begin(), path.end());
            }
            return path;
        }


        unsigned get_first_move (unsigned source,unsigned x) const {
            assert(x < first_out->size()-1);
            std::vector<unsigned>path;
            if(was_node_reached(x)){
                assert(was_node_reached(x));

                unsigned p;
                while(p = predecessor_arc[x], p != invalid_id){
                    path.push_back(x);
                    x = (*tail)[p];
                }
                path.push_back(x);
//            std::reverse(path.begin(), path.end());
                unsigned first_move = path[path.size()-2];
                for(unsigned arc=(*first_out)[source]; arc<(*first_out)[source+1]; ++arc){
                    if((*head)[arc] == first_move){
                        return arc-(*first_out)[source];
                    }
                }
                std::cout<<"arc id not found" << std::endl;
                return -1;
            }else{
                return  -1;
            }
        }

        std::vector<unsigned>get_first_move_table(unsigned source){
            std::vector<unsigned>first_move_table = std::vector<unsigned>();
            for(int i = 0 ; i< first_out->size() - 1; i++){
                if(i == source){
                    first_move_table.push_back(0);
                }else{
                    unsigned first_move = get_first_move(source,i);
                    first_move_table.push_back(first_move+2);
                }
            }
            return first_move_table;
        }

        std::vector<unsigned>get_arc_path_to(unsigned x) const {
            assert(x < first_out->size()-1);
            std::vector<unsigned>path;
            if(was_node_reached(x)){
                unsigned p;
                while(p = predecessor_arc[x], p != invalid_id){
                    path.push_back(p);
                    x = (*tail)[p];
                }
                std::reverse(path.begin(), path.end());
            }
            return path;
        }


        void full_search(unsigned source, const std::vector<unsigned>& weight){
            reset();
            add_source(source);
            while(!is_finished()){
                settle(ScalarGetWeight(weight));
            }
        }


        vector<unsigned int> get_all_distance(unsigned source, const std::vector<unsigned>& weight){
            reset();
            add_source(source);
            while(!is_finished()){
                settle(ScalarGetWeight(weight));
            }
            vector<unsigned> distance_list = vector<unsigned>(first_out->size()-1);
            for(int i = 0; i < distance_list.size(); i++){
                distance_list[i] = get_distance_to(i);
            }
            return distance_list;
        }



//        auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove,  const unsigned &number_of_vertices) {
//            if(was_distance_set.is_set(next_vertex)){
//                if(d <  dist[next_vertex]){
//                    number_of_nodes_generated++;
//                    queue.decrease_key({next_vertex, d});
//                    parent[next_vertex] = number_of_vertices;
//                    dist[next_vertex] = d;
//                    d_allowed[next_vertex] = fmove;
//                }else if(d == dist[next_vertex]){
//                    if(parent[next_vertex] > number_of_vertices){
//                        d_allowed[next_vertex] = fmove;
//                        parent[next_vertex] = number_of_vertices;
//                    }
//                }
//            }else{
//                queue.push({next_vertex,d});
//                parent[next_vertex] = number_of_vertices;
//                dist[next_vertex] = d;
//                was_distance_set.set(next_vertex);
//                d_allowed[next_vertex] = fmove;
//                number_of_nodes_generated ++;
//            }
//        };

        void reach_landmark_sequential(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned &number_of_vertices){


            if(d < dist[next_vertex]){
                number_of_nodes_generated++;
                if(queue.contains_id(next_vertex)){
                    queue.decrease_key({next_vertex, d});
                }else{
                    queue.push({next_vertex, d});
                }
                parent[next_vertex] = number_of_vertices;
                dist[next_vertex] = d;
                d_allowed[next_vertex] = fmove;
            }else if(d == dist[next_vertex]){
                if(parent[next_vertex] > number_of_vertices){
                    d_allowed[next_vertex] = fmove;
                    parent[next_vertex] = number_of_vertices;
                }
            }
        }

        const std::vector<unsigned short>& run_with_reached_landmark_sequential_with_disk_cache(unsigned source_node){
            auto updates_fm_and_distance= [&](const unsigned & landmark_id,const unsigned distance, unsigned fm) {
                number_of_landmark_hit++;
                unsigned* contents = (*landmark_pointer)[landmark_id];
                auto distPtr = dist.begin();
                auto fmPtr = d_allowed.begin();
                for(unsigned i = 0; i < number_of_nodes; i ++ ){
                    unsigned current_distance = distance + *(contents);
                    if (current_distance < *(distPtr)) {
                        *distPtr = current_distance;
                        *fmPtr = fm;
                    }
                    distPtr++;
                    fmPtr++;
                    contents++;
                }
            };
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            number_of_landmark_hit = 0;
            std::fill(dist.begin(), dist.end(), std::numeric_limits<unsigned>::max());
//            memset(&dist[0],-1,dist.size()*4);
            std::fill( d_allowed.begin(), d_allowed.end(), 0xFF);

            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;

            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                unsigned cur_weight = (*weight)[a];
                reach_landmark_sequential((*head)[a], cur_weight, fm_index,1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                if(x.key > dist[x.id]){
                    continue;
                }
                unsigned current_rank = (*rank_mapper)[x.id];

                if(current_rank >= landmark_start) {
                    assert(current_rank - landmark_start < (*rank_mapper).size());
                    updates_fm_and_distance(current_rank - landmark_start, x.key,
                                            d_allowed[x.id]);
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    unsigned cur_weight = (*weight)[a];
                    reach_landmark_sequential((*head)[a], x.key + cur_weight, d_allowed[x.id],parent[x.id]+1);
                }
            }
            d_allowed[source_node] =0xFE;
            return d_allowed;
        }



    private:
        std::vector<unsigned>tentative_distance;
        std::vector<unsigned>predecessor_arc;

        std::vector<unsigned>dist;
        std::vector<std::set<unsigned short>>allowed;
        std::vector<std::set<unsigned >>allowed_unsigned;
        std::vector<unsigned short>d_allowed;
//    std::vector<unsigned short>backward_allowed;
        std::vector<unsigned>parent;
        std::vector<unsigned>up_or_down;
        TimestampFlags was_popped;
        MinIDQueue queue;
        TimestampFlags was_distance_set;

        const std::vector<unsigned>*first_out;
        const std::vector<unsigned>*tail;
        const std::vector<unsigned>*head;
        const std::vector<unsigned>*weight;
        const std::vector<unsigned>*rank_mapper;
        const std::vector<unsigned*>*landmark_pointer;
    };


}

#endif

