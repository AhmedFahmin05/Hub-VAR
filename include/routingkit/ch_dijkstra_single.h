//
// Created by Bojie Shen on 1/11/20.
//

#ifndef ROUTINGKIT_CH_DIJKSTRA_SINGLE_H
#define ROUTINGKIT_CH_DIJKSTRA_SINGLE_H
#include <routingkit/id_queue.h>
#include <routingkit/constants.h>
#include <routingkit/timestamp_flag.h>
#include <vector>
#include <iostream>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sstream>
#include <memory.h>

using namespace std;
namespace RoutingKit{

    class Ch_Dijkstra_Single{
    public:
        unsigned number_of_nodes_expanded;
        unsigned number_of_nodes_generated;
        unsigned number_of_nodes;
        unsigned number_of_landmark;
        unsigned number_of_landmark_hit;
        unsigned landmark_start;
        unsigned long long row_size;
        int file;
        std::ifstream landmark_cache;
        std::vector<unsigned>::iterator lstartPtr;
        std::vector<unsigned>::iterator lendPtr;
        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_nodes = first_out.size()-1;
        }
        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight, const std::vector<unsigned>& rank_mapper):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_nodes = first_out.size()-1;
        }

        Ch_Dijkstra_Single(const std::vector<unsigned>& rank_mapper,
                           const std::vector<unsigned>&up_first_out, const std::vector<unsigned>&up_head, const std::vector<unsigned>&up_weight,
                           const std::vector<unsigned>&down_first_out, const std::vector<unsigned>&down_head, const std::vector<unsigned>&down_weight
        ):
                queue(up_first_out.size()-1),
                was_distance_set(up_first_out.size()-1),
                rank_mapper(&rank_mapper),
                up_first_out(&up_first_out),
                up_head(&up_head),
                up_weight(&up_weight),
                down_first_out(&down_first_out),
                down_head(&down_head),
                down_weight(&down_weight)
        {
            first_move.resize(up_first_out.size()-1);
            up_or_down.resize(up_first_out.size()-1);
            dist.resize(up_first_out.size()-1);
            number_of_nodes = up_first_out.size()-1;
        }
        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper, const std::vector<unsigned*>&landmark_pointer,
                           unsigned landmark_number):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper),landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            number_of_nodes = first_out.size()-1;
        }

        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight, const std::vector<unsigned*>&landmark_pointer,
                           unsigned landmark_number):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            number_of_nodes = first_out.size()-1;
        }

        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper, const std::vector<unsigned*>&landmark_pointer,
                           unsigned landmark_number,unsigned l_start):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper),landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
            landmark_fm.resize(number_of_landmark);
            landmark_reached.resize(number_of_landmark);
            landmark_cached.resize(number_of_nodes);
        }

        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper, const std::vector<unsigned>&landmark_list,
                           unsigned landmark_number,unsigned l_start):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper),landmark_list(&landmark_list){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
            landmark_fm.resize(number_of_landmark);
            landmark_reached.resize(number_of_landmark);
        }

        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper,
                           unsigned landmark_number,unsigned l_start,const string& file_name):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
            landmark_fm.resize(number_of_landmark);
            landmark_reached.resize(number_of_landmark);
            landmark_cached.resize(number_of_nodes);
            landmark_cache.open(file_name, std::ios::binary);
            lstartPtr = landmark_cached.begin();
            lendPtr = landmark_cached.end();
            row_size = number_of_nodes*4;
            file = open(file_name.c_str(), O_RDONLY);
            //preserve row_size memory + 4096;
//
//            preserved_memory= (mmap(nullptr, row_size+4096, PROT_NONE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
//            preserved_memory = (mmap(nullptr, row_size+4096, PROT_READ, MAP_PRIVATE,file, 0));
        }


        Ch_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper,
                           unsigned landmark_number,unsigned l_start,const std::vector<unsigned*>& landmark_pointer):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper),landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
            landmark_fm.resize(number_of_landmark);
            landmark_reached.resize(number_of_landmark);
            landmark_cached.resize(number_of_nodes);
//            landmark_cache.open(file_name, std::ios::binary);
            lstartPtr = landmark_cached.begin();
            lendPtr = landmark_cached.end();
            row_size = number_of_nodes*4;
//            file = open(file_name.c_str(), O_RDONLY);
//            memmap = map;
            //preserve row_size memory + 4096;
//
//            preserved_memory= (mmap(nullptr, row_size+4096, PROT_NONE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0));
//            preserved_memory = (mmap(nullptr, row_size+4096, PROT_READ, MAP_PRIVATE,file, 0));
        }
//


        Ch_Dijkstra_Single(const std::vector<unsigned>& rank_mapper,
                           unsigned landmark_number,unsigned l_start,const std::vector<unsigned*>& landmark_pointer,
                           const std::vector<unsigned>&up_first_out, const std::vector<unsigned>&up_head, const std::vector<unsigned>&up_weight,
                           const std::vector<unsigned>&down_first_out, const std::vector<unsigned>&down_head, const std::vector<unsigned>&down_weight):
                queue(up_first_out.size()-1),
                was_distance_set(up_first_out.size()-1),
                rank_mapper(&rank_mapper),landmark_pointer(&landmark_pointer),
                up_first_out(&up_first_out),
                up_head(&up_head),
                up_weight(&up_weight),
                down_first_out(&down_first_out),
                down_head(&down_head),
                down_weight(&down_weight){

            first_move.resize(up_first_out.size()-1);
            up_or_down.resize(up_first_out.size()-1);
            dist.resize(up_first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = up_first_out.size()-1;
            landmark_fm.resize(number_of_landmark);
            landmark_reached.resize(number_of_landmark);
            landmark_cached.resize(number_of_nodes);
            lstartPtr = landmark_cached.begin();
            lendPtr = landmark_cached.end();
            row_size = number_of_nodes*4;
        }


        void fast_load_vector(const std::string&file_name, std::vector<unsigned>&vec){
            // assume that must exist

            std::ifstream f(file_name, std::ios::binary);

            f.read(reinterpret_cast<char*>(vec.data()), vec.size()*sizeof(unsigned));
//            f.read(reinterpret_cast<char*>(&vec.front()), vec.size()*sizeof(unsigned));
            f.close();
        }

        void reach(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned char &up_down_signal){
            if(was_distance_set.is_set(next_vertex)){
                if(d < dist[next_vertex]){
                    number_of_nodes_generated++;
                    queue.decrease_key({next_vertex, d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    first_move[next_vertex] = fmove;
                }else if(d == dist[next_vertex]){
                    if(up_down_signal == 2 ){
                        // always take up first move ;
                        first_move[next_vertex] = fmove;
                    }
                    up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                }
            }else{

                number_of_nodes_generated++;
                queue.push({next_vertex,d});
                up_or_down[next_vertex] = up_down_signal;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                first_move[next_vertex] = fmove;
            }
        }

        const std::vector<unsigned short>& run(unsigned source_node){
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach((*head)[a], x.key + (*weight)[a], first_move[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]]? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }



        const std::vector<unsigned short>& run_with_up_down(unsigned source_node){
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
//            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
//                reach((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
//                fm_index++;
//            }

            for(unsigned a=(*up_first_out)[source_node]; a<(*up_first_out)[source_node+1]; ++a){
                reach((*up_head)[a], (*up_weight)[a], fm_index++,2);
                fm_index++;
            }
            for(unsigned a=(*down_first_out)[source_node]; a<(*down_first_out)[source_node+1]; ++a){
                reach((*down_head)[a], (*down_weight)[a], fm_index++,1);
                fm_index++;
            }



            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                if(up_or_down[x.id] == 1){
                    // mask 01 indicate moving down
                    for(unsigned a=(*down_first_out)[x.id]; a<(*down_first_out)[x.id+1]; ++a){
                        reach((*down_head)[a], x.key + (*down_weight)[a], first_move[x.id],1);
                    }
                }else{
                    for(unsigned a=(*up_first_out)[x.id]; a<(*up_first_out)[x.id+1]; ++a){
                        reach((*up_head)[a], x.key + (*up_weight)[a], first_move[x.id],2);
                    }
                    for(unsigned a=(*down_first_out)[x.id]; a<(*down_first_out)[x.id+1]; ++a){
                        reach((*down_head)[a], x.key + (*down_weight)[a], first_move[x.id],1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }

        const std::vector<unsigned short>& run_without_mapper(unsigned source_node){
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach((*head)[a], (*weight)[a], fm_index,source_node < (*head)[a] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(x.id  > (*head)[a]){
                            // mask 01 indicate moving down
                            reach((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach((*head)[a], x.key + (*weight)[a], first_move[x.id],x.id< (*head)[a]? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }

        bool pruned_by_landmark (unsigned start, unsigned end, unsigned current_cost){
            unsigned* startPtr = (*landmark_pointer)[start];
            unsigned* endPtr = (*landmark_pointer)[end];
            for(unsigned i = 0 ; i < number_of_landmark; i++){
                if( *(startPtr+i*2) +  *(endPtr+i*2) <= current_cost){
                    return true;
                }
            }
            return false;
        }


        const std::vector<unsigned short>& run_with_landmark(unsigned source_node){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                if(pruned_by_landmark(source_node,x.id,x.key)){
                    first_move[x.id] =0xFF;
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach((*head)[a], x.key + (*weight)[a], first_move[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }

        bool pruned_by_landmark_sequential (unsigned start, unsigned end, unsigned current_cost){
            for(unsigned i = 0 ; i < number_of_landmark; i++){
                unsigned* startPtr = (*landmark_pointer)[i];
                if( *(startPtr+start) +  *(startPtr+end) <= current_cost){
                    return true;
                }
            }
            return false;
        }


        bool pruned_by_recorded_landmark (unsigned start, unsigned end, unsigned current_cost){
            unsigned* startPtr = (*landmark_pointer)[start];
            unsigned* endPtr = (*landmark_pointer)[end];
            for(unsigned i = 0; i < number_of_reached_landmarks;i++){
                if( *(startPtr+landmark_reached[i]) +  *(endPtr+landmark_reached[i]) <= current_cost){
                    return true;
                }
            }
            return false;
        }

        void reach_landmark(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned char &up_down_signal){
            if(was_distance_set.is_set(next_vertex)){
                if(d < dist[next_vertex]){
                    number_of_nodes_generated++;
                    queue.decrease_key({next_vertex, d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    first_move[next_vertex] = fmove;
                }else if(d == dist[next_vertex]){
                    up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                }
            }else{
                number_of_nodes_generated++;
                queue.push({next_vertex,d});
                up_or_down[next_vertex] = up_down_signal;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                first_move[next_vertex] = fmove;
            }
        }

        void reach_landmark_sequential(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned char &up_down_signal){
            if(d < dist[next_vertex]){
                number_of_nodes_generated++;
                if(queue.contains_id(next_vertex)){
                    queue.decrease_key({next_vertex, d});
                }else{
                    queue.push({next_vertex, d});
                }
                up_or_down[next_vertex] = up_down_signal;
                dist[next_vertex] = d;
                first_move[next_vertex] = fmove;
            }else if(d == dist[next_vertex]){
                up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
            }
        }


        const std::vector<unsigned short>& run_with_reached_landmark_sequential(unsigned source_node){
            auto updates_fm_and_distance= [&](const unsigned & landmark_id,const unsigned distance, unsigned fm) {
                number_of_landmark_hit++;

                unsigned* startPtr = (*landmark_pointer)[landmark_id];
                unsigned* endPtr = startPtr + number_of_nodes;
                auto distPtr = dist.begin();
                auto fmPtr = first_move.begin();
                for(; startPtr != endPtr ; startPtr++){
                    unsigned current_distance = distance + *(startPtr);
                    if (current_distance < *(distPtr)) {
                        *distPtr = current_distance;
                        *fmPtr = fm;
                    }
                    distPtr++;
                    fmPtr++;
                }


            };

//            auto pruned_by_reached_landmarks= [&](const unsigned & current_id,const unsigned& current_distance) {
////                for(unsigned i = 0 ; i < number_of_reached_landmarks; i++){
////                    landmark_cached = load_vector<unsigned>("dataset/NY_cached/"+to_string(landmark_reached[i]));
////                    auto startPtr = landmark_cached.begin();
//////                    unsigned* startPtr = (*landmark_pointer)[landmark_reached[i]];
////                    if(*(startPtr + source_node) + *(startPtr + current_id) <= current_distance){
////                        return true;
////                    }
////                }
//                return false;
//            };
            number_of_reached_landmarks = 0 ;
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            number_of_landmark_hit = 0;
            std::fill(dist.begin(), dist.end(), std::numeric_limits<unsigned>::max());
            std::fill( first_move.begin(), first_move.end(), 0xFF);

            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_landmark_sequential((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
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
//                        updates_fm_and_distance(number_of_landmark - 1 - (current_rank - landmark_start), x.key,
//                                                first_move[x.id]);
//                        landmark_reached[number_of_reached_landmarks] =
//                                number_of_landmark - 1 - (current_rank - landmark_start);
                    updates_fm_and_distance(current_rank - landmark_start, x.key,
                                            first_move[x.id]);
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(current_rank   > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_landmark_sequential((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach_landmark_sequential((*head)[a], x.key + (*weight)[a], first_move[x.id],current_rank  < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }



        const std::vector<unsigned short>& run_with_reached_landmark_sequential_with_disk_cache(unsigned source_node){
//            auto updates_fm_and_distance= [&](const unsigned & landmark_id,const unsigned distance, unsigned fm) {
////                vector<unsigned > cached = vector<unsigned>(number_of_nodes);
////                fast_load_vector("dataset/NY_cached/"+to_string(1),cached);
//                number_of_landmark_hit++;
//                landmark_cache.seekg (4*(unsigned long long)number_of_nodes*landmark_id, std::ios::beg);
//                landmark_cache.read(reinterpret_cast<char*>(landmark_cached.data()),row_size);
//                auto distPtr = dist.begin();
//                auto fmPtr = first_move.begin();
//                for(unsigned int & startPtr : landmark_cached){
//                    unsigned current_distance = distance + startPtr;
//                    if (current_distance < *(distPtr)) {
//                        *distPtr = current_distance;
//                        *fmPtr = fm;
//                    }
//                    distPtr++;
//                    fmPtr++;
//                }
//
//
//            };
////
            auto updates_fm_and_distance= [&](const unsigned & landmark_id,const unsigned distance, unsigned fm) {
                number_of_landmark_hit++;
                unsigned* contents = (*landmark_pointer)[landmark_id];
                auto distPtr = dist.begin();
                auto fmPtr = first_move.begin();
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
            number_of_reached_landmarks = 0 ;
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            number_of_landmark_hit = 0;
            std::fill(dist.begin(), dist.end(), std::numeric_limits<unsigned>::max());
//            memset(&dist[0],-1,dist.size()*4);
            std::fill( first_move.begin(), first_move.end(), 0xFF);

            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_landmark_sequential((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
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
                                            first_move[x.id]);
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(current_rank   > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_landmark_sequential((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach_landmark_sequential((*head)[a], x.key + (*weight)[a], first_move[x.id],current_rank  < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;

        }



        const std::vector<unsigned short>& run_with_reached_landmark_sequential_with_disk_cache_with_up_down(unsigned source_node){
            auto updates_fm_and_distance= [&](const unsigned & landmark_id,const unsigned distance, unsigned fm) {
                number_of_landmark_hit++;
                unsigned* contents = (*landmark_pointer)[landmark_id];
                auto distPtr = dist.begin();
                auto fmPtr = first_move.begin();
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
            number_of_reached_landmarks = 0 ;
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            number_of_landmark_hit = 0;
            std::fill(dist.begin(), dist.end(), std::numeric_limits<unsigned>::max());
//            memset(&dist[0],-1,dist.size()*4);
            std::fill( first_move.begin(), first_move.end(), 0xFF);

            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;

            for(unsigned a=(*up_first_out)[source_node]; a<(*up_first_out)[source_node+1]; ++a){
                reach_landmark_sequential((*up_head)[a], (*up_weight)[a], fm_index++,2);
                fm_index++;
            }
            for(unsigned a=(*down_first_out)[source_node]; a<(*down_first_out)[source_node+1]; ++a){
                reach_landmark_sequential((*down_head)[a], (*down_weight)[a], fm_index++,1);
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
                                            first_move[x.id]);
                    continue;
                }
                number_of_nodes_expanded ++;

                if(up_or_down[x.id] == 1){
                    // mask 01 indicate moving down
                    for(unsigned a=(*down_first_out)[x.id]; a<(*down_first_out)[x.id+1]; ++a){
                        reach_landmark_sequential((*down_head)[a], x.key + (*down_weight)[a], first_move[x.id],1);
                    }
                }else{
                    for(unsigned a=(*up_first_out)[x.id]; a<(*up_first_out)[x.id+1]; ++a){
                        reach_landmark_sequential((*up_head)[a], x.key + (*up_weight)[a], first_move[x.id],2);
                    }
                    for(unsigned a=(*down_first_out)[x.id]; a<(*down_first_out)[x.id+1]; ++a){
                        reach_landmark_sequential((*down_head)[a], x.key + (*down_weight)[a], first_move[x.id],1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;

        }

        const std::vector<unsigned short>& run_with_reached_landmark_sequential_without_pointer(unsigned source_node){
            auto updates_fm_and_distance= [&](const unsigned & landmark_id,const unsigned distance, unsigned fm) {
//                unsigned* startPtr = (*landmark_pointer)[landmark_id];
                auto  startPtr = (*landmark_list).begin() + (unsigned long long )landmark_id*number_of_nodes;

                auto endPtr = startPtr + number_of_nodes;
                auto distPtr = dist.begin();
                auto fmPtr = first_move.begin();
                for(; startPtr != endPtr ; startPtr++){
                    unsigned current_distance = distance + *(startPtr);
                    if(current_distance  < *(distPtr)){
                        *distPtr = current_distance;
                        *fmPtr = fm;
                    }
                    distPtr++;
                    fmPtr++;
                }
            };

            auto pruned_by_reached_landmarks= [&](const unsigned & current_id,const unsigned& current_distance) {
                for(unsigned i = 0 ; i < number_of_reached_landmarks; i++){
//                    unsigned* startPtr = (*landmark_pointer)[landmark_reached[i]];
                    auto  startPtr = (*landmark_list).begin() + (unsigned long long )landmark_reached[i]*number_of_nodes;
                    if(*(startPtr + source_node) + *(startPtr + current_id) <= current_distance){
                        return true;
                    }
                }
                return false;
            };

            number_of_reached_landmarks = 0 ;
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill(dist.begin(), dist.end(), std::numeric_limits<unsigned>::max());
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_landmark_sequential((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                if(x.key > dist[x.id]){
                    continue;
                }
                unsigned current_rank = (*rank_mapper)[x.id];

                if(current_rank >= landmark_start) {
                    if(!pruned_by_reached_landmarks(x.id,x.key)) {
                        assert(current_rank - landmark_start < (*rank_mapper).size());
//                        updates_fm_and_distance(number_of_landmark - 1 - (current_rank - landmark_start), x.key,
//                                                first_move[x.id]);
//                        landmark_reached[number_of_reached_landmarks] =
//                                number_of_landmark - 1 - (current_rank - landmark_start);
                        updates_fm_and_distance(current_rank - landmark_start, x.key,
                                                first_move[x.id]);
                        landmark_reached[number_of_reached_landmarks] = current_rank - landmark_start;
                        number_of_reached_landmarks++;
                    }
                    continue;
                }


                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(current_rank   > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_landmark_sequential((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach_landmark_sequential((*head)[a], x.key + (*weight)[a], first_move[x.id],current_rank  < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }

        const std::vector<unsigned short>& run_with_reached_landmark(unsigned source_node){
            number_of_reached_landmarks = 0 ;
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_landmark((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                unsigned current_rank = (*rank_mapper)[x.id];
                if(pruned_by_recorded_landmark(source_node,x.id,x.key)){
                    first_move[x.id] =0xFF;
                    continue;
                }
                if(current_rank >= landmark_start) {
//                    landmark_reached[number_of_reached_landmarks] = number_of_landmark- (current_rank - landmark_start) - 1;
                    landmark_reached[number_of_reached_landmarks] = number_of_landmark- (current_rank - landmark_start) - 1;
                    landmark_fm[number_of_reached_landmarks] = first_move[x.id];
                    number_of_reached_landmarks ++;
                    // set to none, it could be reached through other landmark.
                    first_move[x.id] =0xFF;
                    continue;
                }

                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(current_rank   > (*rank_mapper)[(*head)[a]]){

                            // mask 01 indicate moving down
                            reach_landmark((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach_landmark((*head)[a], x.key + (*weight)[a], first_move[x.id],current_rank  < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }

            first_move[source_node] =0xFE;
            unsigned* startPtr = (*landmark_pointer)[source_node];
            for(unsigned i = 0 ; i < first_move.size();i++){
                if(first_move[i] == 0xFF){
                    unsigned* endPtr = (*landmark_pointer)[i];
                    unsigned min_distance = *(startPtr+landmark_reached[0]) +  *(endPtr+landmark_reached[0]);
                    first_move[i] = landmark_fm[0];
                    for(unsigned j = 1; j < number_of_reached_landmarks;j++){
                        unsigned distance  = *(startPtr+landmark_reached[j]) +  *(endPtr+landmark_reached[j]);
                        if(min_distance > distance ){
                            min_distance = distance;
                            first_move[i] = landmark_fm[j];
                        }
                    }
                }
            }
            return first_move;
        }


        const std::vector<unsigned short>& run_with_landmark_without_mapper(unsigned source_node){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( first_move.begin(), first_move.end(), 0xFF);
            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach((*head)[a], (*weight)[a], fm_index,source_node < (*head)[a] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                if(pruned_by_landmark(source_node,x.id,x.key)){
                    first_move[x.id] =0xFF;
                    continue;
                }
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(x.id > (*head)[a]){
                            // mask 01 indicate moving down
                            reach((*head)[a], x.key + (*weight)[a], first_move[x.id],1);
                        }
                    }else{
                        reach((*head)[a], x.key + (*weight)[a], first_move[x.id],x.id < (*head)[a] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }

        const std::vector<unsigned >& get_distance_table(){
            return dist;
        }

        const std::vector<unsigned short>& get_fm_and_distance(unsigned source_node){
            std::fill( first_move.begin(), first_move.end(),0xFF);
//            std::fill( dist.begin(), dist.end(),std::numeric_limits<unsigned>::max());
            memset(&dist[0],-1,dist.size()*4);

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove, const unsigned char up_down_signal) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d < dist[next_vertex]){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
                        dist[next_vertex] = d;
                        first_move[next_vertex] = fmove;

                    }else if(d == dist[next_vertex]){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    first_move[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_nodes((*head)[a], x.key + (*weight)[a],first_move[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a], first_move[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }



        const std::vector<unsigned short>& get_fm_and_distance_with_up_down(unsigned source_node){
            std::fill( first_move.begin(), first_move.end(),0xFF);
//            std::fill( dist.begin(), dist.end(),std::numeric_limits<unsigned>::max());
            memset(&dist[0],-1,dist.size()*4);

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove, const unsigned char up_down_signal) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d < dist[next_vertex]){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
                        dist[next_vertex] = d;
                        first_move[next_vertex] = fmove;

                    }else if(d == dist[next_vertex]){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    first_move[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*up_first_out)[source_node]; a<(*up_first_out)[source_node+1]; ++a){
                reach_nodes((*up_head)[a], (*up_weight)[a], fm_index++,2);
                fm_index++;
            }
            for(unsigned a=(*down_first_out)[source_node]; a<(*down_first_out)[source_node+1]; ++a){
                reach_nodes((*down_head)[a], (*down_weight)[a], fm_index++,1);
                fm_index++;
            }


            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                if(up_or_down[x.id] == 1){
                    // mask 01 indicate moving down
                    for(unsigned a=(*down_first_out)[x.id]; a<(*down_first_out)[x.id+1]; ++a){
                        reach_nodes((*down_head)[a], x.key + (*down_weight)[a], first_move[x.id],1);
                    }
                }else{
                    for(unsigned a=(*up_first_out)[x.id]; a<(*up_first_out)[x.id+1]; ++a){
                        reach_nodes((*up_head)[a], x.key + (*up_weight)[a], first_move[x.id],2);
                    }
                    for(unsigned a=(*down_first_out)[x.id]; a<(*down_first_out)[x.id+1]; ++a){
                        reach_nodes((*down_head)[a], x.key + (*down_weight)[a], first_move[x.id],1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }


        const std::vector<unsigned>& get_fm_and_distance(unsigned source_node,vector<unsigned short>&first_move_list){
            std::fill( first_move_list.begin(), first_move_list.end(),0xFF);
            std::fill( dist.begin(), dist.end(),std::numeric_limits<unsigned>::max());
            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove, const unsigned char up_down_signal) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d < dist[next_vertex]){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
                        dist[next_vertex] = d;
                        first_move_list[next_vertex] = fmove;

                    }else if(d == dist[next_vertex]){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    first_move_list[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_nodes((*head)[a], x.key + (*weight)[a],first_move_list[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a], first_move_list[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move_list[source_node] =0xFE;
            return dist;
        }


        const vector<unsigned short>& get_fm_and_distance(unsigned source_node,std::vector<unsigned>::iterator startPtr){
            std::fill( first_move.begin(), first_move.end(),0xFF);
//            std::fill( startPtr, endPtr,inf_weight);
            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove, const unsigned char up_down_signal) {
                unsigned& current_d = *(startPtr+next_vertex);
                if(was_distance_set.is_set(next_vertex)){

                    if(d < current_d){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
//                        dist[next_vertex]
                        current_d = d;
                        first_move[next_vertex] = fmove;

                    }else if(d == current_d){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    current_d = d;
                    was_distance_set.set(next_vertex);
                    first_move[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            *(startPtr+source_node) = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_nodes((*head)[a], x.key + (*weight)[a],first_move[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a], first_move[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }



        const vector<unsigned short>& get_fm_and_distance(unsigned source_node,unsigned* startPtr){
            std::fill( first_move.begin(), first_move.end(),0xFF);
//            std::fill( startPtr, endPtr,inf_weight);
            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove, const unsigned char up_down_signal) {
                unsigned& current_d = *(startPtr+next_vertex);
                if(was_distance_set.is_set(next_vertex)){
                    if(d < current_d){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
//                        dist[next_vertex]
                        current_d = d;
                        first_move[next_vertex] = fmove;

                    }else if(d == current_d){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    current_d = d;
                    was_distance_set.set(next_vertex);
                    first_move[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            *(startPtr+source_node) = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_nodes((*head)[a], x.key + (*weight)[a],first_move[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a], first_move[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move[source_node] =0xFE;
            return first_move;
        }

        const std::vector<unsigned>& get_predecessor_and_distance(unsigned source_node, const vector<unsigned short>& predecessor,vector<unsigned short>& backward_first_move,vector<unsigned short>&first_move_list){
            std::fill( backward_first_move.begin(), backward_first_move.end(),0xFF);
            std::fill( first_move_list.begin(), first_move_list.end(),0xFF);
            std::fill( dist.begin(), dist.end(),std::numeric_limits<unsigned>::max());

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& predecessor ,const unsigned short& fmove, const unsigned char up_down_signal) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d < dist[next_vertex]){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
                        dist[next_vertex] = d;
                        backward_first_move[next_vertex] =  predecessor;
                        first_move_list[next_vertex] = fmove;

                    }else if(d == dist[next_vertex]){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    backward_first_move[next_vertex] =  predecessor;
                    first_move_list[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a],predecessor[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if((*rank_mapper)[x.id]  > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            reach_nodes((*head)[a], x.key + (*weight)[a],predecessor[a], first_move_list[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a],predecessor[a], first_move_list[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            backward_first_move[source_node] =0xFE;
            first_move_list[source_node] =0xFE;
            return dist;
        }


        const std::vector<unsigned>& get_predecessor_and_distance_without_mapper(unsigned source_node, const vector<unsigned short>& predecessor,vector<unsigned short>& backward_first_move,vector<unsigned short>&first_move_list){
            std::fill( backward_first_move.begin(), backward_first_move.end(),0xFF);
            std::fill( first_move_list.begin(), first_move_list.end(),0xFF);
            std::fill( dist.begin(), dist.end(),std::numeric_limits<unsigned>::max());

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& predecessor ,const unsigned short& fmove, const unsigned char up_down_signal) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d < dist[next_vertex]){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
                        dist[next_vertex] = d;
                        backward_first_move[next_vertex] =  predecessor;
                        first_move_list[next_vertex] = fmove;

                    }else if(d == dist[next_vertex]){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    was_distance_set.set(next_vertex);
                    backward_first_move[next_vertex] =  predecessor;
                    first_move_list[next_vertex] = fmove;
                }
            };

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a],predecessor[a], fm_index,source_node <(*head)[a] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
//                    unsigned cur_weight = (*weight)[a];
                    if(up_or_down[x.id] == 1){
                        if(x.id >(*head)[a]){
                            // mask 01 indicate moving down
                            reach_nodes((*head)[a], x.key + (*weight)[a],predecessor[a], first_move_list[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a],predecessor[a],first_move_list[x.id],x.id < (*head)[a] ? 2:1);
                    }

                }
            }
            backward_first_move[source_node] =0xFE;
            first_move_list[source_node] =0xFE;
            return dist;
        }

    private:
        std::vector<unsigned>dist;
        std::vector<unsigned short>first_move;
        std::vector<unsigned>up_or_down;
        std::vector<unsigned short>landmark_fm;
        std::vector<unsigned short>landmark_reached;
        std::vector<unsigned>landmark_cached;
        unsigned number_of_reached_landmarks;
//        const void* memmap;

        MinIDQueue queue;
        TimestampFlags was_distance_set;
        const std::vector<unsigned>*first_out;
        const std::vector<unsigned>*head;
        const std::vector<unsigned>*weight;
        const std::vector<unsigned>*rank_mapper;
        const std::vector<unsigned*>*landmark_pointer;
        const std::vector<unsigned>*landmark_list;

        const std::vector<unsigned>*up_first_out;
        const std::vector<unsigned>*up_head;
        const std::vector<unsigned>*up_weight;

        const std::vector<unsigned>*down_first_out;
        const std::vector<unsigned>*down_head;
        const std::vector<unsigned>*down_weight;
    };


}
#endif //ROUTINGKIT_CH_DIJKSTRA_SINGLE_H
