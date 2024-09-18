//
// Created by Bojie Shen on 8/11/20.
//

#ifndef ROUTINGKIT_CH_DIJKSTRA_HL_H
#define ROUTINGKIT_CH_DIJKSTRA_HL_H
#include <routingkit/id_queue.h>
#include <routingkit/constants.h>
#include <routingkit/timestamp_flag.h>
#include <vector>
#include <iostream>

using namespace std;
namespace RoutingKit{

    class Ch_Dijkstra_HL{
    public:
        std::vector<bool>visited;
        std::vector<unsigned>updated;
        std::vector<unsigned>updated_distance;
        unsigned number_of_nodes_expanded;
        unsigned number_of_nodes_generated;
        unsigned number_of_nodes;
        unsigned number_of_landmark;
        unsigned landmark_start;
        unsigned update_number;

        Ch_Dijkstra_HL(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight, const std::vector<unsigned>& rank_mapper):
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

        Ch_Dijkstra_HL(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper, const std::vector<unsigned*>&landmark_pointer,
                           unsigned landmark_number,unsigned l_start):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),
                rank_mapper(&rank_mapper),
                landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            visited.resize(first_out.size()-1);
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
            updated = vector<unsigned >(first_out.size()-1,0);
            update_number =0;
        }


        Ch_Dijkstra_HL(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper, const std::vector<unsigned>&landmark_begin,
                       const std::vector<unsigned>&landmark_list, unsigned landmark_number,unsigned l_start):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),
                rank_mapper(&rank_mapper),
                landmark_list(&landmark_list),landmark_begin(&landmark_begin){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            visited.resize(first_out.size()-1);
            first_move.resize(first_out.size()-1);
            up_or_down.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
            updated = vector<unsigned >(first_out.size()-1,0);
            update_number =0;
            updated_distance.resize(first_out.size()-1);
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

       void updates_fm_and_distance_reached(const unsigned & landmark_id,const unsigned distance, unsigned fm) {
            unsigned* startPtr = (*landmark_pointer)[landmark_id];
            unsigned* endPtr = (*landmark_pointer)[landmark_id+1];
            update_number =0;
            vector<unsigned > reached_landmarks;
            while(startPtr != endPtr){
                unsigned index = *startPtr;
                startPtr++;
                unsigned current_distance = distance + *(startPtr);
                startPtr++;
                if(current_distance  < dist[index]){
                    update_number++;
                    dist[index] = current_distance;
                    first_move[index] = fm;
                    if((*rank_mapper)[index] >= landmark_start){
                        reached_landmarks.push_back(index);
//                        reached_landmarks.push_back(number_of_landmark - 1 - ((*rank_mapper)[index] - landmark_start));
                    }
                }


            }
            if(update_number == 0){
                std::cout<<"err"<<std::endl;
            }
            for(unsigned& l : reached_landmarks){
                updates_fm_and_distance_reached(number_of_landmark - 1 - ((*rank_mapper)[l] - landmark_start),dist[l],fm);

            }
        };

        const std::vector<unsigned >& get_distance_table(){
            return dist;
        }
        const std::vector<unsigned short>& run_with_reached_landmark_sequential(unsigned source_node){
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
                    update_number ++;
                    updates_fm_and_distance_reached(number_of_landmark - 1 - (current_rank - landmark_start), x.key,
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

        void reach(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned char &up_down_signal){
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


        void updates_fm_and_distance(const unsigned & landmark_id,const unsigned distance, unsigned fm,vector<unsigned short>&first_move_list) {
            unsigned i = (*landmark_begin)[landmark_id];
            while(i < (*landmark_begin)[landmark_id+1]){
                unsigned index  = (*landmark_list)[i];
                i++;
                unsigned current_distance = distance +  (*landmark_list)[i];
                i++;
                if(current_distance  < dist[index]){
                    updated[index] = update_number;
                    updated_distance[index] = current_distance;
                    dist[index] = current_distance;
                    first_move_list[index] = fm;
                    if((*rank_mapper)[index] >= landmark_start){
                        updates_fm_and_distance(number_of_landmark - 1 - ((*rank_mapper)[index] - landmark_start),current_distance ,fm, first_move_list);
                    }
                }
            }
        };


        const std::vector<unsigned>& build_landmark(unsigned source_node,vector<unsigned short>&first_move_list){
            auto build_reached= [&](const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned char &up_down_signal) {

                if(d < dist[next_vertex]){
                    number_of_nodes_generated++;
                    if(queue.contains_id(next_vertex)){
                        queue.decrease_key({next_vertex, d});
                    }else{
                        queue.push({next_vertex, d});
                    }
                    up_or_down[next_vertex] = up_down_signal;
                    dist[next_vertex] = d;
                    first_move_list[next_vertex] = fmove;
                }else if(d == dist[next_vertex]){
                    up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                }
            };


            number_of_reached_landmarks = 0 ;
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            update_number++;
            std::fill(dist.begin(), dist.end(), std::numeric_limits<unsigned>::max());
            std::fill(visited.begin(),visited.end(), false);
            std::fill( first_move_list.begin(), first_move_list.end(), 0xFF);
            queue.clear();
            dist[source_node] = 0;
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                build_reached((*head)[a], (*weight)[a], fm_index,(*rank_mapper)[source_node] < (*rank_mapper)[(*head)[a]] ? 2:1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                assert(x.id<dist.size());
                if(x.key > dist[x.id]){
                    continue;
                }else{
                   if(x.key ==updated_distance[x.id] && updated[x.id] == update_number){
                       continue;
                   }
                }
                unsigned current_rank = (*rank_mapper)[x.id];
                visited[x.id] = true;
                if(current_rank >= landmark_start) {
                    updates_fm_and_distance(number_of_landmark - 1 - (current_rank - landmark_start), x.key,
                                            first_move_list[x.id],first_move_list);
                    continue;
                }


                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    if(up_or_down[x.id] == 1){
                        if(current_rank   > (*rank_mapper)[(*head)[a]]){
                            // mask 01 indicate moving down
                            build_reached((*head)[a], x.key + (*weight)[a], first_move_list[x.id],1);
                        }
                    }else{
                        build_reached((*head)[a], x.key + (*weight)[a], first_move_list[x.id],current_rank  < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move_list[source_node] =0xFE;
            return dist;
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


        void get_fm_and_distance(unsigned source_node,std::vector<unsigned>::iterator startPtr,vector<unsigned short>&first_move_list){
            std::fill( first_move_list.begin(), first_move_list.end(),0xFF);

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& fmove, const unsigned char up_down_signal) {
                unsigned& current_d = *(startPtr+next_vertex);
                if(was_distance_set.is_set(next_vertex)){

                    if(d < current_d){
                        number_of_nodes_generated++;
                        queue.decrease_key({next_vertex, d});
                        up_or_down[next_vertex] = up_down_signal;
//                        dist[next_vertex]
                        current_d = d;
                        first_move_list[next_vertex] = fmove;

                    }else if(d == current_d){
                        up_or_down[next_vertex] = up_or_down[next_vertex] | up_down_signal;
                    }
                }else{
                    number_of_nodes_generated++;
                    queue.push({next_vertex,d});
                    up_or_down[next_vertex] = up_down_signal;
                    current_d = d;
                    was_distance_set.set(next_vertex);
                    first_move_list[next_vertex] = fmove;
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
                            reach_nodes((*head)[a], x.key + (*weight)[a],first_move_list[x.id],1);
                        }
                    }else{
                        reach_nodes((*head)[a], x.key + (*weight)[a], first_move_list[x.id],(*rank_mapper)[x.id] < (*rank_mapper)[(*head)[a]] ? 2:1);
                    }

                }
            }
            first_move_list[source_node] =0xFE;
        }



    private:

        std::vector<unsigned>dist;
        std::vector<unsigned short>first_move;
        std::vector<unsigned>up_or_down;
        unsigned number_of_reached_landmarks;


        MinIDQueue queue;
        TimestampFlags was_distance_set;
        const std::vector<unsigned>*first_out;
        const std::vector<unsigned>*head;
        const std::vector<unsigned>*weight;
        const std::vector<unsigned>*rank_mapper;
        const std::vector<unsigned*>*landmark_pointer;
        const std::vector<unsigned>*landmark_list;
        const std::vector<unsigned>*landmark_begin;
    };


}
#endif //ROUTINGKIT_CH_DIJKSTRA_HL_H
