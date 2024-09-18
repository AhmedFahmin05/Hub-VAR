//
// Created by Bojie Shen on 1/11/20.
//

#ifndef ROUTINGKIT_SHORTEST_DIJKSTRA_SINGLE_H
#define ROUTINGKIT_SHORTEST_DIJKSTRA_SINGLE_H
#include <routingkit/id_queue.h>
#include <routingkit/constants.h>
#include <routingkit/timestamp_flag.h>
#include <vector>
#include <iostream>

using namespace std;
namespace RoutingKit{

    class Shortest_Dijkstra_Single{
    public:
        unsigned number_of_nodes_expanded;
        unsigned number_of_nodes_generated;
        unsigned number_of_landmark;

        Shortest_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            first_move.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
        }

        Shortest_Dijkstra_Single(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight, const std::vector<unsigned*>&landmark_pointer,
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
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
            number_of_landmark  = landmark_number;
        }


        void reach(const unsigned & next_vertex, unsigned d, const unsigned short& fmove, const unsigned &number_of_vertices){
            if(was_distance_set.is_set(next_vertex)){
                if(d <  dist[next_vertex]){
                    number_of_nodes_generated++;
                    queue.decrease_key({next_vertex, d});
                    parent[next_vertex] = number_of_vertices;
                    dist[next_vertex] = d;
                    first_move[next_vertex] = fmove;
                }else if(d == dist[next_vertex]){
                    if(parent[next_vertex] > number_of_vertices){
                        first_move[next_vertex] = fmove;
                        parent[next_vertex] = number_of_vertices;
                    }
                }
            }else{
                queue.push({next_vertex,d});
                parent[next_vertex] = number_of_vertices;
                dist[next_vertex] = d;
                was_distance_set.set(next_vertex);
                first_move[next_vertex] = fmove;
                number_of_nodes_generated ++;
            }

        }
        const std::vector<unsigned short>& run(int source_node){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( first_move.begin(), first_move.end(),0xFF);

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach((*head)[a], (*weight)[a], fm_index,1);
                fm_index++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    reach((*head)[a], x.key + (*weight)[a], first_move[x.id],parent[x.id]+1);
                }
            }
            first_move[source_node] =0xFE;

            return first_move;
        }

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

        const std::vector<unsigned short>& run_with_landmark(int source_node){
            number_of_nodes_expanded = 0;
            number_of_nodes_generated = 0;
            std::fill( first_move.begin(), first_move.end(),0xFF);

            queue.clear();
            was_distance_set.reset_all();
            dist[source_node] = 0;
            was_distance_set.set(source_node);
            unsigned short fm_index = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach((*head)[a], (*weight)[a], fm_index,1);
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
                    reach((*head)[a], x.key + (*weight)[a], first_move[x.id],parent[x.id]+1);
                }
            }
            first_move[source_node] =0xFE;

            return first_move;
        }




        const std::vector<unsigned>& get_predecessor_and_distance(int source_node, const vector<unsigned short>& predecessor,vector<unsigned short>& backward_first_move,vector<unsigned short>&first_move_list){
            std::fill( backward_first_move.begin(), backward_first_move.end(),0xFF);
            std::fill( first_move_list.begin(), first_move_list.end(),0xFF);
            std::fill( dist.begin(), dist.end(),inf_weight);

            auto reach_nodes = [&](const unsigned & next_vertex, unsigned d,const unsigned short& predecessor ,const unsigned short& fmove, const unsigned &number_of_vertices) {
                if(was_distance_set.is_set(next_vertex)){
                    if(d <  dist[next_vertex]){
                        queue.decrease_key({next_vertex, d});
                        parent[next_vertex] = number_of_vertices;
                        dist[next_vertex] = d;
                        backward_first_move[next_vertex] =  predecessor;
                        first_move_list[next_vertex] = fmove;
                    }else if(d == dist[next_vertex]){
                        if(parent[next_vertex] > number_of_vertices){
                            backward_first_move[next_vertex] =  predecessor;
                            first_move_list[next_vertex] = fmove;
                            parent[next_vertex] = number_of_vertices;
                        }
                    }
                }else{
                    queue.push({next_vertex,d});
                    parent[next_vertex] = number_of_vertices;
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
            unsigned short fm = 0;
            for(unsigned a=(*first_out)[source_node]; a<(*first_out)[source_node+1]; ++a){
                reach_nodes((*head)[a], (*weight)[a], predecessor[a],fm,1);
                fm++;
            }

            while(!queue.empty()){
                auto x = queue.pop();
                number_of_nodes_expanded ++;
                for(unsigned a=(*first_out)[x.id]; a<(*first_out)[x.id+1]; ++a){
                    reach_nodes((*head)[a], x.key + (*weight)[a], predecessor[a],first_move_list[x.id],parent[x.id]+1);
                }
            }
            backward_first_move[source_node] =0xFE;
            first_move_list[source_node] =0xFE;
            return dist;
        }


    private:
        std::vector<unsigned>dist;
        std::vector<unsigned short>first_move;
        std::vector<unsigned>parent;
        MinIDQueue queue;
        TimestampFlags was_distance_set;
        unsigned number_of_nodes;
        const std::vector<unsigned>*first_out;
        const std::vector<unsigned>*head;
        const std::vector<unsigned>*weight;
        const std::vector<unsigned*>*landmark_pointer;











    };








}









#endif //ROUTINGKIT_SHORTEST_DIJKSTRA_SINGLE_H
