//
// Created by Bojie Shen on 7/2/21.
//

#ifndef ROUTINGKIT_ORI_DIJKSTRA_H
#define ROUTINGKIT_ORI_DIJKSTRA_H

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
#include <fstream>

namespace RoutingKit{
    class Ori_Dijkstra{
    public:
        unsigned number_of_nodes_expanded;
        unsigned number_of_nodes_generated;
        unsigned number_of_nodes;
        unsigned number_of_landmark;
        unsigned number_of_landmark_hit;
        unsigned landmark_start;

        Ori_Dijkstra():first_out(nullptr){}


        Ori_Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head,const std::vector<unsigned>&weight):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
        }

        Ori_Dijkstra(const std::vector<unsigned>&first_out, const std::vector<unsigned>&head, const std::vector<unsigned>&weight,const std::vector<unsigned>& rank_mapper,
                     unsigned landmark_number,unsigned l_start,const std::vector<unsigned*>& landmark_pointer):
                queue(first_out.size()-1),
                was_distance_set(first_out.size()-1),
                first_out(&first_out),
                head(&head),
                weight(&weight),rank_mapper(&rank_mapper),landmark_pointer(&landmark_pointer){
            assert(!first_out.empty());
            assert(first_out.front() == 0);
            assert(first_out.back() == head.size());
            d_allowed.resize(first_out.size()-1);
            parent.resize(first_out.size()-1);
            dist.resize(first_out.size()-1);
            number_of_nodes = first_out.size() -1;
            number_of_landmark  = landmark_number;
            landmark_start = l_start;
            number_of_nodes = first_out.size()-1;
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

        std::vector<unsigned>dist;
        std::vector<unsigned short>d_allowed;
        std::vector<unsigned>parent;
        MinIDQueue queue;
        TimestampFlags was_distance_set;

        const std::vector<unsigned>*first_out;
        const std::vector<unsigned>*head;
        const std::vector<unsigned>*weight;
        const std::vector<unsigned>*rank_mapper;
        const std::vector<unsigned*>*landmark_pointer;
    };


}

#endif //ROUTINGKIT_ORI_DIJKSTRA_H
