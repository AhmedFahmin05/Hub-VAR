#include <routingkit/contraction_hierarchy.h>
#include <routingkit/id_queue.h>
#include <routingkit/sort.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/timer.h>
#include <routingkit/graph_util.h>
#include <routingkit/vector_io.h>

#include <vector>
#include <fstream>
#include <stdexcept>
#include <routingkit/geo_dist.h>
#include "CPD.h"

namespace RoutingKit{

    namespace{

        void sort_arcs_and_remove_multi_and_loop_arcs(
                unsigned node_count, std::vector<unsigned>&tail, std::vector<unsigned>&head, std::vector<unsigned>&weight, std::vector<unsigned>&input_arc_id,
                const std::function<void(std::string)>&log_message
        ){

            long long timer = 0;  // initialize to avoid warning, not needed
            if(log_message){
                timer = -get_micro_time();
                log_message("Start removing loops and multi arcs from input.");
            }

            {
                auto p = compute_inverse_sort_permutation_first_by_tail_then_by_head_and_apply_sort_to_tail(node_count, tail, head);
                head = apply_inverse_permutation(p, std::move(head));
                weight = apply_inverse_permutation(p, std::move(weight));
                input_arc_id = apply_inverse_permutation(p, std::move(input_arc_id));
            }

            unsigned arc_count = head.size();

            if(arc_count != 0){

                unsigned out = 0;
                for(unsigned in = 0; in < arc_count; ++in){
                    if(tail[in] != head[in]){
                        tail[out] = tail[in];
                        head[out] = head[in];
                        weight[out] = weight[in];
                        input_arc_id[out] = input_arc_id[in];
                        ++out;
                    }
                }
                arc_count = out;
            }

            if(arc_count != 0)
            {
                unsigned out = 1;
                for(unsigned in = 1; in < arc_count; ++in){
                    if(tail[in-1] != tail[in] || head[in-1] != head[in]){
                        tail[out] = tail[in];
                        head[out] = head[in];
                        weight[out] = weight[in];
                        input_arc_id[out] = input_arc_id[in];
                        ++out;
                    }else{
                        if(weight[in] < weight[out-1]){
                            weight[out-1] = weight[in];
                            input_arc_id[out-1] = input_arc_id[in];
                        }
                    }
                }
                arc_count = out;
            }

            tail.erase(tail.begin()+arc_count, tail.end());
            head.erase(head.begin()+arc_count, head.end());
            weight.erase(weight.begin()+arc_count, weight.end());
            input_arc_id.erase(input_arc_id.begin()+arc_count, input_arc_id.end());

            if(log_message){
                timer += get_micro_time();
                log_message("The input map size: "+ to_string(head.size()));
                log_message("Finished removing loops and multi arcs from input. Needed "+std::to_string(timer)+"musec time.");
            }
        }




        class ShorterPathTest{
        public:
            ShorterPathTest(){}
            ShorterPathTest(const Graph&graph, unsigned max_pop_count):
                    max_pop_count(max_pop_count), graph(&graph),
                    forward_tentative_distance(graph.node_count()), backward_tentative_distance(graph.node_count()),
                    forward_queue(graph.node_count()), backward_queue(graph.node_count()),
                    was_forward_pushed(graph.node_count()), was_backward_pushed(graph.node_count())
            {}

        private:
// Uncomment the following lines to get some very expensive but very detailed asserts
            void assert_no_witness_found(unsigned len)const{
//			for(unsigned x = 0; x<graph->node_count(); ++x){
//				if(was_forward_pushed.is_set(x) && was_backward_pushed.is_set(x)){
//					assert(forward_tentative_distance[x] + backward_tentative_distance[x] > len);
//				}
//			}
            }

            void assert_forward_tentative_distances_correct()const{
//			for(unsigned x = 0; x<graph->node_count(); ++x){
//				if(was_forward_pushed.is_set(x)){
//					if(forward_tentative_distance[x] != 0){
//						bool witness_found = false;
//						for(unsigned yx = 0; yx < graph->in_deg(x); ++yx){
//							unsigned y = graph->in(x, yx).node;
//							unsigned w = graph->in(x, yx).weight;
//							if(forward_tentative_distance[y] + w == forward_tentative_distance[x])
//								witness_found = true;
//						}
//						assert(witness_found && "forward_tentative_distance array contains garbage");
//					}
//				}
//			}
            }

            void assert_backward_tentative_distances_correct()const{
//			for(unsigned x = 0; x<graph->node_count(); ++x){
//				if(was_backward_pushed.is_set(x)){
//					if(backward_tentative_distance[x] != 0){
//						bool witness_found = false;
//						for(unsigned xy = 0; xy < graph->out_deg(x); ++xy){
//							unsigned y = graph->out(x, xy).node;
//							unsigned w = graph->out(x, xy).weight;
//							if(backward_tentative_distance[y] + w == backward_tentative_distance[x])
//								witness_found = true;
//						}
//						assert(witness_found && "backward_tentative_distance array contains garbage");
//					}
//				}
//			}
            }

            void assert_witness_found(unsigned len)const{
//			assert_forward_tentative_distances_correct();
//			assert_backward_tentative_distances_correct();
//			for(unsigned x = 0; x<graph->node_count(); ++x){
//				if(was_forward_pushed.is_set(x) && was_backward_pushed.is_set(x)){
//					if(forward_tentative_distance[x] + backward_tentative_distance[x] <= len)
//						return;
//				}
//			}
//			assert(false && "no witness exists but algorithm claimed that there is one");
            }

            template<class GetOutDeg, class GetOutArc>
            bool forward_settle(
                    MinIDQueue&forward_queue,
                    TimestampFlags&was_forward_pushed,
                    const TimestampFlags&was_backward_pushed,
                    std::vector<unsigned>&forward_tentative_distance,
                    const std::vector<unsigned>&backward_tentative_distance,
                    const GetOutDeg&graph_out_deg,
                    const GetOutArc&graph_out,
                    unsigned bypass,
                    unsigned len
            ){
                auto p = forward_queue.pop();

                unsigned popped_node = p.id;
                unsigned distance_to_popped_node = p.key;

                assert(forward_tentative_distance[popped_node] == distance_to_popped_node);
                assert(was_forward_pushed.is_set(popped_node));

                if(was_backward_pushed.is_set(popped_node)){
                    if(distance_to_popped_node + backward_tentative_distance[popped_node] <= len){
                        assert_witness_found(len);
                        return true;
                    }
                }

                bool witness_found = false;

                for(unsigned out_arc = 0; out_arc < graph_out_deg(popped_node); ++out_arc){
                    unsigned next_node = graph_out(popped_node, out_arc).node;

                    if(next_node == bypass)
                        continue;

                    unsigned next_node_distance = distance_to_popped_node + graph_out(popped_node, out_arc).weight;

                    if(was_forward_pushed.is_set(next_node)){
                        if(next_node_distance < forward_tentative_distance[next_node]){
                            forward_queue.decrease_key({next_node, next_node_distance});
                            forward_tentative_distance[next_node] = next_node_distance;

                            if(was_backward_pushed.is_set(next_node)){
                                if(next_node_distance + backward_tentative_distance[next_node] <= len){
                                    assert_witness_found(len);
                                    witness_found = true;
                                }
                            }
                        }
                    } else {
                        was_forward_pushed.set(next_node);
                        forward_tentative_distance[next_node] = next_node_distance;
                        forward_queue.push({next_node, next_node_distance});

                        if(was_backward_pushed.is_set(next_node)){
                            if(next_node_distance + backward_tentative_distance[next_node] <= len){
                                assert_witness_found(len);
                                witness_found = true;
                            }
                        }

                    }
                }
                if(!witness_found)
                    assert_no_witness_found(len);
                return witness_found;
            }

            unsigned bypass_node;
        public:
            void pin_source(unsigned s, unsigned new_bypass_node){
                was_forward_pushed.reset_all();
                forward_queue.clear();
                forward_queue.push({s, 0});
                forward_tentative_distance[s] = 0;
                was_forward_pushed.set(s);
                bypass_node = new_bypass_node;
            }

            bool does_shorter_or_equal_path_to_target_exist(unsigned t, unsigned len){
                was_backward_pushed.reset_all();
                backward_queue.clear();
                backward_queue.push({t, 0});
                backward_tentative_distance[t] = 0;
                was_backward_pushed.set(t);

                unsigned pop_count = 0;

                if(was_forward_pushed.is_set(t))
                    if(forward_tentative_distance[t] <= len)
                        return true;

                assert_no_witness_found(len);

                while(!forward_queue.empty() && !backward_queue.empty()){
                    if(forward_queue.peek().key + backward_queue.peek().key > len){
                        return false;
                    }

                    if(forward_queue.peek().key <= backward_queue.peek().key){
                        if(
                                forward_settle(
                                        forward_queue,
                                        was_forward_pushed, was_backward_pushed,
                                        forward_tentative_distance, backward_tentative_distance,
                                        [&](unsigned x){return graph->out_deg(x);},
                                        [&](unsigned x, unsigned a){return graph->out(x, a);},
                                        bypass_node,
                                        len
                                )
                                ){
                            assert_witness_found(len);
                            return true;
                        } else {
                            assert_no_witness_found(len);
                        }
                    } else {
                        if(
                                forward_settle(
                                        backward_queue,
                                        was_backward_pushed, was_forward_pushed,
                                        backward_tentative_distance, forward_tentative_distance,
                                        [&](unsigned x){return graph->in_deg(x);},
                                        [&](unsigned x, unsigned a){return graph->in(x, a);},
                                        bypass_node,
                                        len
                                )
                                ){
                            assert_witness_found(len);
                            return true;
                        } else {
                            assert_no_witness_found(len);
                        }
                    }

                    ++pop_count;

                    if(pop_count > max_pop_count){
                        return false;
                    }
                }
                return false;
            }

            bool does_shorter_or_equal_path_exist(unsigned s, unsigned t, unsigned len, unsigned bypass){
                if(s == t)
                    return true;

                was_forward_pushed.reset_all();
                was_backward_pushed.reset_all();

                forward_queue.clear();
                backward_queue.clear();

                forward_queue.push({s, 0});
                backward_queue.push({t, 0});

                forward_tentative_distance[s] = 0;
                backward_tentative_distance[t] = 0;

                was_forward_pushed.set(s);
                was_backward_pushed.set(t);

                unsigned pop_count = 0;

                assert_no_witness_found(len);

                while(!forward_queue.empty() && !backward_queue.empty()){

                    if(forward_queue.peek().key + backward_queue.peek().key > len){
                        return false;
                    }

                    if(forward_queue.peek().key <= backward_queue.peek().key){
                        if(
                                forward_settle(
                                        forward_queue,
                                        was_forward_pushed, was_backward_pushed,
                                        forward_tentative_distance, backward_tentative_distance,
                                        [&](unsigned x){return graph->out_deg(x);},
                                        [&](unsigned x, unsigned a){return graph->out(x, a);},
                                        bypass,
                                        len
                                )
                                ){
                            assert_witness_found(len);
                            return true;
                        } else {
                            assert_no_witness_found(len);
                        }
                    } else {
                        if(
                                forward_settle(
                                        backward_queue,
                                        was_backward_pushed, was_forward_pushed,
                                        backward_tentative_distance, forward_tentative_distance,
                                        [&](unsigned x){return graph->in_deg(x);},
                                        [&](unsigned x, unsigned a){return graph->in(x, a);},
                                        bypass,
                                        len
                                )
                                ){
                            assert_witness_found(len);
                            return true;
                        } else {
                            assert_no_witness_found(len);
                        }
                    }

                    ++pop_count;

                    if(pop_count > max_pop_count)
                        return false;
                }
                return false;
            }

            unsigned get_max_pop_count()const{
                return max_pop_count;
            }
        private:
            unsigned max_pop_count;
            const Graph*graph;
            std::vector<unsigned>forward_tentative_distance;
            std::vector<unsigned>backward_tentative_distance;
            MinIDQueue forward_queue;
            MinIDQueue backward_queue;
            TimestampFlags was_forward_pushed;
            TimestampFlags was_backward_pushed;
        };

        unsigned estimate_node_importance(const Graph&graph, ShorterPathTest&shorter_path_test, unsigned node){
            unsigned level = graph.level(node);

            unsigned added_arc_count = 0;
            unsigned added_hop_count = 0;

            for(unsigned in_arc = 0; in_arc < graph.in_deg(node); ++in_arc){
                unsigned in_node = graph.in(node, in_arc).node;
                shorter_path_test.pin_source(in_node, node);
                for(unsigned out_arc = 0; out_arc < graph.out_deg(node); ++out_arc){
                    unsigned out_node = graph.out(node, out_arc).node;
                    if(in_node != out_node){
                        if(
                                !shorter_path_test.does_shorter_or_equal_path_to_target_exist(
                                        out_node,
                                        graph.in(node, in_arc).weight + graph.out(node, out_arc).weight
                                )
                                ){
                            ++added_arc_count;
                            added_hop_count += graph.in(node, in_arc).hop_length;
                            added_hop_count += graph.out(node, out_arc).hop_length;
                        }
                    }
                }
            }

            unsigned removed_arc_count = 1;
            removed_arc_count += graph.in_deg(node);
            removed_arc_count += graph.out_deg(node);

            unsigned removed_hop_count = 1;
            for(unsigned in_arc = 0; in_arc < graph.in_deg(node); ++in_arc)
                removed_hop_count += graph.in(node, in_arc).hop_length;
            for(unsigned out_arc = 0; out_arc < graph.out_deg(node); ++out_arc)
                removed_hop_count += graph.out(node, out_arc).hop_length;

            return 1 + 1000*level + (1000*added_arc_count) / removed_arc_count + (1000*added_hop_count) / removed_hop_count;
        }

        void contract_node(Graph&graph, ShorterPathTest&shorter_path_test, unsigned node_being_contracted){
            for(unsigned in_arc = 0; in_arc < graph.in_deg(node_being_contracted); ++in_arc){
                unsigned in_node = graph.in(node_being_contracted, in_arc).node;
                shorter_path_test.pin_source(in_node, node_being_contracted);
                for(unsigned out_arc = 0; out_arc < graph.out_deg(node_being_contracted); ++out_arc){
                    unsigned out_node = graph.out(node_being_contracted, out_arc).node;
                    if(in_node != out_node){
                        if(
                                !shorter_path_test.does_shorter_or_equal_path_to_target_exist(
                                        out_node,
                                        graph.in(node_being_contracted, in_arc).weight + graph.out(node_being_contracted, out_arc).weight
                                )
                                ){
                            graph.add_arc_or_reduce_arc_weight(
                                    in_node, node_being_contracted, out_node,
                                    graph.in(node_being_contracted, in_arc).weight     + graph.out(node_being_contracted, out_arc).weight,
                                    graph.in(node_being_contracted, in_arc).hop_length + graph.out(node_being_contracted, out_arc).hop_length
                            );
                        }
                    }
                }
            }

            graph.remove_all_incident_arcs(node_being_contracted);

            assert(graph.out_deg(node_being_contracted) == 0);
            assert(graph.in_deg(node_being_contracted) == 0);
        }

    }

    namespace {

        struct ContractionHierarchyExtraInfo{
            struct Side{
                std::vector<unsigned>mid_node;
                std::vector<unsigned>tail;
            };

            Side forward, backward;
        };

        void build_ch_and_order(
                Graph&graph,
                ContractionHierarchy&ch,
                ContractionHierarchyExtraInfo&ch_extra,
                unsigned max_pop_count,
                const std::function<void(std::string)>&log_message
        ){
            long long timer = 0;  // initialize to avoid warning, not needed
            long long last_log_message_time = 0;  // initialize to avoid warning, not needed
            if(log_message){
                last_log_message_time = get_micro_time();
                timer = -last_log_message_time;
                log_message("Start building queue.");
            }

            const unsigned node_count = graph.node_count();

            ShorterPathTest shorter_path_test(graph, max_pop_count);

            ch.rank.resize(node_count);
            ch.order.resize(node_count);
            MinIDQueue queue(node_count);


            for(unsigned i=0; i<node_count; ++i){
                queue.push({i, estimate_node_importance(graph, shorter_path_test, i)});

                if(log_message){
                    long long current_time = get_micro_time();
                    if(current_time - last_log_message_time > 1000000){
                        last_log_message_time = current_time;
                        log_message("Added "+std::to_string(i+1) + " of " + std::to_string(node_count) + " nodes to the queue. Running for "+std::to_string(timer+current_time)+"musec.");
                    }
                }
            }

            if(log_message){
                timer += get_micro_time();
                log_message("Finished building queue. Needed "+std::to_string(timer)+"musec time.");
                log_message("Start contracting nodes.");
                timer = -get_micro_time();
            }

            std::vector<unsigned>neighbor_list;
            std::vector<bool>is_neighbor(node_count, false);

            unsigned contracted_node_count = 0;

            while(!queue.empty()){
                unsigned node_being_contracted = queue.pop().id;

                ch.rank[node_being_contracted] = contracted_node_count;
                ch.order[contracted_node_count] = node_being_contracted;

                // Mark the neighbors
                for(unsigned in_arc = 0; in_arc < graph.in_deg(node_being_contracted); ++in_arc){
                    unsigned x = graph.in(node_being_contracted, in_arc).node;
                    assert(node_being_contracted != x);
                    if(!is_neighbor[x]){
                        neighbor_list.push_back(x);
                        is_neighbor[x] = true;
                    }
                }

                for(unsigned out_arc = 0; out_arc < graph.out_deg(node_being_contracted); ++out_arc){
                    unsigned x = graph.out(node_being_contracted, out_arc).node;
                    assert(node_being_contracted != x);
                    if(!is_neighbor[x]){
                        neighbor_list.push_back(x);
                        is_neighbor[x] = true;
                    }
                }

                // Add the arcs to the search graph

                for(unsigned out_arc = 0; out_arc < graph.out_deg(node_being_contracted); ++out_arc){
                    ch_extra.forward.tail.push_back(node_being_contracted);

                    const auto&a = graph.out(node_being_contracted, out_arc);
                    if(ch.forward.head.size() == invalid_id)
                        throw std::runtime_error("CH may contain at most 2^32-1 shortcuts per direction");
                    ch.forward.head.push_back(a.node);
                    ch.forward.weight.push_back(a.weight);
                    ch_extra.forward.mid_node.push_back(a.mid_node);
                }

                for(unsigned in_arc = 0; in_arc < graph.in_deg(node_being_contracted); ++in_arc){
                    ch_extra.backward.tail.push_back(node_being_contracted);

                    const auto&a = graph.in(node_being_contracted, in_arc);
                    if(ch.backward.head.size() == invalid_id)
                        throw std::runtime_error("CH may contain at most 2^32-1 shortcuts per direction");
                    ch.backward.head.push_back(a.node);
                    ch.backward.weight.push_back(a.weight);
                    ch_extra.backward.mid_node.push_back(a.mid_node);
                }

                unsigned neighbor_level = graph.level(node_being_contracted)+1;
                unsigned out_deg = graph.out_deg(node_being_contracted);
                unsigned in_deg = graph.in_deg(node_being_contracted);

                contract_node(graph, shorter_path_test, node_being_contracted);

                for(auto x:neighbor_list){
                    is_neighbor[x] = false;
                    graph.raise_level(x, neighbor_level);
                    unsigned new_key = estimate_node_importance(graph, shorter_path_test, x);
                    assert(queue.contains_id(x));
                    unsigned old_key = queue.get_key(x);
                    if(old_key < new_key)
                        queue.increase_key({x, new_key});
                    else if(old_key > new_key)
                        queue.decrease_key({x, new_key});
                }

                neighbor_list.clear();

                ++contracted_node_count;

                if(log_message){

                    long long current_time = get_micro_time();
                    if(current_time - last_log_message_time > 1000000){
                        last_log_message_time = current_time;
                        log_message("Contracted "+std::to_string(contracted_node_count) + " of " + std::to_string(node_count) + ". The in degree of last node was " + std::to_string(in_deg)+ " and out degree was " + std::to_string(out_deg)+". Running for "+std::to_string(timer+current_time)+"musec.");
                    }
                }
            }

            ch.forward.head.shrink_to_fit();
            ch.forward.weight.shrink_to_fit();
            ch_extra.forward.mid_node.shrink_to_fit();
            ch_extra.forward.tail.shrink_to_fit();

            ch.backward.head.shrink_to_fit();
            ch.backward.weight.shrink_to_fit();
            ch_extra.backward.mid_node.shrink_to_fit();
            ch_extra.backward.tail.shrink_to_fit();

            if(log_message){
                timer += get_micro_time();
                log_message("Finished contracting nodes. Needed "+std::to_string(timer)+"musec.");
            }

        }


        void build_ch_given_rank(
                Graph&graph,
                ContractionHierarchy&ch,
                ContractionHierarchyExtraInfo&ch_extra,
                const std::vector<unsigned>&rank,
                unsigned max_pop_count,
                const std::function<void(std::string)>&log_message
        ){
            unsigned node_count = graph.node_count();

            long long timer = 0;  // initialize to avoid warning, not needed
            long long last_log_message_time = 0;  // initialize to avoid warning, not needed
            if(log_message){
                last_log_message_time = get_micro_time();
                timer = -last_log_message_time;
                log_message("Start building contraction hierarchy with given rank.");
            }

            ShorterPathTest shorter_path_test(graph, max_pop_count);
            ch.rank = rank;

            ch.order = invert_permutation(rank);

            for(unsigned i=0; i < node_count; ++i){
                unsigned node_being_contracted = ch.order[i];

                for(unsigned out_arc = 0; out_arc < graph.out_deg(node_being_contracted); ++out_arc){
                    ch_extra.forward.tail.push_back(node_being_contracted);

                    const auto&a = graph.out(node_being_contracted, out_arc);
                    if(ch.forward.head.size() == invalid_id)
                        throw std::runtime_error("CH may contain at most 2^32-1 shortcuts per direction");
                    ch.forward.head.push_back(a.node);
                    ch.forward.weight.push_back(a.weight);
                    ch_extra.forward.mid_node.push_back(a.mid_node);
                }

                for(unsigned in_arc = 0; in_arc < graph.in_deg(node_being_contracted); ++in_arc){
                    ch_extra.backward.tail.push_back(node_being_contracted);

                    const auto&a = graph.in(node_being_contracted, in_arc);
                    if(ch.backward.head.size() == invalid_id)
                        throw std::runtime_error("CH may contain at most 2^32-1 shortcuts per direction");
                    ch.backward.head.push_back(a.node);
                    ch.backward.weight.push_back(a.weight);
                    ch_extra.backward.mid_node.push_back(a.mid_node);
                }

                unsigned out_deg = graph.out_deg(node_being_contracted);
                unsigned in_deg = graph.in_deg(node_being_contracted);

                contract_node(graph, shorter_path_test, node_being_contracted);

                if(log_message){
                    long long current_time = get_micro_time();
                    if(current_time - last_log_message_time > 1000000){
                        last_log_message_time = current_time;
                        log_message("Contracted "+std::to_string(i+1) + " of " + std::to_string(node_count) + ". The in degree of last node was " + std::to_string(in_deg)+ " and out degree was " + std::to_string(out_deg)+". Running for "+std::to_string(timer+current_time)+"musec.");
                    }
                }
            }

            ch.forward.head.shrink_to_fit();
            ch.forward.weight.shrink_to_fit();
            ch_extra.forward.mid_node.shrink_to_fit();
            ch_extra.forward.tail.shrink_to_fit();

            ch.backward.head.shrink_to_fit();
            ch.backward.weight.shrink_to_fit();
            ch_extra.backward.mid_node.shrink_to_fit();
            ch_extra.backward.tail.shrink_to_fit();

            if(log_message){
                timer += get_micro_time();
                log_message("Finished contracting nodes. Needed "+std::to_string(timer)+"musec.");
            }
        }


        void make_internal_nodes_and_rank_coincide(
                ContractionHierarchy&ch,
                ContractionHierarchyExtraInfo&ch_extra,
                const std::function<void(std::string)>&log_message
        ){
            long long timer = 0; // initialize to avoid warning, not needed
            if(log_message){
                timer = -get_micro_time();
                log_message("Start reordering nodes by rank.");
            }

            inplace_apply_permutation_to_elements_of(ch.rank, ch.forward.head);
            inplace_apply_permutation_to_elements_of(ch.rank, ch_extra.forward.tail);
            inplace_apply_permutation_to_possibly_invalid_elements_of(ch.rank, ch_extra.forward.mid_node);
            inplace_apply_permutation_to_elements_of(ch.rank, ch.backward.head);
            inplace_apply_permutation_to_elements_of(ch.rank, ch_extra.backward.tail);
            inplace_apply_permutation_to_possibly_invalid_elements_of(ch.rank, ch_extra.backward.mid_node);

#ifndef NDEBUG
            for(unsigned i=0; i<ch_extra.forward.tail.size(); ++i)
                assert(ch_extra.forward.tail[i] < ch.forward.head[i]);
            for(unsigned i=0; i<ch_extra.backward.tail.size(); ++i)
                assert(ch_extra.backward.tail[i] < ch.backward.head[i]);
#endif

            if(log_message){
                timer += get_micro_time();
                log_message("Finished reordering nodes by rank. Needed "+std::to_string(timer)+"musec.");
            }
        }

        void sort_ch_arcs_and_build_first_out_arrays(
                ContractionHierarchy&ch,
                ContractionHierarchyExtraInfo&ch_extra,
                const std::function<void(std::string)>&log_message
        ){
            long long timer = 0; // initialize to avoid warning, not needed
            if(log_message){
                timer = -get_micro_time();
                log_message("Start sorting arcs.");
            }

            unsigned node_count = ch.rank.size();

            {
                auto r = compute_inverse_sort_permutation_first_by_tail_then_by_head_and_apply_sort_to_tail(node_count, ch_extra.forward.tail, ch.forward.head);

                ch.forward.head = apply_inverse_permutation(r, ch.forward.head);
                ch.forward.weight = apply_inverse_permutation(r, ch.forward.weight);
                ch_extra.forward.mid_node = apply_inverse_permutation(r, ch_extra.forward.mid_node);

                ch.forward.first_out = invert_vector(ch_extra.forward.tail, node_count);
            }

            {
                auto r = compute_inverse_sort_permutation_first_by_tail_then_by_head_and_apply_sort_to_tail(node_count, ch_extra.backward.tail, ch.backward.head);
                ch.backward.head = apply_inverse_permutation(r, ch.backward.head);
                ch.backward.weight = apply_inverse_permutation(r, ch.backward.weight);
                ch_extra.backward.mid_node = apply_inverse_permutation(r,ch_extra.backward.mid_node);

                ch.backward.first_out = invert_vector(ch_extra.backward.tail, node_count);
            }

            if(log_message){
                timer += get_micro_time();
                log_message("Finished sorting arcs. Needed "+std::to_string(timer)+"musec.");
            }
        }

        void optimize_order_for_cache(
                ContractionHierarchy&ch,
                const ContractionHierarchyExtraInfo&ch_extra,
                const std::function<void(std::string)>&log_message
        ){
            long long timer = 0; // initialize to avoid warning, not needed
            if(log_message){
                timer = -get_micro_time();
                log_message("Start optimizing order for cache.");
            }

            unsigned node_count = ch.rank.size();
            unsigned forward_arc_count = ch.forward.head.size();
            unsigned backward_arc_count = ch.backward.head.size();

            std::vector<bool>is_in_bottom_level(node_count, true);
            for(unsigned a=0; a<forward_arc_count; ++a)
                is_in_bottom_level[ch.forward.head[a]] = false;
            for(unsigned a=0; a<backward_arc_count; ++a)
                is_in_bottom_level[ch.backward.head[a]] = false;
            std::vector<unsigned>new_order(node_count);

            unsigned new_order_end = node_count;

            std::vector<bool>is_in_new_order(node_count, false);

            MinIDQueue q(node_count);
            for(unsigned r=0; r<node_count; ++r){
                if(is_in_bottom_level[r]){
                    unsigned search_space_end = new_order_end;

                    q.push({r, ch.rank[r]});
                    assert(!is_in_new_order[r]);
                    is_in_new_order[r] = true;

                    while(!q.empty()){
                        unsigned x = q.pop().id;

                        --new_order_end;
                        new_order[new_order_end] = x;

                        auto on_node = [&](unsigned y){
                            assert(!is_in_bottom_level[y]);
                            if(!is_in_new_order[y]){
                                is_in_new_order[y] = true;
                                q.push({y, ch.rank[y]});
                            }
                        };

                        for(unsigned xy = ch.forward.first_out[x]; xy < ch.forward.first_out[x+1]; ++xy)
                            on_node(ch.forward.head[xy]);

                        for(unsigned xy = ch.backward.first_out[x]; xy < ch.backward.first_out[x+1]; ++xy)
                            on_node(ch.backward.head[xy]);
                    }
                    std::reverse(new_order.begin()+new_order_end, new_order.begin()+search_space_end);
                }
            }

            assert(new_order_end == 0);
            assert(is_permutation(new_order));

            ch.rank = invert_permutation(new_order);
            ch.order = std::move(new_order);

            if(log_message){
                timer += get_micro_time();
                log_message("Finished optimizing order for cache. Needed "+std::to_string(timer)+"musec.");
            }
        }

        void build_unpacking_information(
                unsigned node_count,
                const std::vector<unsigned>&tail,
                const std::vector<unsigned>&head,
                const std::vector<unsigned>&input_arc_id,
                ContractionHierarchy&ch,
                const ContractionHierarchyExtraInfo&ch_extra,
                const std::function<void(std::string)>&log_message
        ){
            assert(is_sorted_using_less(tail));

            long long timer = 0;  // initialize to avoid warning, not needed
            if(log_message){
                log_message("Start building path unpacking information.");
                timer = -get_micro_time();
            }

            ch.forward.shortcut_first_arc = std::vector<unsigned>(ch.forward.head.size());
            ch.forward.shortcut_second_arc = std::vector<unsigned>(ch.forward.head.size());
            ch.forward.is_shortcut_an_original_arc = BitVector(ch.forward.head.size(), BitVector::uninitialized);
            ch.backward.shortcut_first_arc = std::vector<unsigned>(ch.backward.head.size());
            ch.backward.shortcut_second_arc = std::vector<unsigned>(ch.backward.head.size());
            ch.backward.is_shortcut_an_original_arc = BitVector(ch.backward.head.size(), BitVector::uninitialized);

            auto first_out = invert_vector(tail, node_count);

            for(unsigned x=0; x<node_count; ++x){
                for(unsigned xy=ch.forward.first_out[x]; xy<ch.forward.first_out[x+1]; ++xy){
                    unsigned y=ch.forward.head[xy];
                    unsigned z=ch_extra.forward.mid_node[xy];
                    if(z == invalid_id){
                        ch.forward.is_shortcut_an_original_arc.set(xy);

                        auto a = find_arc_given_sorted_head(first_out, head, ch.order[x], ch.order[y]);
                        ch.forward.shortcut_first_arc[xy] = input_arc_id[a];
                        ch.forward.shortcut_second_arc[xy] = head[a];
                    }else{
                        ch.forward.is_shortcut_an_original_arc.reset(xy);
                        ch.forward.shortcut_first_arc[xy] = find_arc_given_sorted_head(ch.backward.first_out, ch.backward.head, z, x);
                        ch.forward.shortcut_second_arc[xy] = find_arc_given_sorted_head(ch.forward.first_out, ch.forward.head, z, y);

                        assert(ch.forward.weight[xy] == ch.backward.weight[ch.forward.shortcut_first_arc[xy]] + ch.forward.weight[ch.forward.shortcut_second_arc[xy]]);

                    }
                }
            }

            for(unsigned x=0; x<node_count; ++x){
                for(unsigned xy=ch.backward.first_out[x]; xy<ch.backward.first_out[x+1]; ++xy){
                    unsigned y=ch.backward.head[xy];
                    unsigned z=ch_extra.backward.mid_node[xy];
                    if(z == invalid_id){
                        ch.backward.is_shortcut_an_original_arc.set(xy);
                        auto a = find_arc_given_sorted_head(first_out, head, ch.order[y], ch.order[x]);
                        ch.backward.shortcut_first_arc[xy] = input_arc_id[a];
                        ch.backward.shortcut_second_arc[xy] = head[a];
                    }else{
                        ch.backward.is_shortcut_an_original_arc.reset(xy);
                        ch.backward.shortcut_first_arc[xy] = find_arc_given_sorted_head(ch.backward.first_out, ch.backward.head, z, y);
                        ch.backward.shortcut_second_arc[xy] = find_arc_given_sorted_head(ch.forward.first_out, ch.forward.head, z, x);
                    }
                }
            }

#ifndef NDEBUG
            unsigned input_arc_count = max_element_of(input_arc_id, 0u)+1;

            for(unsigned a=0; a<ch.forward.head.size(); ++a){
                if(!ch.forward.is_shortcut_an_original_arc.is_set(a)){
                    assert(ch.forward.shortcut_first_arc[a] < ch.backward.head.size());
                    assert(ch.forward.shortcut_second_arc[a] < ch.forward.head.size());
                }else{
                    assert(ch.forward.shortcut_first_arc[a] < input_arc_count);
                    assert(ch.forward.shortcut_second_arc[a] < node_count);
                }
            }
            for(unsigned a=0; a<ch.backward.head.size(); ++a){
                if(!ch.backward.is_shortcut_an_original_arc.is_set(a)){
                    assert(ch.backward.shortcut_first_arc[a] < ch.backward.head.size());
                    assert(ch.backward.shortcut_second_arc[a] < ch.forward.head.size());
                }else{
                    assert(ch.backward.shortcut_first_arc[a] < input_arc_count);
                    assert(ch.backward.shortcut_second_arc[a] < node_count);
                }
            }
#endif


            if(log_message){
                timer += get_micro_time();
                log_message("Finished building path unpacking information. Needed "+std::to_string(timer)+"musec.");
                log_message("Contraction Hierarchy is fully constructed.");
            }
        }

        void log_input_graph_statistics(unsigned node_count, const std::vector<unsigned>&tail, const std::vector<unsigned>&head, const std::function<void(std::string)>&log_message){
            if(log_message){
                log_message("Input graph has "+std::to_string(node_count)+" nodes and "+std::to_string(tail.size())+" arcs.");
                std::vector<unsigned>deg(node_count, 0);
                for(unsigned i=0; i<tail.size(); ++i)
                    ++deg[tail[i]];
                unsigned max_out_degree = max_element_of(deg);
                std::fill(deg.begin(), deg.end(), 0);
                for(unsigned i=0; i<tail.size(); ++i)
                    ++deg[head[i]];
                unsigned max_in_degree = max_element_of(deg);
                log_message("The input's maximum in-degree is "+std::to_string(max_in_degree)+" and its maximum out-degree is "+std::to_string(max_out_degree)+".");

            }
        }

        void log_contraction_hierarchy_statistics(const ContractionHierarchy&ch, const std::function<void(std::string)>&log_message){
            if(log_message){
                log_message("CH has "+std::to_string(ch.forward.head.size())+" forward arcs.");
                log_message("CH has "+std::to_string(ch.backward.head.size())+" backward arcs.");
            }
        }
    }


    ContractionHierarchy ContractionHierarchy::build(
            unsigned node_count, std::vector<unsigned>tail, std::vector<unsigned>head, std::vector<unsigned>weight,
            const std::function<void(std::string)>&log_message, unsigned max_pop_count
    ){
        assert(tail.size() == head.size());
        assert(tail.size() == weight.size());
        assert(max_element_of(tail) < node_count);
        assert(max_element_of(head) < node_count);


        ContractionHierarchy ch;
        ContractionHierarchyExtraInfo ch_extra;

        log_input_graph_statistics(node_count, tail, head, log_message);

        std::vector<unsigned>input_arc_id = identity_permutation(head.size());

        {
            sort_arcs_and_remove_multi_and_loop_arcs(node_count, tail, head, weight, input_arc_id, log_message);
        }

        {
            Graph graph(node_count, tail, head, weight);
            build_ch_and_order(graph, ch, ch_extra, max_pop_count, log_message);
        }

        {
            // This optimizes the order in a postprocessing step
            sort_ch_arcs_and_build_first_out_arrays(ch, ch_extra, log_message);
            optimize_order_for_cache(ch, ch_extra, log_message);
        }

        {
            make_internal_nodes_and_rank_coincide(ch, ch_extra, log_message);
            sort_ch_arcs_and_build_first_out_arrays(ch, ch_extra, log_message);
        }

        build_unpacking_information(node_count, tail, head, input_arc_id, ch, ch_extra, log_message);

        log_contraction_hierarchy_statistics(ch, log_message);

        return ch;
    }

    ContractionHierarchy ContractionHierarchy::build_given_rank(
            std::vector<unsigned>rank,
            std::vector<unsigned>tail, std::vector<unsigned>head, std::vector<unsigned>weight,
            const std::function<void(std::string)>&log_message, unsigned max_pop_count
    ){
        unsigned node_count = rank.size();

        assert(tail.size() == head.size());
        assert(tail.size() == weight.size());
        assert(max_element_of(tail) < node_count);
        assert(max_element_of(head) < node_count);



        ContractionHierarchy ch;
        ContractionHierarchyExtraInfo ch_extra;


        log_input_graph_statistics(node_count, tail, head, log_message);

        std::vector<unsigned>input_arc_id = identity_permutation(head.size());

        {
            sort_arcs_and_remove_multi_and_loop_arcs(node_count, tail, head, weight, input_arc_id, log_message);
        }


        {
            Graph graph(node_count, tail, head, weight);
            build_ch_given_rank(graph, ch, ch_extra, rank, max_pop_count, log_message);
        }

        {
            make_internal_nodes_and_rank_coincide(ch, ch_extra, log_message);
            sort_ch_arcs_and_build_first_out_arrays(ch, ch_extra, log_message);
        }

        build_unpacking_information(node_count, tail, head, input_arc_id, ch, ch_extra, log_message);

        log_contraction_hierarchy_statistics(ch, log_message);

        return ch; // NVRO
    }

    ContractionHierarchy ContractionHierarchy::build_given_order(
            std::vector<unsigned>order,
            std::vector<unsigned>tail, std::vector<unsigned>head, std::vector<unsigned>weight,
            const std::function<void(std::string)>&log_message, unsigned max_pop_count
    ){
        return build_given_rank(invert_permutation(order), tail, head, weight, log_message, max_pop_count);
    }

    void check_contraction_hierarchy_for_errors(const ContractionHierarchy&ch){
        unsigned node_count = ch.rank.size();

        if(ch.rank != invert_permutation(ch.order))
            throw std::runtime_error("CH is invalid because: ch.rank != invert_permutation(ch.order)");

        if(ch.forward.first_out.size() != node_count+1)
            throw std::runtime_error("CH is invalid because: ch.forward.first_out.size() != node_count+1");
        if(ch.backward.first_out.size() != node_count+1)
            throw std::runtime_error("CH is invalid because: ch.backward.first_out.size() != node_count+1");

        unsigned forward_arc_count = ch.forward.first_out.back();

        if(ch.forward.first_out.front() != 0)
            throw std::runtime_error("CH is invalid because: ch.forward.first_out.front() != 0");
        if(!is_sorted_using_less(ch.forward.first_out))
            throw std::runtime_error("CH is invalid because: !is_sorted_using_less(ch.forward.first_out)");
        if(ch.forward.head.size() != forward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.forward.head.size() != forward_arc_count");
        if(ch.forward.weight.size() != forward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.forward.weight.size() != forward_arc_count");
        if(ch.forward.shortcut_first_arc.size() != forward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.forward.shortcut_first_arc.size() != forward_arc_count");
        if(ch.forward.shortcut_second_arc.size() != forward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.forward.shortcut_second_arc.size() != forward_arc_count");
        if(ch.forward.is_shortcut_an_original_arc.size() != forward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.forward.is_shortcut_an_original_arc.size() != forward_arc_count");
        if(!ch.forward.head.empty() && max_element_of(ch.forward.head) >= node_count)
            throw std::runtime_error("CH is invalid because: !ch.forward.head.empty() && max_element_of(ch.forward.head) >= node_count");

        unsigned backward_arc_count = ch.backward.first_out.back();

        if(ch.backward.first_out.front() != 0)
            throw std::runtime_error("CH is invalid because: ch.backward.first_out.front() != 0");
        if(!is_sorted_using_less(ch.backward.first_out))
            throw std::runtime_error("CH is invalid because: !is_sorted_using_less(ch.backward.first_out)");
        if(ch.backward.head.size() != backward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.backward.head.size() != backward_arc_count");
        if(ch.backward.weight.size() != backward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.backward.weight.size() != backward_arc_count");
        if(ch.backward.shortcut_first_arc.size() != backward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.backward.shortcut_first_arc.size() != backward_arc_count");
        if(ch.backward.shortcut_second_arc.size() != backward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.backward.shortcut_second_arc.size() != backward_arc_count");
        if(ch.backward.is_shortcut_an_original_arc.size() != backward_arc_count)
            throw std::runtime_error("CH is invalid because: ch.backward.is_shortcut_an_original_arc.size() != backward_arc_count");
        if(!ch.backward.head.empty() && max_element_of(ch.backward.head) >= node_count)
            throw std::runtime_error("CH is invalid because: !ch.backward.head.empty() && max_element_of(ch.backward.head) >= node_count");

        for(unsigned x=0; x<node_count; ++x){
            for(unsigned xy=ch.forward.first_out[x]; xy<ch.forward.first_out[x+1]; ++xy){
                unsigned y=ch.forward.head[xy];
                if(y <= x)
                    throw std::runtime_error("CH is invalid because: forward graph contains downward arc "+std::to_string(x)+" -> "+std::to_string(y));
            }
            for(unsigned xy=ch.backward.first_out[x]; xy<ch.backward.first_out[x+1]; ++xy){
                unsigned y=ch.backward.head[xy];
                if(y <= x)
                    throw std::runtime_error("CH is invalid because: backward graph contains downward arc "+std::to_string(x)+" -> "+std::to_string(y));
            }
        }

        for(unsigned xy=0; xy<forward_arc_count; ++xy){
            if(!ch.forward.is_shortcut_an_original_arc.is_set(xy)){
                if(ch.forward.shortcut_first_arc[xy] >= backward_arc_count)
                    throw std::runtime_error("CH is invalid because: ch.forward.shortcut_first_arc["+std::to_string(xy)+"] >= backward_arc_count");
                if(ch.forward.shortcut_second_arc[xy] >= forward_arc_count)
                    throw std::runtime_error("CH is invalid because: ch.forward.shortcut_second_arc["+std::to_string(xy)+"] >= forward_arc_count");
                if(ch.forward.shortcut_second_arc[xy] >= xy)
                    throw std::runtime_error("CH is invalid because: ch.forward.shortcut_second_arc["+std::to_string(xy)+"] >= "+std::to_string(xy));
                if(ch.forward.weight[xy] != ch.backward.weight[ch.forward.shortcut_first_arc[xy]] + ch.forward.weight[ch.forward.shortcut_second_arc[xy]])
                    throw std::runtime_error("CH is invalid because: ch.forward.weight[xy] != ch.backward.weight[ch.forward.shortcut_first_arc[xy]] + ch.forward.weight[ch.forward.shortcut_second_arc[xy]]");
            } else {
                if(ch.forward.shortcut_first_arc[xy] == invalid_id)
                    throw std::runtime_error("CH is invalid because: ch.forward.shortcut_first_arc["+std::to_string(xy)+"] == invalid_id for an original arc");
                if(ch.forward.shortcut_second_arc[xy] >= node_count)
                    throw std::runtime_error("CH is invalid because: ch.forward.shortcut_second_arc["+std::to_string(xy)+"] >= node_count for an original arc");
            }
        }

        for(unsigned xy=0; xy<backward_arc_count; ++xy){
            if(!ch.backward.is_shortcut_an_original_arc.is_set(xy)){
                if(ch.backward.shortcut_first_arc[xy] >= backward_arc_count)
                    throw std::runtime_error("CH is invalid because: ch.backward.shortcut_first_arc["+std::to_string(xy)+"] >= backward_arc_count");
                if(ch.backward.shortcut_second_arc[xy] >= forward_arc_count)
                    throw std::runtime_error("CH is invalid because: ch.backward.shortcut_second_arc["+std::to_string(xy)+"] >= forward_arc_count");
                if(ch.backward.shortcut_first_arc[xy] >= xy)
                    throw std::runtime_error("CH is invalid because: ch.backward.shortcut_first_arc["+std::to_string(xy)+"] >= "+std::to_string(xy));
                if(ch.backward.weight[xy] != ch.backward.weight[ch.backward.shortcut_first_arc[xy]] + ch.forward.weight[ch.backward.shortcut_second_arc[xy]])
                    throw std::runtime_error("CH is invalid because: ch.backward.weight[xy] != ch.backward.weight[ch.backward.shortcut_first_arc[xy]] + ch.forward.weight[ch.backward.shortcut_second_arc[xy]]");
            } else {
                if(ch.backward.shortcut_first_arc[xy] == invalid_id)
                    throw std::runtime_error("CH is invalid because: ch.backward.shortcut_first_arc["+std::to_string(xy)+"] == invalid_id for an original arc");
                if(ch.backward.shortcut_second_arc[xy] >= node_count)
                    throw std::runtime_error("CH is invalid because: ch.backward.shortcut_second_arc["+std::to_string(xy)+"] >= node_count for an original arc");
            }
        }

    }


    namespace {
        const unsigned long long ch_magic_number = 0x436f6e7448696572ull;

        struct CHFileHeader{
            unsigned long long magic_number;
            unsigned node_count;
            unsigned forward_arc_count;
            unsigned backward_arc_count;
        };
    }
    void ContractionHierarchy::print_graph(){
        for(int i = 0;i < rank.size();i++){
            std::cout<<"node_id: "<<i+1<< " rank: "<<rank[i]+1<<std::endl;
        }
        for(int i = 0;i < rank.size();i++){
            std::cout<< "node_id: "<<i+1<<std::endl;
            for(unsigned arc = forward.first_out[rank[i]]; arc < forward.first_out[rank[i]+1]; ++arc){
                std::cout<<i+1<<" "<<order[forward.head[arc]]+1 <<std::endl;
            }
            for(unsigned arc = backward.first_out[rank[i]]; arc < backward.first_out[rank[i]+1]; ++arc){
                std::cout<<i+1<<" "<<order[backward.head[arc]]+1 <<std::endl;
            }
        }



    }

    ContractionHierarchy ContractionHierarchy::read(std::istream&in){
        return read(
                [&](char*p, unsigned long long l){
                    if(!in.read(p, l))
                        throw std::runtime_error("std::istream::read failed while reading a contraction hierarchy");
                }
        );
    }

    ContractionHierarchy ContractionHierarchy::read(std::istream&in, unsigned long long file_size){
        return read(
                [&](char*p, unsigned long long l){
                    if(!in.read(p, l))
                        throw std::runtime_error("std::istream::read failed while reading a contraction hierarchy");
                },
                file_size
        );
    }

    void ContractionHierarchy::write(std::ostream&out) const {
        write(
                [&](const char*p, unsigned long long l){
                    if(!out.write(p, l))
                        throw std::runtime_error("std::ostream::write failed while reading a contraction hierarchy");
                }
        );
    }

    void ContractionHierarchy::convert_to_graph(vector<unsigned>&first_out, vector<unsigned>& head, vector<unsigned>&weight ){
        int number_of_node = order.size();
        int start_node_id =  0 ;
        top_n_percentage_node_order =  std::vector<unsigned>(number_of_node-start_node_id);
        top_n_percentage_node_mapper=  std::vector<unsigned>(number_of_node);
        std::fill(top_n_percentage_node_mapper.begin(), top_n_percentage_node_mapper.end(), 0);

        for(int i = start_node_id; i<number_of_node; i++){
            top_n_percentage_node_order[i-start_node_id] = order[i];
            //start with id = 1 ;
            top_n_percentage_node_mapper[i] = i-start_node_id + 1;
        }
        top_n_percentage_graph = Graph (number_of_node-start_node_id);
        for(int i = start_node_id; i<number_of_node; i++){
            //insert forward
            for(int j = forward.first_out[i]; j < forward.first_out[i+1];++j ) {
                unsigned head = forward.head[j];
                if (top_n_percentage_node_mapper[head] != 0){
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,forward.weight[j],j,
                                                       1);
                    top_n_percentage_graph.add_in_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,forward.weight[j],j,1);
                }
            }
            //insert backward
            for(int j = backward.first_out[i]; j < backward.first_out[i+1];++j ) {
                unsigned head = backward.head[j];
                if (top_n_percentage_node_mapper[head] != 0) {
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,backward.weight[j],j,0);
                    top_n_percentage_graph.add_in_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,backward.weight[j],j,0);
                }
            }
        }
        vector<unsigned> arc_id;
        vector<unsigned> is_forward;
        top_n_percentage_graph.graph_to_vector(first_out,head,weight,arc_id,is_forward);
    }
    void ContractionHierarchy::construct_top_n_percentage_graph(double percentage){
        int number_of_node = order.size();
        int start_node_id =  number_of_node * (1-percentage) ;
        top_n_percentage_node_order =  std::vector<unsigned>(number_of_node-start_node_id);
        top_n_percentage_node_mapper=  std::vector<unsigned>(number_of_node);
        std::fill(top_n_percentage_node_mapper.begin(), top_n_percentage_node_mapper.end(), 0);

        for(int i = start_node_id; i<number_of_node; i++){
            top_n_percentage_node_order[i-start_node_id] = order[i];
            //start with id = 1 ;
            top_n_percentage_node_mapper[i] = i-start_node_id + 1;
        }
        top_n_percentage_graph = Graph (number_of_node-start_node_id);
        for(int i = start_node_id; i<number_of_node; i++){
            //insert forward
            for(int j = forward.first_out[i]; j < forward.first_out[i+1];++j ) {
                unsigned head = forward.head[j];
                if (top_n_percentage_node_mapper[head] != 0){
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,forward.weight[j],j,1);
                    top_n_percentage_graph.add_in_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,forward.weight[j],j,1);
                }
            }
            //insert backward
            for(int j = backward.first_out[i]; j < backward.first_out[i+1];++j ) {
                unsigned head = backward.head[j];
                if (top_n_percentage_node_mapper[head] != 0) {
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,backward.weight[j],j,0);
                    top_n_percentage_graph.add_in_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,backward.weight[j],j,0);
                }
            }
        }
        reordering_top_n_percentage_graph_based_on_DFS(start_node_id);
        for (int i = 0; i < top_n_percentage_rank_mapper.size(); i ++){
            top_n_percentage_rank_mapper[i] = top_n_percentage_rank_mapper[i] + start_node_id;
        }

    }

    void ContractionHierarchy::construct_top_n_percentage_graph_without_reorder(double percentage){
        int number_of_node = order.size();
        int start_node_id =  number_of_node * (1-percentage) ;
        top_n_percentage_node_order =  std::vector<unsigned>(number_of_node-start_node_id);
        top_n_percentage_node_mapper=  std::vector<unsigned>(number_of_node);
        std::fill(top_n_percentage_node_mapper.begin(), top_n_percentage_node_mapper.end(), 0);

        for(int i = start_node_id; i<number_of_node; i++){
            top_n_percentage_node_order[i-start_node_id] = order[i];
            //start with id = 1 ;
            top_n_percentage_node_mapper[i] = i-start_node_id + 1;
        }
        top_n_percentage_graph = Graph (number_of_node-start_node_id);
        for(int i = start_node_id; i<number_of_node; i++){
            //insert forward
            for(int j = forward.first_out[i]; j < forward.first_out[i+1];++j ) {
                unsigned head = forward.head[j];
                if (top_n_percentage_node_mapper[head] != 0){
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,forward.weight[j],j,1);
                    top_n_percentage_graph.add_in_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,forward.weight[j],j,1);
                }
            }
            //insert backward
            for(int j = backward.first_out[i]; j < backward.first_out[i+1];++j ) {
                unsigned head = backward.head[j];
                if (top_n_percentage_node_mapper[head] != 0) {
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,backward.weight[j],j,0);
                    top_n_percentage_graph.add_in_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,backward.weight[j],j,0);
                }
            }
        }
//    reordering_top_n_percentage_graph_based_on_DFS(start_node_id);
//    for (int i = 0; i < top_n_percentage_rank_mapper.size(); i ++){
//        top_n_percentage_rank_mapper[i] = top_n_percentage_rank_mapper[i] + start_node_id;
//    }

    }


    std::vector<unsigned> ContractionHierarchy::get_dfs_ordering(){
        std::vector<unsigned> DFS_ordering = top_n_percentage_graph.get_DFS_ordering();
        if(is_permutation(DFS_ordering)){
            std::cout<<"DFS ordering is correct" << std::endl;
        }
        std::vector<unsigned> invert_DFS = invert_permutation(DFS_ordering);
        return invert_DFS;
    }


    void ContractionHierarchy::construct_top_n_percentage_plain_graph(double percentage){
        int number_of_node = order.size();
        int start_node_id =  number_of_node * (1-percentage) ;
        top_n_percentage_node_order =  std::vector<unsigned>(number_of_node-start_node_id);
        top_n_percentage_node_mapper=  std::vector<unsigned>(number_of_node);
        std::fill(top_n_percentage_node_mapper.begin(), top_n_percentage_node_mapper.end(), 0);

        for(int i = start_node_id; i<number_of_node; i++){
            top_n_percentage_node_order[i-start_node_id] = order[i];
            //start with id = 1 ;
            top_n_percentage_node_mapper[i] = i-start_node_id + 1;
        }
        top_n_percentage_graph = Graph (number_of_node-start_node_id);
        for(int i = start_node_id; i<number_of_node; i++){
            //insert forward
            for(int j = forward.first_out[i]; j < forward.first_out[i+1];++j ) {
                unsigned head = forward.head[j];
                if (top_n_percentage_node_mapper[head] != 0){
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,forward.weight[j],j,
                                                       1);
                    top_n_percentage_graph.add_in_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,forward.weight[j],j,1);
                }
            }
            //insert backward
            for(int j = backward.first_out[i]; j < backward.first_out[i+1];++j ) {
                unsigned head = backward.head[j];
                if (top_n_percentage_node_mapper[head] != 0) {
                    // only add top_n_percentage_node
                    top_n_percentage_graph.add_out_arc(top_n_percentage_node_mapper[head]-1,i-start_node_id,backward.weight[j],j,0);
                    top_n_percentage_graph.add_in_arc(i-start_node_id,top_n_percentage_node_mapper[head]-1,backward.weight[j],j,0);
                }
            }
        }
        top_n_percentage_rank_mapper =  std::vector<unsigned>(number_of_node-start_node_id);
        for (int i = 0; i < top_n_percentage_rank_mapper.size(); i ++){
            top_n_percentage_rank_mapper[i] = i + start_node_id;
        }

    }

    std::vector<unsigned> ContractionHierarchy::reordering(unsigned start_node_id,const vector<unsigned>& DFS_ordering ,const vector<unsigned>& invert_DFS ){
        std::vector<unsigned> new_order = std::vector<unsigned>(top_n_percentage_node_order.size());
        std::vector<unsigned> new_mapper = std::vector<unsigned>(top_n_percentage_node_mapper.size());
        top_n_percentage_rank_mapper = DFS_ordering;
        std::fill(new_mapper.begin(), new_mapper.end(), 0);
        for(unsigned i = 0; i < new_order.size(); i++){
            new_order[i] = top_n_percentage_node_order[DFS_ordering[i]];
        }


        for(unsigned i = start_node_id; i<new_mapper.size(); i++ ){
            new_mapper[i] = invert_DFS[i-start_node_id] + 1;
        }

        Graph new_graph = Graph(top_n_percentage_node_order.size());
        for(unsigned i = 0; i < new_order.size(); i++){
            for(unsigned j = 0; j < top_n_percentage_graph.out_deg(i); j++){
                auto out_arc = top_n_percentage_graph.out(i,j);
                new_graph.add_out_arc(invert_DFS[i],invert_DFS[out_arc.node],out_arc.weight,out_arc.arc_id,out_arc.is_forward);
            }
            for(unsigned j = 0; j < top_n_percentage_graph.in_deg(i); j++){
                auto in_arc = top_n_percentage_graph.in(i,j);
                new_graph.add_in_arc(invert_DFS[i],invert_DFS[in_arc.node],in_arc.weight,in_arc.arc_id,in_arc.is_forward);
            }
        }


        top_n_percentage_node_order = new_order;
        top_n_percentage_node_mapper = new_mapper;
        top_n_percentage_graph = new_graph;
        return DFS_ordering;

    }

    std::vector<unsigned> ContractionHierarchy::reordering_top_n_percentage_graph_based_on_DFS(unsigned start_node_id){
        std::vector<unsigned> DFS_ordering = top_n_percentage_graph.get_DFS_ordering();
        if(is_permutation(DFS_ordering)){
            std::cout<<"DFS ordering is correct" << std::endl;
        }
        std::vector<unsigned> invert_DFS = invert_permutation(DFS_ordering);
        std::vector<unsigned> new_order = std::vector<unsigned>(top_n_percentage_node_order.size());
        std::vector<unsigned> new_mapper = std::vector<unsigned>(top_n_percentage_node_mapper.size());
        top_n_percentage_rank_mapper = DFS_ordering;
        std::fill(new_mapper.begin(), new_mapper.end(), 0);
        for(unsigned i = 0; i < new_order.size(); i++){
            new_order[i] = top_n_percentage_node_order[DFS_ordering[i]];
        }


        for(unsigned i = start_node_id; i<new_mapper.size(); i++ ){
            new_mapper[i] = invert_DFS[i-start_node_id] + 1;
        }

        Graph new_graph = Graph(top_n_percentage_node_order.size());
        for(unsigned i = 0; i < new_order.size(); i++){
            for(unsigned j = 0; j < top_n_percentage_graph.out_deg(i); j++){
                auto out_arc = top_n_percentage_graph.out(i,j);
                new_graph.add_out_arc(invert_DFS[i],invert_DFS[out_arc.node],out_arc.weight,out_arc.arc_id,out_arc.is_forward);
            }
            for(unsigned j = 0; j < top_n_percentage_graph.in_deg(i); j++){
                auto in_arc = top_n_percentage_graph.in(i,j);
                new_graph.add_in_arc(invert_DFS[i],invert_DFS[in_arc.node],in_arc.weight,in_arc.arc_id,in_arc.is_forward);
            }
        }


        top_n_percentage_node_order = new_order;
        top_n_percentage_node_mapper = new_mapper;
        top_n_percentage_graph = new_graph;
        return DFS_ordering;

    }
    unsigned long long ContractionHierarchy::get_size(){
        unsigned long long  size = 0;
        size = size + forward.weight.size()* sizeof(forward.weight[0]);
        size = size + forward.head.size()* sizeof(forward.head[0]);
        size = size + forward.first_out.size()* sizeof(forward.first_out[0]);
        size = size + forward.shortcut_first_arc.size()* sizeof(forward.shortcut_first_arc[0]);
        size = size + forward.shortcut_second_arc.size()* sizeof(forward.shortcut_second_arc[0]);
        size = size + forward.is_shortcut_an_original_arc.size()/8;

        size = size + backward.weight.size()* sizeof(backward.weight[0]);
        size = size + backward.head.size()* sizeof(backward.head[0]);
        size = size + backward.first_out.size()* sizeof(backward.first_out[0]);
        size = size + backward.shortcut_first_arc.size()* sizeof(backward.shortcut_first_arc[0]);
        size = size + backward.shortcut_second_arc.size()* sizeof(backward.shortcut_second_arc[0]);
        size = size + backward.is_shortcut_an_original_arc.size()/8;
        size = size + order.size()* sizeof(order[0]);
        size = size + rank.size()*sizeof(rank[0]);
        return size;
    }

    ContractionHierarchy ContractionHierarchy::load_file(const std::string&file_name){
        ContractionHierarchy ch;
        open_file_for_loading(file_name, [&](std::istream&in, unsigned long long file_size){ch = read(in, file_size);});
        return ch;
    }


//double ContractionHierarchy::get_Euclidean_distance()
    void ContractionHierarchy::save_file(const std::string&file_name) const {
        open_file_for_saving(file_name, [&](std::ostream&out){write(out);});
    }

    namespace{
        void check_header(CHFileHeader header){
            if(header.magic_number != ch_magic_number)
                throw std::runtime_error("CH file magic number broken. Is this really a CH file?");
        }

        ContractionHierarchy finish_read(std::function<void(char*, unsigned long long)>in, CHFileHeader header){
            ContractionHierarchy ch;
            ch.rank = read_vector<unsigned>(in, header.node_count);
            ch.order = invert_permutation(ch.rank);

            ch.forward.first_out = read_vector<unsigned>(in, header.node_count+1);
            ch.forward.head = read_vector<unsigned>(in, header.forward_arc_count);
            ch.forward.weight = read_vector<unsigned>(in, header.forward_arc_count);
            ch.forward.is_shortcut_an_original_arc = read_bit_vector(in, header.forward_arc_count);
            ch.forward.shortcut_first_arc = read_vector<unsigned>(in, header.forward_arc_count);
            ch.forward.shortcut_second_arc = read_vector<unsigned>(in, header.forward_arc_count);

            ch.backward.first_out = read_vector<unsigned>(in, header.node_count+1);
            ch.backward.head = read_vector<unsigned>(in, header.backward_arc_count);
            ch.backward.weight = read_vector<unsigned>(in, header.backward_arc_count);
            ch.backward.is_shortcut_an_original_arc = read_bit_vector(in, header.backward_arc_count);
            ch.backward.shortcut_first_arc = read_vector<unsigned>(in, header.backward_arc_count);
            ch.backward.shortcut_second_arc = read_vector<unsigned>(in, header.backward_arc_count);

            return ch; // NVRO
        }
    }

    ContractionHierarchy ContractionHierarchy::read(std::function<void(char*, unsigned long long)>in, unsigned long long file_size){
        CHFileHeader header = read_value<CHFileHeader>(in);
        check_header(header);
        unsigned long long expected_file_size = (
                sizeof(CHFileHeader)
                + sizeof(unsigned)*(
                        header.node_count
                        + (
                                header.node_count+1 +
                                4*header.forward_arc_count

                        )
                        + (
                                header.node_count+1 +
                                4*header.backward_arc_count

                        )
                )
                + ((header.backward_arc_count+511)/512) * 64
                + ((header.forward_arc_count+511)/512) * 64
        );
        if(expected_file_size != file_size)
            throw std::runtime_error("CH file has a different size than specified in the header. This file is corrupt.");
        return finish_read(in, header);
    }

    ContractionHierarchy ContractionHierarchy::read(std::function<void(char*, unsigned long long)>in){
        CHFileHeader header = read_value<CHFileHeader>(in);
        check_header(header);
        return finish_read(in, header);
    }

    void ContractionHierarchy::write(std::function<void(const char*, unsigned long long)>out) const {
        CHFileHeader header;
        header.magic_number = ch_magic_number;
        header.node_count = forward.first_out.size()-1;
        header.forward_arc_count = forward.head.size();
        header.backward_arc_count = backward.head.size();

        write_value(out, header);
        write_vector(out, rank);

        write_vector(out, forward.first_out);
        write_vector(out, forward.head);
        write_vector(out, forward.weight);
        write_bit_vector(out, forward.is_shortcut_an_original_arc);
        write_vector(out, forward.shortcut_first_arc);
        write_vector(out, forward.shortcut_second_arc);

        write_vector(out, backward.first_out);
        write_vector(out, backward.head);
        write_vector(out, backward.weight);
        write_bit_vector(out, backward.is_shortcut_an_original_arc);
        write_vector(out, backward.shortcut_first_arc);
        write_vector(out, backward.shortcut_second_arc);
    }

    void ContractionHierarchy::load_cpd(FILE*f){
        top_n_percentage_cpd.load(f);
    }
    void ContractionHierarchy::load_cpd2(FILE*f){
        new_cpd.load(f);
    }

//    void ContractionHierarchy::load_cpd2(FILE*f){
//        correct_cpd.load(f);
//    }
    ContractionHierarchyQuery::ContractionHierarchyQuery(const ContractionHierarchy&ch):
            ch(&ch),
            was_forward_pushed(ch.node_count()), was_backward_pushed(ch.node_count()),
            forward_queue(ch.node_count()), backward_queue(ch.node_count()),
            forward_tentative_distance(ch.node_count()), backward_tentative_distance(ch.node_count()),
            forward_predecessor_type(ch.node_count()), backward_predecessor_type(ch.node_count()),
            forward_predecessor_node(ch.node_count()), backward_predecessor_node(ch.node_count()),
            forward_predecessor_arc(ch.node_count()), backward_predecessor_arc(ch.node_count()),
            shortest_path_meeting_node(invalid_id), path_meeting_nodes(),
            was_node_cached(ch.top_n_percentage_node_order.size()),distance_cache(ch.top_n_percentage_node_order.size()),

//    forward_caching(200,vector<unsigned>(ch.top_n_percentage_node_order.size())),
//    backward_caching(200,vector<unsigned>(ch.top_n_percentage_node_order.size())),
//    forward_time_flag(200,TimestampFlags(ch.top_n_percentage_node_order.size())),
//    backward_time_flag(200,TimestampFlags(ch.top_n_percentage_node_order.size())),
            cpd_distance_vector(1000),cpd_path_vector(1000),
            forward_cpd_nodes(200),backward_cpd_nodes(200),
            state(ContractionHierarchyQuery::InternalState::initialized)
    {}

    ContractionHierarchyQuery&ContractionHierarchyQuery::reset(){
        assert(ch && "query object must have an attached CH");

        was_forward_pushed.reset_all();
        forward_queue.clear();
        was_backward_pushed.reset_all();
        backward_queue.clear();
        was_node_cached.reset_all();

        shortest_path_meeting_node = invalid_id;
        path_meeting_nodes.clear();

        state = ContractionHierarchyQuery::InternalState::initialized;
        return *this;
    }

    ContractionHierarchyQuery&ContractionHierarchyQuery::reset(const ContractionHierarchy&new_ch){
        if(forward_tentative_distance.size() == new_ch.node_count()){
            reset();
            ch = &new_ch;
        } else {
            *this = ContractionHierarchyQuery(new_ch);
        }
        return *this;
    }

    ContractionHierarchyQuery&ContractionHierarchyQuery::add_source(unsigned external_s, unsigned dist_to_s){
        assert(ch && "query object must have an attached CH");
        assert(external_s < ch->node_count() && "node out of bounds");
        assert(state == ContractionHierarchyQuery::InternalState::initialized || state == ContractionHierarchyQuery::InternalState::target_pinned);

        unsigned s = ch->rank[external_s];
        source = s;
        if(!forward_queue.contains_id(s)){
            forward_queue.push({s, dist_to_s});
            forward_tentative_distance[s] = dist_to_s;
            forward_predecessor_node[s] = invalid_id;
        }else{
            if(dist_to_s < forward_tentative_distance[s]){
                forward_tentative_distance[s] = dist_to_s;
                forward_queue.decrease_key({s, dist_to_s});
            }
        }
        was_forward_pushed.set(s);
        return *this;
    }

    ContractionHierarchyQuery&ContractionHierarchyQuery::add_target(unsigned external_t, unsigned dist_to_t){
        assert(ch && "query object must have an attached CH");
        assert(external_t < ch->node_count() && "node out of bounds");
        assert(state == ContractionHierarchyQuery::InternalState::initialized || state == ContractionHierarchyQuery::InternalState::source_pinned);

        unsigned t = ch->rank[external_t];
        target =t;
        if(!backward_queue.contains_id(t)){
            backward_queue.push({t, dist_to_t});
            backward_tentative_distance[t] = dist_to_t;
            backward_predecessor_node[t] = invalid_id;
        }else{
            if(dist_to_t < backward_tentative_distance[t]){
                backward_tentative_distance[t] = dist_to_t;
                backward_queue.decrease_key({t, dist_to_t});
            }
        }

        was_backward_pushed.set(t);
        return *this;
    }




    namespace{

        template<class SetPred>
        void forward_expand_upward_ch_arcs_of_node(
                unsigned node,
                unsigned distance_to_node,
                const std::vector<unsigned>&forward_first_out,
                const std::vector<unsigned>&forward_head,
                const std::vector<unsigned>&forward_weight,
                TimestampFlags&was_forward_pushed,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&forward_tentative_distance,
                const SetPred&set_predecessor, unsigned& number_of_nodes_generated
        ){
            for(unsigned arc = forward_first_out[node]; arc < forward_first_out[node+1]; ++arc){
                unsigned h = forward_head[arc], d = distance_to_node + forward_weight[arc];
                if(was_forward_pushed.is_set(h)){
                    if(d < forward_tentative_distance[h]){
                        if(forward_queue.contains_id(h)){
                            forward_queue.decrease_key({h, d});
                        }else{
                            forward_queue.push({h, d});
                            number_of_nodes_generated++;
                        }
//                    forward_queue.decrease_key({h, d});
                        forward_tentative_distance[h] = d;
                        set_predecessor(h, node, arc);
                    }
                } else if(d < inf_weight){
                    forward_queue.push({h, d});
                    forward_tentative_distance[h] = d;
                    was_forward_pushed.set(h);
                    set_predecessor(h, node, arc);
                    number_of_nodes_generated++;
                }
            }
        }


        template<class SetPred,class GetLandmark>
        void forward_expand_upward_ch_arcs_of_node_with_landmark(
                unsigned node,
                unsigned distance_to_node,
                const std::vector<unsigned>&forward_first_out,
                const std::vector<unsigned>&forward_head,
                const std::vector<unsigned>&forward_weight,
                TimestampFlags&was_forward_pushed,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&forward_tentative_distance,
                const SetPred&set_predecessor,const GetLandmark&get_landmark_distance, const unsigned& shortest_distance, unsigned& number_of_nodes_generated
        ){
            for(unsigned arc = forward_first_out[node]; arc < forward_first_out[node+1]; ++arc){
                unsigned h = forward_head[arc], d = distance_to_node + forward_weight[arc];

                if(was_forward_pushed.is_set(h)){
                    unsigned forward_distance = forward_tentative_distance[h];
                    if(d < forward_distance ){
                        unsigned f = d + get_landmark_distance(h);
                        if(f < 1.5 * shortest_distance) {
                            if(forward_queue.contains_id(h)){
                                forward_queue.decrease_key({h, f,d});
                                forward_tentative_distance[h] = d;
                                set_predecessor(h, node, arc);
                                number_of_nodes_generated++;
                            }else{
                                forward_queue.push({h, f,d});
                                forward_tentative_distance[h] = d;
                                set_predecessor(h, node, arc);
                                number_of_nodes_generated++;
                            }
                        }
                    }
                } else {
                    unsigned f = d + get_landmark_distance(h);
                    if(f < 1.5 * shortest_distance) {
                        forward_queue.push({h, f,d});
                        forward_tentative_distance[h] = d;
                        was_forward_pushed.set(h);
                        set_predecessor(h, node, arc);
                        number_of_nodes_generated++;
                    }
                }

            }
        }



        template<class SetPred>
        void forward_expand_upward_ch_arcs_of_node(
                unsigned node,
                unsigned distance_to_node,
                const std::vector<unsigned>&forward_first_out,
                const std::vector<unsigned>&forward_head,
                const std::vector<unsigned>&forward_weight,
                TimestampFlags&was_forward_pushed,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&forward_tentative_distance,
                const SetPred&set_predecessor
        ){
            for(unsigned arc = forward_first_out[node]; arc < forward_first_out[node+1]; ++arc){
                unsigned h = forward_head[arc], d = distance_to_node + forward_weight[arc];
                if(was_forward_pushed.is_set(h)){
                    if(d < forward_tentative_distance[h]){
                        if(forward_queue.contains_id(h)){
                            forward_queue.decrease_key({h, d});
                        }else{
                            forward_queue.push({h, d});
                        }
//                    forward_queue.decrease_key({h, d});
                        forward_tentative_distance[h] = d;
                        set_predecessor(h, node, arc);
                    }
                } else if(d < inf_weight){
                    forward_queue.push({h, d});
                    forward_tentative_distance[h] = d;
                    was_forward_pushed.set(h);
                    set_predecessor(h, node, arc);
                }
            }
        }

        template<class SetPred>
        void
        forward_expand_upward_ch_arcs_of_node_count(
                unsigned node,
                unsigned distance_to_node,
                const std::vector<unsigned>&forward_first_out,
                const std::vector<unsigned>&forward_head,
                const std::vector<unsigned>&forward_weight,
                TimestampFlags&was_forward_pushed,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&forward_tentative_distance,
                const SetPred&set_predecessor,unsigned& number_of_node_generated
        ){
            for(unsigned arc = forward_first_out[node]; arc < forward_first_out[node+1]; ++arc){
                unsigned h = forward_head[arc], d = distance_to_node + forward_weight[arc];
                if(was_forward_pushed.is_set(h)){
                    if(d < forward_tentative_distance[h]){
                        forward_queue.decrease_key({h, d});
                        forward_tentative_distance[h] = d;
                        set_predecessor(h, node, arc);
                        number_of_node_generated++;
                    }
                } else if(d < inf_weight){
                    forward_queue.push({h, d});
                    forward_tentative_distance[h] = d;
                    was_forward_pushed.set(h);
                    set_predecessor(h, node, arc);
                    number_of_node_generated++;
                }
            }
        }






        bool forward_can_stall_at_node(
                unsigned node,
                const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
                const TimestampFlags&was_forward_pushed,
                std::vector<unsigned>&forward_tentative_distance, const std::vector<unsigned>&backward_tentative_distance
        ){
            for(unsigned arc = backward_first_out[node]; arc < backward_first_out[node+1]; ++arc){
                unsigned x = backward_head[arc];
                if(was_forward_pushed.is_set(x)){
                    if(forward_tentative_distance[x] + backward_weight[arc] <= forward_tentative_distance[node])
                        return true;
                }
            }
            return false;
        }

        bool forward_can_stall_at_node_all_path(
                unsigned node,
                const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
                const TimestampFlags&was_forward_pushed,
                std::vector<unsigned>&forward_tentative_distance, const std::vector<unsigned>&backward_tentative_distance
        ){
            for(unsigned arc = backward_first_out[node]; arc < backward_first_out[node+1]; ++arc){
                unsigned x = backward_head[arc];
                if(was_forward_pushed.is_set(x)){
                    if(forward_tentative_distance[x] + backward_weight[arc] < forward_tentative_distance[node])
                        return true;
                }
            }
            return false;
        }

        void forward_settle_node(
                unsigned&shortest_path_length,
                unsigned&shortest_path_meeting_node,
                const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
                const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
                TimestampFlags&was_forward_pushed, const TimestampFlags&was_backward_pushed,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&forward_tentative_distance, const std::vector<unsigned>&backward_tentative_distance,
                std::vector<unsigned>&forward_predecessor_node, std::vector<unsigned>&forward_predecessor_arc
        ){


            auto p = forward_queue.pop();
            auto popped_node = p.id;
            auto distance_to_popped_node = p.key;

            if(was_backward_pushed.is_set(popped_node)){
                if(shortest_path_length > distance_to_popped_node + backward_tentative_distance[popped_node]){
                    shortest_path_length = distance_to_popped_node + backward_tentative_distance[popped_node];
                    shortest_path_meeting_node = popped_node;
                }
            }

            if(
                    !forward_can_stall_at_node(
                            popped_node,
                            backward_first_out, backward_head, backward_weight,
                            was_forward_pushed,
                            forward_tentative_distance, backward_tentative_distance
                    )
                    )
                forward_expand_upward_ch_arcs_of_node(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        was_forward_pushed, forward_queue,
                        forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc){
                            forward_predecessor_node[x] = pred_node;
                            forward_predecessor_arc[x] = pred_arc;
                        }
                );
        }


        void full_forward_search(
                const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
                TimestampFlags&was_forward_pushed,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&forward_tentative_distance,
                std::vector<unsigned>&forward_predecessor_node, std::vector<unsigned>&forward_predecessor_arc
        ){
            while(!forward_queue.empty()){
                auto p = forward_queue.pop();
                auto popped_node = p.id;
                auto distance_to_popped_node = p.key;

                forward_expand_upward_ch_arcs_of_node(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        was_forward_pushed, forward_queue,
                        forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc){
                            forward_predecessor_node[x] = pred_node;
                            forward_predecessor_arc[x] = pred_arc;
                        }
                );
            }
        }




    }

    ContractionHierarchyQuery& ContractionHierarchyQuery::run(){
        assert(ch && "query object must have an attached CH");
        assert(!forward_queue.empty() && "must add at least one source before calling run");
        assert(!backward_queue.empty() && "must add at least one target before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);

        shortest_path_length = inf_weight;
        shortest_path_meeting_node = invalid_id;
        number_of_nodes_generated = 0;
        number_of_nodes_expanded =0;
        bool forward_next = true;


        for(;;){
            bool forward_finished = false;
            if(forward_queue.empty())
                forward_finished = true;
            else if(forward_queue.peek().key >= shortest_path_length)
                forward_finished = true;

            bool backward_finished = false;
            if(backward_queue.empty())
                backward_finished = true;
            else if(backward_queue.peek().key >= shortest_path_length)
                backward_finished = true;

            if(forward_finished && backward_finished)
                break;

            if(forward_finished)
                forward_next = false;
            if(backward_finished)
                forward_next = true;

            if(forward_next){
                forward_settle_node(
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        was_forward_pushed, was_backward_pushed,
                        forward_queue,
                        forward_tentative_distance, backward_tentative_distance,
                        forward_predecessor_node, forward_predecessor_arc, target
                );
                forward_next = false;
            } else {
                forward_settle_node(
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        was_backward_pushed, was_forward_pushed,
                        backward_queue,
                        backward_tentative_distance, forward_tentative_distance,
                        backward_predecessor_node, backward_predecessor_arc,source
                );
                forward_next = true;
            }
        }
        // cout << "shortest_path_meeting_node  " << shortest_path_meeting_node << endl;
        state = ContractionHierarchyQuery::InternalState::run;
        return *this;
    }

    ContractionHierarchyQuery& ContractionHierarchyQuery::run_with_landmark(const unsigned shortest_path_distance){
        assert(ch && "query object must have an attached CH");
        assert(!forward_queue.empty() && "must add at least one source before calling run");
        assert(!backward_queue.empty() && "must add at least one target before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);

        // shortest_path_length =  get_min_landmard_upperband(source,target);
        //cout << shortest_path_distance << endl;
        //cout << get_min_landmard_upperband(source,target) << endl;
        //cout << get_max_landmard_distance(source,target) << endl;
        shortest_path_length = shortest_path_distance;
        shortest_path_meeting_node = invalid_id;
        number_of_nodes_generated = 0;
        number_of_nodes_expanded =0;
        bool forward_next = true;

        for(;;){
            bool forward_finished = false;
            if(forward_queue.empty())
                forward_finished = true;
            else if(forward_queue.peek().key >= 1.5 * shortest_path_length)
                forward_finished = true;

            bool backward_finished = false;
            if(backward_queue.empty())
                backward_finished = true;
            else if(backward_queue.peek().key >= 1.5 * shortest_path_length)
                backward_finished = true;

            if(forward_finished && backward_finished)
                break;

            if(forward_finished)
                forward_next = false;
            if(backward_finished)
                forward_next = true;

            if(forward_next){
                forward_settle_node_with_landmark(
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        was_forward_pushed, was_backward_pushed,
                        forward_queue,
                        forward_tentative_distance, backward_tentative_distance,
                        forward_predecessor_node, forward_predecessor_arc,target
                );
                forward_next = false;
            } else {
                forward_settle_node_with_landmark(
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        was_backward_pushed, was_forward_pushed,
                        backward_queue,
                        backward_tentative_distance, forward_tentative_distance,
                        backward_predecessor_node, backward_predecessor_arc,source
                );
                forward_next = true;
            }
        }
        // cout << shortest_path_length << endl;
        state = ContractionHierarchyQuery::InternalState::run;
        return *this;
    }



    ContractionHierarchyQuery& ContractionHierarchyQuery::run_to_get_all_paths(){
        assert(ch && "query object must have an attached CH");
        assert(!forward_queue.empty() && "must add at least one source before calling run");
        assert(!backward_queue.empty() && "must add at least one target before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);

        shortest_path_length = inf_weight;
        path_meeting_nodes.clear();
        number_of_nodes_generated = 0;
        number_of_nodes_expanded =0;
        bool forward_next = true;

        for(;;){
            bool forward_finished = false;
            if(forward_queue.empty())
                forward_finished = true;
            else if(forward_queue.peek().key > shortest_path_length)
                forward_finished = true;

            bool backward_finished = false;
            if(backward_queue.empty())
                backward_finished = true;
            else if(backward_queue.peek().key > shortest_path_length)
                backward_finished = true;

            if(forward_finished && backward_finished)
                break;

            if(forward_finished)
                forward_next = false;
            if(backward_finished)
                forward_next = true;

            if(forward_next){
                forward_settle_node_all_path(
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        was_forward_pushed, was_backward_pushed,
                        forward_queue,
                        forward_tentative_distance, backward_tentative_distance,
                        forward_predecessor_node, forward_predecessor_arc,target
                );
                forward_next = false;
            } else {
                forward_settle_node_all_path(
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        was_backward_pushed, was_forward_pushed,
                        backward_queue,
                        backward_tentative_distance, forward_tentative_distance,
                        backward_predecessor_node, backward_predecessor_arc,source
                );
                forward_next = true;
            }
        }

        state = ContractionHierarchyQuery::InternalState::run;
        return *this;
    }
//

    bool ContractionHierarchyQuery::get_cpd_path(unsigned forward_start,unsigned backward_end,
                                                 vector<unsigned>&arc_path,vector<unsigned>&is_forward){
        unsigned cur_id = forward_start;
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, backward_end);
        if(next_move == 0xFF){
            return false;
        }
        unsigned next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
        for(;;){
            arc_path.push_back( ch->top_n_percentage_arc_id[next_id]);
            is_forward.push_back( ch->top_n_percentage_is_forward[next_id]);
            cur_id = ch->top_n_percentage_head[next_id];
            if(cur_id == backward_end){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, backward_end);
            next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
        }
        return true;
    }

    unsigned ContractionHierarchyQuery::run_cpd_search(unsigned source, unsigned target,vector<unsigned>&path){
        number_of_first_move_calls = 0;
        target = ch->top_n_percentage_node_mapper[target] ;
        unsigned cur_id = ch->top_n_percentage_node_mapper[source];
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
        number_of_first_move_calls++;
        if(next_move == 0xFF){
            return inf_weight;
        }
        path[0] = source;
        unsigned number_of_nodes = 1;
        unsigned cost = 0;
        for (;;){
            path[number_of_nodes] = ch->top_n_percentage_node_invert_mapper[cur_id];
            cost = cost +  ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id]+next_move];
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id]+next_move];
            if(cur_id == target){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
            number_of_first_move_calls++;
            number_of_nodes++;
        }
        path[number_of_nodes] =ch->top_n_percentage_node_invert_mapper[target];
        return cost;
    }

    unsigned ContractionHierarchyQuery::run_full_cpd_search(unsigned source, unsigned target,vector<unsigned>&path){
        number_of_first_move_calls = 0;
        target = ch->top_n_percentage_node_mapper[ch->rank[target]]-1 ;
        unsigned cur_id = ch->top_n_percentage_node_mapper[ch->rank[source]]-1;
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
        number_of_first_move_calls++;
        if(next_move == 0xFF){
            return inf_weight;
        }
        path[0] = source;
        unsigned number_of_nodes = 1;
        unsigned cost = 0;
        for (;;){
            path[number_of_nodes] = ch->top_n_percentage_node_order[cur_id];
            cost = cost +  ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id]+next_move];
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id]+next_move];
            if(cur_id == target){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
            number_of_first_move_calls++;
            number_of_nodes++;
        }
        path[number_of_nodes] = ch->top_n_percentage_node_order[target];
        return cost;
    }


    unsigned ContractionHierarchyQuery::run_full_cpd_search(unsigned source, unsigned target){
        number_of_first_move_calls = 0;
        target = ch->top_n_percentage_node_mapper[ch->rank[target]]-1 ;
        unsigned cur_id = ch->top_n_percentage_node_mapper[ch->rank[source]]-1;
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
        number_of_first_move_calls++;
        if(next_move == 0xFF){
            return inf_weight;
        }
        unsigned number_of_nodes = 1;
        unsigned cost = 0;
        for (;;){
            cost = cost +  ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id]+next_move];
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id]+next_move];
            if(cur_id == target){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
            number_of_first_move_calls++;
            number_of_nodes++;
        }
        return cost;
    }




    void ContractionHierarchyQuery::run_cpd_get_path(unsigned source, unsigned target,vector<unsigned>&path){
        number_of_first_move_calls = 0;
        target = ch->top_n_percentage_node_mapper[target] ;
        unsigned cur_id = ch->top_n_percentage_node_mapper[source];
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
        number_of_first_move_calls++;
        if(next_move == 0xFF){
            return;
        }
        path[0] = source;
        unsigned number_of_nodes = 1;
        for (;;){
            path[number_of_nodes] = ch->top_n_percentage_node_invert_mapper[cur_id];
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id]+next_move];
            if(cur_id == target){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
            number_of_first_move_calls++;
            number_of_nodes++;
        }
        path[number_of_nodes] =ch->top_n_percentage_node_invert_mapper[target];
    }

    unsigned ContractionHierarchyQuery::get_bi_cpd_path(unsigned forward_start,unsigned backward_end,
                                                        vector<unsigned>&node_path){
        bool is_forward = true;
        int signal = 0;
        vector<unsigned> forward_path =  vector<unsigned>();
        vector<unsigned> backward_path =  vector<unsigned>();
        unsigned cost = 0;

        unsigned start = ch->rank[ch->top_n_percentage_node_order[forward_start]];
        unsigned back =ch->rank[ch->top_n_percentage_node_order[backward_end]];


        for(;;){
            if(signal== -2){
                break;
            }
            if(is_forward){
                cost = cost + bi_cpd_path(forward_start, backward_end,forward_path,signal);
                is_forward= false;
            }else{
                cost = cost + bi_cpd_path(backward_end, forward_start,backward_path,signal);
                is_forward= true;
            }
        }


        if(forward_path.size() > 0 &&backward_path.size()>0){
            node_path.push_back(start);
            for(int i =0; i < forward_path.size(); i++){
                node_path.push_back(forward_path[i]);
            }
            if(forward_path[forward_path.size() -1] == backward_path[backward_path.size()-1]){
                for(int i = backward_path.size()-2; i>=0; i--){
                    node_path.push_back(backward_path[i]);
                }
            }
            node_path.push_back(back);
        }else if(forward_path.size() == 0) {
            for (int i = backward_path.size() - 1; i >= 0; i--) {
                node_path.push_back(backward_path[i]);
            }
            node_path.push_back(back);
        }else if(backward_path.size() == 0) {
            node_path.push_back(start);
            for(int i =0; i < forward_path.size(); i++){
                node_path.push_back(forward_path[i]);
            }
        }
        return cost;

    }

    unsigned ContractionHierarchyQuery::bi_cpd_path(unsigned& forward_start,unsigned& backward_end,
                                                    vector<unsigned>&node_path, int & signal){
        auto retrieve_next_move = [&](const int& source, const int& target) {
            if(source == target ){
                return -2;
            }
            const int& first_move = ch->new_cpd.get_first_move(source, target);
            return first_move -2;
        };

        unsigned cost = 0;
        int cur_id = forward_start;
        for(;;){
            int next_move = retrieve_next_move( cur_id,backward_end);
            if(next_move == -1){
                forward_start = cur_id;
                signal = -1;
                return cost;
            }else if(next_move == -2){
                forward_start = backward_end;
                //reached target
                signal = -2;
                return cost;
            }else{
                cost = cost + ch->top_n_percentage_graph.out(cur_id,next_move).weight;
//                unsigned previous_id = cur_id;
                cur_id = ch->top_n_percentage_graph.out(cur_id,next_move).node;
//                if(ch->rank[ch->top_n_percentage_node_order[previous_id]]>ch->rank[ch->top_n_percentage_node_order[cur_id]]){
//                    std::cout<< "error" << std::endl;
//                }
                node_path.push_back(ch->rank[ch->top_n_percentage_node_order[cur_id]]);
            }
        }
    }

    /*void ContractionHierarchyQuery::forward_settle_node(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            const unsigned target
    ){


        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];
        if(ch_was_backward_pushed.is_set(popped_node)){
            cout << popped_node << endl;
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }

        if(
                !forward_can_stall_at_node(
                        popped_node,
                        backward_first_out, backward_head, backward_weight,
                        ch_was_forward_pushed,
                        ch_forward_tentative_distance, ch_backward_tentative_distance
                )
                ) {
            number_of_nodes_expanded++;
            forward_expand_upward_ch_arcs_of_node_count(
                    popped_node, distance_to_popped_node,
                    forward_first_out, forward_head, forward_weight,
                    ch_was_forward_pushed, ch_forward_queue,
                    ch_forward_tentative_distance,
                    [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                        ch_forward_predecessor_node[x] = pred_node;
                        ch_forward_predecessor_arc[x] = pred_arc;
                    }, number_of_nodes_generated
            );
        }


    }*/

    void ContractionHierarchyQuery::forward_settle_node(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            const unsigned target
    ){

        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];
        if(ch_was_backward_pushed.is_set(popped_node)){
            //Search overlap
            path_meeting_nodes.insert(popped_node);
            /*if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }*/
        }
        number_of_nodes_expanded++;
        forward_expand_upward_ch_arcs_of_node_count(
                popped_node, distance_to_popped_node,
                forward_first_out, forward_head, forward_weight,
                ch_was_forward_pushed, ch_forward_queue,
                ch_forward_tentative_distance,
                [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                    ch_forward_predecessor_node[x] = pred_node;
                    ch_forward_predecessor_arc[x] = pred_arc;
                    }, number_of_nodes_generated
            );

    }
    void ContractionHierarchyQuery::forward_settle_node_with_landmark(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            const unsigned target
    ){


        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];
        if(ch_was_backward_pushed.is_set(popped_node)){
            //Search overlap
            // cout << popped_node << endl;
            path_meeting_nodes.insert(popped_node);
            /*if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }*/
        }



            number_of_nodes_expanded++;
            forward_expand_upward_ch_arcs_of_node_with_landmark(
                    popped_node, distance_to_popped_node,
                    forward_first_out, forward_head, forward_weight,
                    ch_was_forward_pushed, ch_forward_queue,
                    ch_forward_tentative_distance,
                    [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                        ch_forward_predecessor_node[x] = pred_node;
                        ch_forward_predecessor_arc[x] = pred_arc;
                    }, [&](unsigned start) {
                        return get_max_landmard_distance(start, target);
                    }, shortest_path_length, number_of_nodes_generated
            );

    }
    void ContractionHierarchyQuery::forward_settle_node_all_path(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            const unsigned target
    ){


        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];
        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length >= distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];

                path_meeting_nodes.insert(popped_node);
            }
        }

        if(
                !forward_can_stall_at_node_all_path(
                        popped_node,
                        backward_first_out, backward_head, backward_weight,
                        ch_was_forward_pushed,
                        ch_forward_tentative_distance, ch_backward_tentative_distance
                )
                )

//            forward_expand_upward_ch_arcs_of_node(
//                    popped_node, distance_to_popped_node,
//                    forward_first_out, forward_head, forward_weight,
//                    ch_was_forward_pushed, ch_forward_queue,
//                    ch_forward_tentative_distance,
//                    [&](unsigned x, unsigned pred_node, unsigned pred_arc){
//                        ch_forward_predecessor_node[x] = pred_node;
//                        ch_forward_predecessor_arc[x] = pred_arc;
//                    }
//            );
            number_of_nodes_expanded ++;
        forward_expand_upward_ch_arcs_of_node(
                popped_node, distance_to_popped_node,
                forward_first_out, forward_head, forward_weight,
                ch_was_forward_pushed, ch_forward_queue,
                ch_forward_tentative_distance,
                [&](unsigned x, unsigned pred_node, unsigned pred_arc){
                    ch_forward_predecessor_node[x] = pred_node;
                    ch_forward_predecessor_arc[x] = pred_arc;
                }
        );

    }


    void ContractionHierarchyQuery::forward_settle_node_with_CPD(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            std::vector<unsigned>&ch_forward_cpd_nodes, const std::vector<unsigned>&ch_backward_cpd_nodes,
            unsigned & ch_number_of_forward_cpd_nodes, const unsigned & ch_number_of_backward_cpd_nodes,
            unsigned& ch_best_cpd_forward,unsigned& ch_best_cpd_backward, unsigned& target

    ){

        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
//        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];
        auto distance_to_popped_node = p.g;
        if(distance_to_popped_node > ch_forward_tentative_distance[popped_node] || p.key >=shortest_path_length){
            return;
        }
        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }

        if (
                !forward_can_stall_at_node(
                        popped_node,
                        backward_first_out, backward_head, backward_weight,
                        ch_was_forward_pushed,
                        ch_forward_tentative_distance, ch_backward_tentative_distance
                )
                )
        {
            if (popped_node >=start_id ) {
                if(distance_to_popped_node+get_max_landmard_distance(popped_node,target)  < shortest_path_length){

//                    update_cpd_distance(popped_node, ch_backward_cpd_nodes, ch_number_of_backward_cpd_nodes,
//                                        distance_to_popped_node,ch_backward_tentative_distance,ch_best_cpd_forward,ch_best_cpd_backward);
                    //pushed into CPD node
//                        update_cpd_distance(popped_node, ch_backward_cpd_nodes, ch_number_of_backward_cpd_nodes,
//                                            distance_to_popped_node,ch_best_cpd_forward,ch_best_cpd_backward,
//                                            ch_forward_tentative_distance,ch_backward_tentative_distance,
//                                            ch_was_forward_pushed,ch_was_backward_pushed);
                    was_node_cached.reset_all();
                    for(int i = 0 ; i < ch_number_of_backward_cpd_nodes; i++){
                        unsigned backward_node = ch_backward_cpd_nodes[i];
                        if(distance_to_popped_node + get_max_landmard_distance(popped_node,backward_node) + ch_backward_tentative_distance[backward_node] > shortest_path_length){
                            continue;
                        }
                        // always calculate distance from backward to forward !
                        unsigned current_distance =ch_backward_tentative_distance[backward_node]+
                                                   get_cpd_distance_with_cache(ch->top_n_percentage_node_mapper[popped_node] -1,
                                                                               ch->top_n_percentage_node_mapper[backward_node] -1,
                                                                               ch_forward_tentative_distance,ch_was_forward_pushed,distance_to_popped_node
                                                           ,ch_forward_predecessor_node,ch_forward_predecessor_arc, backward_predecessor_type
                                                   );

                        number_of_path_extraction ++;
                        if(shortest_path_length > current_distance){
                            shortest_path_meeting_node = invalid_id;
                            ch_best_cpd_forward = popped_node;
                            ch_best_cpd_backward = backward_node;
                            shortest_path_length = current_distance;
                        }
                    }
                    ch_forward_cpd_nodes[ch_number_of_forward_cpd_nodes] = popped_node;
                    ch_number_of_forward_cpd_nodes++;
                }

//                timer1.stop();
//                cpd_time = cpd_time + timer1.elapsed_time_micro();
            } else {
                // only expand non cpd node
//                forward_expand_upward_ch_arcs_of_node(
//                        popped_node, distance_to_popped_node,
//                        forward_first_out, forward_head, forward_weight,
//                        ch_was_forward_pushed, ch_forward_queue,
//                        ch_forward_tentative_distance,
//                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
//                            ch_forward_predecessor_node[x] = pred_node;
//                            ch_forward_predecessor_arc[x] = pred_arc;
//                        },number_of_nodes_generated
//                );
                number_of_nodes_expanded ++;
                forward_expand_upward_ch_arcs_of_node_with_landmark(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        ch_was_forward_pushed, ch_forward_queue,
                        ch_forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                            ch_forward_predecessor_node[x] = pred_node;
                            ch_forward_predecessor_arc[x] = pred_arc;
                        }, [&](unsigned start) {
                            return get_max_landmard_distance(start,target);
                        },shortest_path_length,number_of_nodes_generated
                );
            }
        }
    }


    void ContractionHierarchyQuery::forward_settle_node_with_CPD(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            std::vector<bool>&ch_forward_predecessor_type,
            std::vector<unsigned>&ch_forward_cpd_nodes, const std::vector<unsigned>&ch_backward_cpd_nodes,
            unsigned & ch_number_of_forward_cpd_nodes, const unsigned & ch_number_of_backward_cpd_nodes,
            unsigned& target,const bool is_forward

    ){

        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
        auto distance_to_popped_node = p.g;
        if(distance_to_popped_node > ch_forward_tentative_distance[popped_node] || p.key >=shortest_path_length){
            return;
        }
        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }
        if ( !forward_can_stall_at_node( popped_node,
                                         backward_first_out, backward_head, backward_weight,
                                         ch_was_forward_pushed,
                                         ch_forward_tentative_distance, ch_backward_tentative_distance )
                )
        {
            if (popped_node >= start_id) {
                    for (int i = 0; i < ch_number_of_backward_cpd_nodes; i++) {
                        unsigned backward_node = ch_backward_cpd_nodes[i];
                        if (distance_to_popped_node + get_max_landmard_distance(popped_node, backward_node) +
                            ch_backward_tentative_distance[backward_node] >= shortest_path_length) {
                            continue;
                        }
                        number_of_path_extraction++;
                        unsigned current_distance = ch_backward_tentative_distance[backward_node] +
                                                    get_cpd_distance_with_cache(
                                                            ch->top_n_percentage_node_mapper[popped_node] - 1,
                                                            ch->top_n_percentage_node_mapper[backward_node] - 1,
                                                            ch_forward_tentative_distance, ch_was_forward_pushed,
                                                            distance_to_popped_node, ch_forward_predecessor_node,
                                                            ch_forward_predecessor_arc, ch_forward_predecessor_type
                                                    );
                        if (shortest_path_length > current_distance) {
                            shortest_path_meeting_node = backward_node;
                            shortest_path_length = current_distance;
                        }
                    }
                    ch_forward_cpd_nodes[ch_number_of_forward_cpd_nodes] = popped_node;
                    ch_number_of_forward_cpd_nodes++;
            } else {
                // only expand non cpd node
                number_of_nodes_expanded ++;
                forward_expand_upward_ch_arcs_of_node_with_landmark(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        ch_was_forward_pushed, ch_forward_queue,
                        ch_forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                            ch_forward_predecessor_node[x] = pred_node;
                            ch_forward_predecessor_arc[x] = pred_arc;
                            ch_forward_predecessor_type[x]  = is_forward;
                        }, [&](unsigned start) {
                            return get_max_landmard_distance(start,target);
                        },shortest_path_length,number_of_nodes_generated
                );
            }
        }
    }



    void ContractionHierarchyQuery::backward_settle_node_with_CPD(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            std::vector<bool>&ch_forward_predecessor_type,
            std::vector<unsigned>&ch_backward_predecessor_node, std::vector<unsigned>&ch_backward_predecessor_arc,
            std::vector<bool>&ch_backward_predecessor_type,
            std::vector<unsigned>&ch_forward_cpd_nodes, const std::vector<unsigned>&ch_backward_cpd_nodes,
            unsigned & ch_number_of_forward_cpd_nodes, const unsigned & ch_number_of_backward_cpd_nodes,
            unsigned& target,const bool is_forward

    ){

        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
        auto distance_to_popped_node = p.g;
        if(distance_to_popped_node > ch_forward_tentative_distance[popped_node] || p.key >=shortest_path_length){
            return;
        }
        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }
        if ( !forward_can_stall_at_node( popped_node,
                                         backward_first_out, backward_head, backward_weight,
                                         ch_was_forward_pushed,
                                         ch_forward_tentative_distance, ch_backward_tentative_distance )
                )
        {
            if (popped_node >= start_id) {
                for (int i = 0; i < ch_number_of_backward_cpd_nodes; i++) {
                    unsigned backward_node = ch_backward_cpd_nodes[i];
                    if (distance_to_popped_node + get_max_landmard_distance(popped_node, backward_node) +
                        ch_backward_tentative_distance[backward_node] >= shortest_path_length) {
                        continue;
                    }
                    number_of_path_extraction++;
                    unsigned current_distance = distance_to_popped_node+
                                                get_cpd_distance_with_cache(
                                                        ch->top_n_percentage_node_mapper[backward_node] - 1,
                                                        ch->top_n_percentage_node_mapper[popped_node] - 1,
                                                        ch_backward_tentative_distance, ch_was_backward_pushed,
                                                        ch_backward_tentative_distance[backward_node] , ch_backward_predecessor_node,
                                                        ch_backward_predecessor_arc, ch_backward_predecessor_type
                                                );
                    if (shortest_path_length > current_distance) {
                        shortest_path_meeting_node = popped_node;
                        shortest_path_length = current_distance;
                    }
                }
                ch_forward_cpd_nodes[ch_number_of_forward_cpd_nodes] = popped_node;
                ch_number_of_forward_cpd_nodes++;
            } else {
                // only expand non cpd node
                number_of_nodes_expanded ++;
                forward_expand_upward_ch_arcs_of_node_with_landmark(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        ch_was_forward_pushed, ch_forward_queue,
                        ch_forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                            ch_forward_predecessor_node[x] = pred_node;
                            ch_forward_predecessor_arc[x] = pred_arc;
                            ch_forward_predecessor_type[x]  = is_forward;
                        }, [&](unsigned start) {
                            return get_max_landmard_distance(start,target);
                        },shortest_path_length,number_of_nodes_generated
                );
            }
        }
    }





    void ContractionHierarchyQuery::forward_settle_node_with_CPD(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            std::vector<bool>&ch_forward_predecessor_type,
            std::vector<unsigned>&ch_forward_cpd_nodes, const std::vector<unsigned>&ch_backward_cpd_nodes,
            unsigned & ch_number_of_forward_cpd_nodes, const unsigned & ch_number_of_backward_cpd_nodes,
            unsigned& ch_best_cpd_forward,unsigned& ch_best_cpd_backward, unsigned& target,const bool is_forward

    ){

        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
        auto distance_to_popped_node = p.g;
        if(distance_to_popped_node > ch_forward_tentative_distance[popped_node] || p.key >=shortest_path_length){
            return;
        }
        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }

        if ( !forward_can_stall_at_node( popped_node,
                                         backward_first_out, backward_head, backward_weight,
                                         ch_was_forward_pushed,
                                         ch_forward_tentative_distance, ch_backward_tentative_distance )
                )
        {
            if (popped_node >= start_id ) {
                if(distance_to_popped_node+get_max_landmard_distance(popped_node,target)  < shortest_path_length) {
                    for (int i = 0; i < ch_number_of_backward_cpd_nodes; i++) {
                        unsigned backward_node = ch_backward_cpd_nodes[i];
                        if (distance_to_popped_node + get_max_landmard_distance(popped_node, backward_node) +
                            ch_backward_tentative_distance[backward_node] >= shortest_path_length) {
                            continue;
                        }
                        number_of_path_extraction++;
                        // always calculate distance from backward to forward !
                        unsigned current_distance = ch_backward_tentative_distance[backward_node] +
                                                    get_cpd_distance_with_cache(
                                                            ch->top_n_percentage_node_mapper[popped_node] - 1,
                                                            ch->top_n_percentage_node_mapper[backward_node] - 1,
                                                            ch_forward_tentative_distance, ch_was_forward_pushed,
                                                            distance_to_popped_node, ch_forward_predecessor_node,
                                                            ch_forward_predecessor_arc, ch_forward_predecessor_type
                                                    );
                        if (shortest_path_length > current_distance) {
                            shortest_path_meeting_node = backward_node;
//                            ch_best_cpd_forward = popped_node;
//                            ch_best_cpd_backward = backward_node;
                            shortest_path_length = current_distance;
                        }
                    }
                    ch_forward_cpd_nodes[ch_number_of_forward_cpd_nodes] = popped_node;
                    ch_number_of_forward_cpd_nodes++;
                }
            } else {
                // only expand non cpd node
                number_of_nodes_expanded ++;
                forward_expand_upward_ch_arcs_of_node_with_landmark(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        ch_was_forward_pushed, ch_forward_queue,
                        ch_forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                            ch_forward_predecessor_node[x] = pred_node;
                            ch_forward_predecessor_arc[x] = pred_arc;
                            ch_forward_predecessor_type[x]  = is_forward;
                        }, [&](unsigned start) {
                            return get_max_landmard_distance(start,target);
                        },shortest_path_length,number_of_nodes_generated
                );
            }
        }
    }

    void ContractionHierarchyQuery::forward_settle_node_with_CPD_bi_caching(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            std::vector<unsigned>&ch_forward_cpd_nodes, const std::vector<unsigned>&ch_backward_cpd_nodes,
            unsigned & ch_number_of_forward_cpd_nodes, const unsigned & ch_number_of_backward_cpd_nodes,
            unsigned& ch_best_cpd_forward,unsigned& ch_best_cpd_backward,unsigned target,
            std::vector<vector<unsigned>>& f_cache, std::vector<vector<unsigned>>& b_cache,
            std::vector<TimestampFlags>&f_flag,std::vector<TimestampFlags>&b_flag

    ){

        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];

        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }

        if (
                !forward_can_stall_at_node(
                        popped_node,
                        backward_first_out, backward_head, backward_weight,
                        ch_was_forward_pushed,
                        ch_forward_tentative_distance, ch_backward_tentative_distance
                )
                )
        {
            if (ch->top_n_percentage_node_mapper[popped_node] != 0 ) {


                if(distance_to_popped_node+get_max_landmard_distance(popped_node,target)  < shortest_path_length){

//                    update_cpd_distance(popped_node, ch_backward_cpd_nodes, ch_number_of_backward_cpd_nodes,
//                                        distance_to_popped_node,ch_backward_tentative_distance,ch_best_cpd_forward,ch_best_cpd_backward);
                    //pushed into CPD node
//                        update_cpd_distance(popped_node, ch_backward_cpd_nodes, ch_number_of_backward_cpd_nodes,
//                                            distance_to_popped_node,ch_best_cpd_forward,ch_best_cpd_backward,
//                                            ch_forward_tentative_distance,ch_backward_tentative_distance,
//                                            ch_was_forward_pushed,ch_was_backward_pushed);
                    was_node_cached.reset_all();
                    for(int i = 0 ; i < ch_number_of_backward_cpd_nodes; i++){
                        unsigned backward_node = ch_backward_cpd_nodes[i];
                        if(distance_to_popped_node + get_max_landmard_distance(popped_node,backward_node) + ch_backward_tentative_distance[backward_node] > shortest_path_length){
                            continue;
                        }
                        unsigned d1 = get_cpd_distance_with_bi_cache(ch->top_n_percentage_node_mapper[popped_node] -1,ch->top_n_percentage_node_mapper[backward_node] -1,b_cache[i],b_flag[i]);
//                        unsigned d2 = get_cpd_distance_with_cache(ch->top_n_percentage_node_mapper[backward_node] -1,ch->top_n_percentage_node_mapper[popped_node] -1);
//                        if(d1 != d2 ){
//                            std::cout<<"eer"<<std::endl;
//                        }
                        // always calculate distance from backward to forward !
                        unsigned current_distance = distance_to_popped_node + d1
                                                    + ch_backward_tentative_distance[backward_node];
                        number_of_path_extraction ++;
                        if(shortest_path_length > current_distance){
                            shortest_path_meeting_node = invalid_id;
                            ch_best_cpd_forward = popped_node;
                            ch_best_cpd_backward = backward_node;
                            shortest_path_length = current_distance;
                        }
                    }
                    ch_forward_cpd_nodes[ch_number_of_forward_cpd_nodes] = popped_node;
//                    f_cache[ch_number_of_forward_cpd_nodes] = vector<unsigned>(forward_first_out.size());
                    ch_number_of_forward_cpd_nodes++;
                }

//                timer1.stop();
//                cpd_time = cpd_time + timer1.elapsed_time_micro();
            } else {
                // only expand non cpd node
//                forward_expand_upward_ch_arcs_of_node(
//                        popped_node, distance_to_popped_node,
//                        forward_first_out, forward_head, forward_weight,
//                        ch_was_forward_pushed, ch_forward_queue,
//                        ch_forward_tentative_distance,
//                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
//                            ch_forward_predecessor_node[x] = pred_node;
//                            ch_forward_predecessor_arc[x] = pred_arc;
//                        },number_of_nodes_generated
//                );
                number_of_nodes_expanded ++;
                forward_expand_upward_ch_arcs_of_node_with_landmark(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        ch_was_forward_pushed, ch_forward_queue,
                        ch_forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                            ch_forward_predecessor_node[x] = pred_node;
                            ch_forward_predecessor_arc[x] = pred_arc;
                        }, [&](unsigned start) {
                            return get_max_landmard_distance(start,target);
                        },shortest_path_length,number_of_nodes_generated
                );
            }
        }
    }

    void ContractionHierarchyQuery::forward_settle_node_until_find_CPD(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            std::vector<unsigned>&ch_forward_cpd_nodes, const std::vector<unsigned>&ch_backward_cpd_nodes,
            unsigned & ch_number_of_forward_cpd_nodes, const unsigned & ch_number_of_backward_cpd_nodes,
            unsigned& ch_best_cpd_forward,unsigned& ch_best_cpd_backward, unsigned& target

    ){

//        forward_queue.push({h, d});
//        forward_tentative_distance[h] = d;
//        was_forward_pushed.set(h);
//        set_predecessor(h, node, arc);
        for(;;) {
            if(ch_forward_queue.empty()){
                break;
            }

            auto p = ch_forward_queue.pop();
            auto popped_node = p.id;
            auto distance_to_popped_node = p.key;


            if (ch_was_backward_pushed.is_set(popped_node)) {
                if (shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]) {
                    shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                    shortest_path_meeting_node = popped_node;
                }
            }

            if (
                    !forward_can_stall_at_node(
                            popped_node,
                            backward_first_out, backward_head, backward_weight,
                            ch_was_forward_pushed,
                            ch_forward_tentative_distance, ch_backward_tentative_distance
                    )
                    ) {
                if (ch->top_n_percentage_node_mapper[popped_node] != 0) {
                    timer1.start();
                    if (distance_to_popped_node + get_max_landmard_distance(ch->order[popped_node], target) <
                        shortest_path_length) {

                        update_cpd_distance(popped_node, ch_backward_cpd_nodes, ch_number_of_backward_cpd_nodes,
                                            distance_to_popped_node, ch_backward_tentative_distance,
                                            ch_best_cpd_forward, ch_best_cpd_backward);
//
                        ch_forward_cpd_nodes[ch_number_of_forward_cpd_nodes] = popped_node;
                        ch_number_of_forward_cpd_nodes++;
                        timer1.stop();
                        cpd_time = cpd_time + timer1.elapsed_time_micro();
                        break;
                    }
                    timer1.stop();
                    cpd_time = cpd_time + timer1.elapsed_time_micro();
                    break;

                } else {
                    // only expand non cpd node
                    forward_expand_upward_ch_arcs_of_node(
                            popped_node, distance_to_popped_node,
                            forward_first_out, forward_head, forward_weight,
                            ch_was_forward_pushed, ch_forward_queue,
                            ch_forward_tentative_distance,
                            [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                                ch_forward_predecessor_node[x] = pred_node;
                                ch_forward_predecessor_arc[x] = pred_arc;
                            }
                    );
                }
            }
        }
    }



    void ContractionHierarchyQuery::forward_settle_node_with_direction_cpd(
            const std::vector<unsigned>&forward_first_out, const std::vector<unsigned>&forward_head, const std::vector<unsigned>&forward_weight,
            const std::vector<unsigned>&backward_first_out, const std::vector<unsigned>&backward_head, const std::vector<unsigned>&backward_weight,
            TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
            MinIDQueue&ch_forward_queue,
            std::vector<unsigned>&ch_forward_tentative_distance, std::vector<unsigned>&ch_backward_tentative_distance,
            std::vector<unsigned>&ch_forward_predecessor_node, std::vector<unsigned>&ch_forward_predecessor_arc,
            const unsigned target
    ){


        auto p = ch_forward_queue.pop();
        auto popped_node = p.id;
//        auto distance_to_popped_node = p.key;
        auto distance_to_popped_node = ch_forward_tentative_distance[popped_node];
        if(ch_was_backward_pushed.is_set(popped_node)){
            if(shortest_path_length > distance_to_popped_node + ch_backward_tentative_distance[popped_node]){
                shortest_path_length = distance_to_popped_node + ch_backward_tentative_distance[popped_node];
                shortest_path_meeting_node = popped_node;
            }
        }

        if(
                !forward_can_stall_at_node(
                        popped_node,
                        backward_first_out, backward_head, backward_weight,
                        ch_was_forward_pushed,
                        ch_forward_tentative_distance, ch_backward_tentative_distance
                )
                )
        {
            if (popped_node >= start_id) {
                unsigned fm = ch->top_n_percentage_cpd.get_first_move(popped_node-start_id,ch->top_n_percentage_node_mapper[target]-1);
                number_of_first_move_calls ++;
                if(fm == 2){
                    number_of_nodes_expanded++;
                    forward_expand_upward_ch_arcs_of_node_with_landmark(
                            popped_node, distance_to_popped_node,
                            forward_first_out, forward_head, forward_weight,
                            ch_was_forward_pushed, ch_forward_queue,
                            ch_forward_tentative_distance,
                            [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                                ch_forward_predecessor_node[x] = pred_node;
                                ch_forward_predecessor_arc[x] = pred_arc;
                            }, [&](unsigned start) {
                                return get_max_landmard_distance(start, target);
                            }, shortest_path_length, number_of_nodes_generated
                    );
                }else{
                    number_of_path_extraction++;
                }

            }else{
                number_of_nodes_expanded++;
                forward_expand_upward_ch_arcs_of_node_with_landmark(
                        popped_node, distance_to_popped_node,
                        forward_first_out, forward_head, forward_weight,
                        ch_was_forward_pushed, ch_forward_queue,
                        ch_forward_tentative_distance,
                        [&](unsigned x, unsigned pred_node, unsigned pred_arc) {
                            ch_forward_predecessor_node[x] = pred_node;
                            ch_forward_predecessor_arc[x] = pred_arc;
                        }, [&](unsigned start) {
                            return get_max_landmard_distance(start, target);
                        }, shortest_path_length, number_of_nodes_generated
                );
            }
        }

    }


    unsigned ContractionHierarchyQuery::bi_cpd(unsigned& forward_start,unsigned& backward_end,int& signal
    ){
        auto retrieve_next_move = [&](const int& source, const int& target) {
            if(source == target ){
                return -2;
            }
            const int& first_move = ch->top_n_percentage_cpd.get_first_move(source, target);
            number_of_first_move_calls ++;
            return first_move -2;
        };

        unsigned cost = 0;
        int cur_id = forward_start;
        for(;;){
            int next_move = retrieve_next_move( cur_id,backward_end);
            if(next_move == -1){
                forward_start = cur_id;
                signal = -1;
                return cost;
            }else if(next_move == -2){
                forward_start = backward_end;
                //reached target
                signal = -2;
                return cost;
            }else{
                cost = cost + ch->top_n_percentage_graph.out(cur_id,next_move).weight;

                cur_id = ch->top_n_percentage_graph.out(cur_id,next_move).node;

            }


        }


    }


    unsigned ContractionHierarchyQuery::bi_cpd_with_cache(unsigned& forward_start,unsigned& backward_end,int& signal,
                                                          std::vector<unsigned>&ch_forward_tentative_distance,
                                                          TimestampFlags&ch_was_forward_pushed,
                                                          const unsigned& forward_distance
    ){
        auto retrieve_next_move = [&](const int& source, const int& target) {
            if(source == target ){
                return -2;
            }
            const int& first_move = ch->top_n_percentage_cpd.get_first_move(source, target);
            number_of_first_move_calls ++;
            return first_move -2;
        };

        unsigned cost = 0;
        int cur_id = forward_start;
        for(;;){
            int next_move = retrieve_next_move( cur_id,backward_end);
            if(next_move == -1){
                forward_start = cur_id;
                signal = -1;
                return cost;
            }else if(next_move == -2){
                forward_start = backward_end;
                //reached target
                signal = -2;
                return cost;
            }else{
                cost = cost + ch->top_n_percentage_graph.out(cur_id,next_move).weight;

                cur_id = ch->top_n_percentage_graph.out(cur_id,next_move).node;
                unsigned cur_rank = ch->top_n_percentage_rank_mapper[cur_id];
                unsigned cur_distance = cost+forward_distance;


                if(ch_was_forward_pushed.is_set(cur_rank)){
                    if(ch_forward_tentative_distance[cur_rank] >= cur_distance){
                        ch_forward_tentative_distance[cur_rank] = cur_distance;
                    }else{
                        number_of_path_extraction --;
                        signal = -2;
                        return inf_weight;
                    }
                }else{
                    ch_was_forward_pushed.set(cur_rank);
                    ch_forward_tentative_distance[cur_rank]  = cur_distance;
                }
            }


        }


    }
    unsigned ContractionHierarchyQuery::get_bi_direction_cpd_distance(unsigned forward_start,unsigned backward_end
    ){
        bool is_forward = true;
        int signal = 0;

        unsigned cost = 0;
        for(;;){
            if(signal== -2){
                break;
            }
            if(is_forward){
                cost = cost + bi_cpd(forward_start, backward_end,signal);
                is_forward= false;
            }else{
                cost = cost + bi_cpd(backward_end, forward_start,signal);
                is_forward= true;
            }
        }
        return cost;


    }

    unsigned ContractionHierarchyQuery::get_bi_cpd_distance_with_cache(unsigned forward_start,unsigned backward_end,
                                                                       std::vector<unsigned>&ch_forward_tentative_distance,std::vector<unsigned>&ch_backward_tentative_distance,
                                                                       TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
                                                                       const unsigned& forward_distance, const unsigned& backward_distance
    ){
        bool is_forward = true;
        int signal = 0;

        unsigned cost = 0;
        unsigned bi_direction_distance =0;
        for(;;){
            if(signal== -2){
                break;
            }
            if(is_forward){
                cost = cost  +  bi_cpd_with_cache(forward_start, backward_end,signal,ch_forward_tentative_distance,ch_was_forward_pushed,forward_distance);
                is_forward= false;
//                bi_direction_distance = bi_cpd_with_cache(forward_start, backward_end,signal,ch_forward_tentative_distance,ch_was_forward_pushed,forward_distance);
//                if(bi_direction_distance == inf_weight){
//                    return inf_weight;
//                }else{
//                    cost = cost + bi_direction_distance;
//                    is_forward= false;
//                }

            }else{
                cost = cost  + bi_cpd_with_cache(backward_end, forward_start,signal,ch_backward_tentative_distance,ch_was_backward_pushed,backward_distance);
                is_forward= true;
//                        bi_direction_distance = bi_cpd_with_cache(backward_end, forward_start,signal,ch_backward_tentative_distance,ch_was_backward_pushed,backward_distance);
//
//                if(bi_direction_distance == inf_weight){
//                    return inf_weight;
//                }else{
//                    cost = cost + bi_direction_distance;
//                    is_forward= true;
//                }

            }
        }
        return cost;


    }

    unsigned ContractionHierarchyQuery::get_cpd_distance_with_cache(unsigned forward_start,unsigned backward_end
    ){
        unsigned cur_id = forward_start;
        if(forward_start == backward_end){
            return 0;
        }
        int number_of_node = 0;
        if(was_node_cached.is_set(cur_id)){
//            std::cout<<"here"<<std::endl;
            return distance_cache[cur_id];
        }
        cpd_distance_vector[number_of_node] = 0;
        cpd_path_vector[number_of_node] = cur_id;

        number_of_node ++ ;
        unsigned  final_distance  = 0 ;
        unsigned next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
        number_of_first_move_calls++;
        if(next_move== 0xFF)
            return inf_weight;
        for (;;){
//            final_distance =final_distance + ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];
            cpd_distance_vector[number_of_node] = ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];;
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id] + next_move];
            cpd_path_vector[number_of_node] = cur_id;
            if(cur_id == backward_end){
                break;
            }
            if(was_node_cached.is_set(cur_id)){
//                std::cout<<"here"<<std::endl;
                final_distance = distance_cache[cur_id];
                break;
            }
            number_of_node ++;
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
            number_of_first_move_calls++;
        }
        number_of_node --;
        while( number_of_node >= 0){
            // caching backward !
            assert( number_of_node <cpd_distance_vector.size() &&  number_of_node>=0);
            final_distance = final_distance + cpd_distance_vector[number_of_node+1];
            const unsigned& node= cpd_path_vector[number_of_node];
            was_node_cached.set(node);
            distance_cache[node] = final_distance;
            number_of_node--;
        }
        return final_distance;
    }


    unsigned ContractionHierarchyQuery::get_cpd_distance_with_cache(unsigned forward_start,
                                                                    unsigned backward_end,
                                                                    std::vector<unsigned>&ch_forward_tentative_distance,
                                                                    TimestampFlags&ch_was_forward_pushed,
                                                                    const unsigned& forward_distance
    ){
        unsigned cur_id = forward_start;
        if(forward_start == backward_end){
            return forward_distance ;
        }
        int number_of_node = 0;
        number_of_node ++ ;
        unsigned forward_cost = forward_distance;
        unsigned next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
        unsigned cur_rank;
        number_of_first_move_calls++;
        if(next_move== 0xFF)
            return inf_weight;
        for (;;){
            forward_cost = forward_cost + ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id] + next_move];
            if(cur_id == backward_end){
                break;
            }
            cur_rank = ch->top_n_percentage_rank_mapper[cur_id];
            if(ch_was_forward_pushed.is_set(cur_rank)){
                if(ch_forward_tentative_distance[cur_rank] >= forward_cost){
                    ch_forward_tentative_distance[cur_rank] = forward_cost;
                }else{
                    number_of_path_extraction --;
                    return inf_weight;
                }
            }
            else{
                ch_was_forward_pushed.set(cur_rank);
                ch_forward_tentative_distance[cur_rank]  = forward_cost;
            }

            number_of_node ++;
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
            number_of_first_move_calls++;
        }
        return forward_cost;
    }


    unsigned ContractionHierarchyQuery::get_cpd_distance_with_cache(unsigned forward_start,
                                                                    unsigned backward_end,
                                                                    std::vector<unsigned>&ch_forward_tentative_distance,
                                                                    TimestampFlags&ch_was_forward_pushed,
                                                                    const unsigned& forward_distance,
                                                                    std::vector<unsigned>& forward_node_predecessor,
                                                                    std::vector<unsigned>& forward_arc_predecessor,
                                                                    std::vector<bool>& forward_type_predecessor
    ){
        unsigned cur_id = forward_start;
        if(forward_start == backward_end){
            return forward_distance ;
        }
        unsigned forward_cost = forward_distance;
        unsigned next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
        unsigned cur_rank = ch->top_n_percentage_rank_mapper[cur_id];
        number_of_first_move_calls++;
        if(next_move== 0xFF)
            return inf_weight;
        for (;;){
            unsigned id =ch->top_n_percentage_first_out[cur_id] + next_move;
            forward_cost = forward_cost + ch->top_n_percentage_weight[id];
            unsigned pre_node = cur_rank;
            cur_id = ch->top_n_percentage_head[id];
            cur_rank = ch->top_n_percentage_rank_mapper[cur_id];
            if(ch_was_forward_pushed.is_set(cur_rank)){
                if(ch_forward_tentative_distance[cur_rank] >forward_cost){
                    forward_arc_predecessor[cur_rank] = ch->top_n_percentage_arc_id[id];
                    forward_node_predecessor[cur_rank] = pre_node;
                    forward_type_predecessor[cur_rank] = ch->top_n_percentage_is_forward[id];
                    ch_forward_tentative_distance[cur_rank] = forward_cost;
                }else if(ch_forward_tentative_distance[cur_rank] < forward_cost){
                    number_of_path_extraction --;
                    return inf_weight;
                }
            }
            else{
                forward_arc_predecessor[cur_rank] = ch->top_n_percentage_arc_id[id];
                forward_node_predecessor[cur_rank] = pre_node;
                forward_type_predecessor[cur_rank] =  ch->top_n_percentage_is_forward[id];
                ch_was_forward_pushed.set(cur_rank);
                ch_forward_tentative_distance[cur_rank]  = forward_cost;
            }

            if(cur_id == backward_end){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
            number_of_first_move_calls++;
        }
        return forward_cost;
    }


//    unsigned ContractionHierarchyQuery::get_cpd_distance_with_cache(unsigned forward_start,unsigned backward_end
//    ){
//        unsigned cur_id = forward_start;
//        if(forward_start == backward_end){
//            return 0;
//        }
//        int number_of_node = 0;
//        if(was_node_cached.is_set(cur_id)){
////            std::cout<<"here"<<std::endl;
//            return distance_cache[cur_id];
//        }
//        cpd_distance_vector[number_of_node] = 0;
//        cpd_path_vector[number_of_node] = cur_id;
//
//        number_of_node ++ ;
//        unsigned  final_distance  = 0 ;
//        unsigned next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
//        number_of_first_move_calls++;
//        if(next_move== 0xFF)
//            return inf_weight;
//        for (;;){
////            final_distance =final_distance + ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];
//            cpd_distance_vector[number_of_node] = ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];;
//            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id] + next_move];
//            cpd_path_vector[number_of_node] = cur_id;
//            if(cur_id == backward_end){
//                break;
//            }
//            if(was_node_cached.is_set(cur_id)){
////                std::cout<<"here"<<std::endl;
//                final_distance = distance_cache[cur_id];
//                break;
//            }
//            number_of_node ++;
//            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
//            number_of_first_move_calls++;
//        }
//        number_of_node --;
//        while( number_of_node >= 0){
//            // caching backward !
//            assert( number_of_node <cpd_distance_vector.size() &&  number_of_node>=0);
//            final_distance = final_distance + cpd_distance_vector[number_of_node+1];
//            const unsigned& node= cpd_path_vector[number_of_node];
//            was_node_cached.set(node);
//            distance_cache[node] = final_distance;
//            number_of_node--;
//        }
//        return final_distance;
//    }


    unsigned ContractionHierarchyQuery::get_cpd_distance_with_bi_cache(unsigned forward_start,unsigned backward_end,vector<unsigned>& cache
    ){
        unsigned cur_id = forward_start;
        if(forward_start == backward_end){
            return 0;
        }
        int number_of_node = 0;
        if(cache[cur_id] !=0){
            return cache[cur_id];
        }
        cpd_distance_vector[number_of_node] = 0;
        cpd_path_vector[number_of_node] = cur_id;

        number_of_node ++ ;
        unsigned  final_distance  = 0 ;
        unsigned next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
        number_of_first_move_calls++;
        if(next_move== 0xFF)
            return inf_weight;
        for (;;){
//            final_distance =final_distance + ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];
            cpd_distance_vector[number_of_node] = ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];;
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id] + next_move];
            cpd_path_vector[number_of_node] = cur_id;
            if(cur_id == backward_end){
                break;
            }
            if(cache[cur_id] !=0){
                final_distance = cache[cur_id];
                break;
            }
            number_of_node ++;
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
            number_of_first_move_calls++;
        }
        number_of_node --;
        while( number_of_node >= 0){
            // caching backward !
            assert( number_of_node <cpd_distance_vector.size() &&  number_of_node>=0);
            final_distance = final_distance + cpd_distance_vector[number_of_node+1];
            const unsigned& node= cpd_path_vector[number_of_node];
            cache[node]= final_distance;
            number_of_node--;
        }
        return final_distance;
    }


    unsigned ContractionHierarchyQuery::get_cpd_distance_with_bi_cache(unsigned forward_start,unsigned backward_end,vector<unsigned>& cache,TimestampFlags& flag
    ){
        unsigned cur_id = forward_start;
        if(forward_start == backward_end){
            return 0;
        }
        int number_of_node = 0;
        if(flag.is_set(cur_id)){
            return cache[cur_id];
        }
        cpd_distance_vector[number_of_node] = 0;
        cpd_path_vector[number_of_node] = cur_id;

        number_of_node ++ ;
        unsigned  final_distance  = 0 ;
        unsigned next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
        number_of_first_move_calls++;
        if(next_move== 0xFF)
            return inf_weight;
        for (;;){
//            final_distance =final_distance + ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];
            cpd_distance_vector[number_of_node] = ch->top_n_percentage_weight[ch->top_n_percentage_first_out[cur_id] + next_move];;
            cur_id = ch->top_n_percentage_head[ch->top_n_percentage_first_out[cur_id] + next_move];
            cpd_path_vector[number_of_node] = cur_id;
            if(cur_id == backward_end){
                break;
            }
            if(flag.is_set(cur_id)){
                final_distance = cache[cur_id];
                break;
            }
            number_of_node ++;
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
            number_of_first_move_calls++;
        }
        number_of_node --;
        while( number_of_node >= 0){
            // caching backward !
            assert( number_of_node <cpd_distance_vector.size() &&  number_of_node>=0);
            final_distance = final_distance + cpd_distance_vector[number_of_node+1];
            const unsigned& node= cpd_path_vector[number_of_node];
            flag.set(node);
            cache[node]= final_distance;
            number_of_node--;
        }
        return final_distance;
    }

//    unsigned ContractionHierarchyQuery::get_cpd_distance_with_cache(unsigned forward_start,unsigned backward_end,
//            std::vector<unsigned>&ch_forward_tentative_distance,
//            TimestampFlags&ch_was_forward_pushed,
//            const unsigned& forward_distance
//    ){
//        unsigned cur_id = forward_start;
//        if(forward_start == backward_end){
//            return 0;
//        }
//        unsigned cost = 0;
//        int number_of_node = 0;
//        if(was_node_cached.is_set(cur_id)){
//            return distance_cache[cur_id];
//        }
//        cpd_distance_vector[number_of_node] = 0;
//        cpd_path_vector[number_of_node] = cur_id;
//
//        number_of_node ++ ;
//        unsigned  final_distance  = 0 ;
//        unsigned cur_rank = invalid_id;
//        unsigned next_move;
//        for (;;){
//            if(cur_id == backward_end){
//                break;
//            }else{
//                next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
//            }
//            if(next_move == 1){
//                return inf_weight;
//            }else{
//                next_move = next_move  -2 ;
//                cost = cost + ch->top_n_percentage_graph.out(cur_id,next_move).weight;
//                cpd_distance_vector[number_of_node -1] = ch->top_n_percentage_graph.out(cur_id,next_move).weight;
//                cur_id = ch->top_n_percentage_graph.out(cur_id,next_move).node;
//                if(was_node_cached.is_set(cur_id)){
//                    final_distance = distance_cache[cur_id];
//                    break;
//                }
//                cur_rank = ch->top_n_percentage_rank_mapper[cur_id];
//                if(ch_was_forward_pushed.is_set(cur_rank)){
//                    if(ch_forward_tentative_distance[cur_rank] < cost+forward_distance){
//                        number_of_path_extraction --;
//                        return inf_weight;
//                    }
//                }
//                cpd_distance_vector[number_of_node] = 0;
//                cpd_path_vector[number_of_node] = cur_id;
//            }
//            number_of_node ++;
//        }
//        number_of_node =  number_of_node -1;
//        while( number_of_node >= 0){
//            // caching backward !
//            assert( number_of_node <cpd_distance_vector.size() &&  number_of_node>=0);
//            final_distance = final_distance + cpd_distance_vector[number_of_node];
//            const unsigned& node= cpd_path_vector[number_of_node];
//            was_node_cached.set(node);
//            distance_cache[node] = final_distance;
//            number_of_node--;
//        }
//        return final_distance;
//    }

    unsigned ContractionHierarchyQuery::get_cpd_distance_with_cache(unsigned forward_start,unsigned backward_end,
                                                                    std::vector<unsigned>&ch_forward_tentative_distance,std::vector<unsigned>&ch_backward_tentative_distance,
                                                                    TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed,
                                                                    const unsigned& forward_distance, const unsigned& backward_distance
    ){
        unsigned cur_id = forward_start;
        if(forward_start == backward_end){
            return 0;
        }
        unsigned cost = 0;
        int number_of_node = 0;
        if(was_node_cached.is_set(cur_id)){
            return distance_cache[cur_id];
        }
        cpd_distance_vector[number_of_node] = 0;
        cpd_path_vector[number_of_node] = cur_id;

        number_of_node ++ ;
        unsigned  final_distance  = 0 ;
        unsigned cur_rank = invalid_id;
        unsigned cur_distance = inf_weight;
        unsigned next_move;
        for (;;){
            if(cur_id == backward_end){
                break;
            }else{
                next_move = ch->top_n_percentage_cpd.get_first_move(cur_id,backward_end);
            }
            if(next_move == 1){
                return inf_weight;
            }else{
                next_move = next_move  -2 ;
                cost = cost + ch->top_n_percentage_graph.out(cur_id,next_move).weight;
                cpd_distance_vector[number_of_node -1] = ch->top_n_percentage_graph.out(cur_id,next_move).weight;
                cur_id = ch->top_n_percentage_graph.out(cur_id,next_move).node;
                if(was_node_cached.is_set(cur_id)){
                    final_distance = distance_cache[cur_id];
                    break;
                }
                cur_rank = ch->top_n_percentage_rank_mapper[cur_id];
                cur_distance = cost+forward_distance;

                if(ch_was_forward_pushed.is_set(cur_rank)){
                    if(ch_forward_tentative_distance[cur_rank] >= cur_distance){
                        ch_forward_tentative_distance[cur_rank] = cur_distance;
                    }else{
                        return inf_weight;
                    }
                }else{
                    ch_was_forward_pushed.set(cur_rank);
                    ch_forward_tentative_distance[cur_rank]  = cur_distance;
                }


                cpd_distance_vector[number_of_node] = 0;
                cpd_path_vector[number_of_node] = cur_id;
            }
            number_of_node ++;
        }
        number_of_node =  number_of_node -1;
        while( number_of_node >= 0){
            // caching backward !
            assert( number_of_node <cpd_distance_vector.size() &&  number_of_node>=0);
            final_distance = final_distance + cpd_distance_vector[number_of_node];
            const unsigned& node= cpd_path_vector[number_of_node];
            was_node_cached.set(node);
            distance_cache[node] = final_distance;

            cur_rank = ch->top_n_percentage_rank_mapper[node];
            cur_distance = final_distance +backward_distance;
            if(ch_was_backward_pushed.is_set(cur_rank)){
                if(ch_backward_tentative_distance[cur_rank] > cur_distance){
                    ch_backward_tentative_distance[cur_rank] = cur_distance;
                }
            }else{
                ch_was_backward_pushed.set(cur_rank);
                ch_backward_tentative_distance[cur_rank]  = cur_distance;
            }
            number_of_node--;
        }
        return final_distance;
    }


//    unsigned ContractionHierarchyQuery::get_max_landmard_distance(unsigned start, unsigned end){
//        unsigned max = 0 ;
//        for(int i = 0 ; i < ch->number_of_landmark; i++){
//
//            unsigned start_distance = ch->landmark_list[start+i*number_of_top_n_percentage_nodes];
//            unsigned end_distance =ch->landmark_list[end+i*number_of_top_n_percentage_nodes];
//            unsigned landmard_distance = abs((long long int)start_distance - end_distance);
//            if(max < landmard_distance){
//                max = landmard_distance;
//            }
//        }
//        return max;
//    }
//
//
//    unsigned ContractionHierarchyQuery::get_min_landmard_upperband(unsigned start, unsigned end){
//        unsigned min = inf_weight ;
//        for(int i = 0 ; i < ch->number_of_landmark; i++){
//
//            unsigned start_distance = ch->landmark_list[start+i*number_of_top_n_percentage_nodes];
//            unsigned end_distance =ch->landmark_list[end+i*number_of_top_n_percentage_nodes];
//            unsigned landmard_distance = start_distance + end_distance;
//            if(min > landmard_distance){
//                min = landmard_distance;
//            }
//        }
//        return min+1;
//    }

    unsigned ContractionHierarchyQuery::get_max_landmard_distance(unsigned start, unsigned end){
        unsigned max = 0 ;
        auto startPtr = ch->landmark_list.begin()+ start*  ch->number_of_landmark;
        auto endPtr =ch->landmark_list.begin()+ end*  ch->number_of_landmark;
        for(int i = 0 ; i < ch->number_of_landmark; i++){
            unsigned landmard_distance = abs((long long int)*(startPtr) - *(endPtr));
            if(max < landmard_distance){
                max = landmard_distance;
            }
            startPtr++;
            endPtr++;
        }
        return max;
    }


    unsigned ContractionHierarchyQuery::get_min_landmard_upperband(unsigned start, unsigned end){
        unsigned min = inf_weight ;
        auto startPtr = ch->landmark_list.begin()+ start*  ch->number_of_landmark;
        auto endPtr =ch->landmark_list.begin()+ end*  ch->number_of_landmark;
        for(int i = 0 ; i < ch->number_of_landmark; i++){
            unsigned landmard_distance = *(startPtr) + *(endPtr);
            if(min > landmard_distance){
                min = landmard_distance;
            }
            startPtr++;
            endPtr++;
        }
        return min+1;
    }

    void ContractionHierarchyQuery::update_cpd_distance(const unsigned forward_node, const std::vector<unsigned>&ch_backward_cpd_nodes, const unsigned&
    ch_number_of_backward_nodes, const unsigned distance_to_popped_node, const std::vector<unsigned>&ch_backward_tentative_distance,
                                                        unsigned& ch_best_cpd_forward,unsigned& ch_best_cpd_backward){

        was_node_cached.reset_all();
        unsigned forward_start = ch->top_n_percentage_node_mapper[forward_node] -1;
        for(int i = 0 ; i < ch_number_of_backward_nodes; i++){
            unsigned backward_node = ch_backward_cpd_nodes[i];
            unsigned backward_end = ch->top_n_percentage_node_mapper[backward_node] -1;
            unsigned distance_to_backward = ch_backward_tentative_distance[backward_node];
            unsigned max_landmark_distance = get_max_landmard_distance(ch->top_n_percentage_node_order[forward_start],ch->top_n_percentage_node_order[backward_end]);
            if(distance_to_popped_node + max_landmark_distance + distance_to_backward > shortest_path_length){
                continue;
            }
            // always calculate distance from backward to forward !
            unsigned cpd_distance = get_cpd_distance_with_cache(backward_end,forward_start);
            number_of_path_extraction ++;
            unsigned current_distance = distance_to_popped_node + cpd_distance + distance_to_backward;
            if(shortest_path_length > current_distance){
                shortest_path_meeting_node = invalid_id;
                ch_best_cpd_forward = forward_node;
                ch_best_cpd_backward = backward_node;
                shortest_path_length = current_distance;
            }
        }

    }




    void ContractionHierarchyQuery::update_cpd_distance(const unsigned forward_node, const std::vector<unsigned>&ch_backward_cpd_nodes, const unsigned&
    ch_number_of_backward_nodes, const unsigned distance_to_popped_node, unsigned& ch_best_cpd_forward,unsigned& ch_best_cpd_backward,
                                                        std::vector<unsigned>&ch_forward_tentative_distance,std::vector<unsigned>&ch_backward_tentative_distance,
                                                        TimestampFlags&ch_was_forward_pushed, TimestampFlags&ch_was_backward_pushed){
        was_node_cached.reset_all();
        unsigned forward_start = ch->top_n_percentage_node_mapper[forward_node] -1;
        for(int i = 0 ; i < ch_number_of_backward_nodes; i++){
            unsigned backward_node = ch_backward_cpd_nodes[i];
            unsigned backward_end = ch->top_n_percentage_node_mapper[backward_node] -1;
            unsigned distance_to_backward = ch_backward_tentative_distance[backward_node];
            unsigned max_landmark_distance = get_max_landmard_distance(forward_node,backward_node);
            if(distance_to_popped_node + max_landmark_distance + distance_to_backward > shortest_path_length){
                continue;
            }
            // always calculate distance from backward to forward !
            unsigned cpd_distance = get_cpd_distance_with_cache(backward_end,forward_start);
            unsigned current_distance = distance_to_popped_node + cpd_distance + distance_to_backward;
            number_of_path_extraction ++;
            if(shortest_path_length > current_distance){
                shortest_path_meeting_node = invalid_id;
                ch_best_cpd_forward = forward_node;
                ch_best_cpd_backward = backward_node;
                shortest_path_length = current_distance;
            }
        }
    }



    ContractionHierarchyQuery& ContractionHierarchyQuery::run_with_cpd(){
        assert(ch && "query object must have an attached CH");
        assert(!forward_queue.empty() && "must add at least one source before calling run");
        assert(!backward_queue.empty() && "must add at least one target before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);

        shortest_path_length = get_min_landmard_upperband(source,target);
        shortest_path_meeting_node = invalid_id;
        number_of_first_move_calls = 0;
        bool forward_next = true;
        number_of_forward_cpd_nodes=0;
        number_of_backward_cpd_nodes = 0;
        number_of_path_extraction = 0;
        number_of_nodes_generated = 0;
        number_of_nodes_expanded =0 ;
//        cpd_time = 0 ;
//        best_cpd_forward = invalid_id;
//        best_cpd_backward = invalid_id;
        for(;;) {
            bool forward_finished = false;
            if (forward_queue.empty())
                forward_finished = true;
            else if (forward_queue.peek().key >= shortest_path_length)
                forward_finished = true;

            bool backward_finished = false;
            if (backward_queue.empty())
                backward_finished = true;
            else if (backward_queue.peek().key >= shortest_path_length)
                backward_finished = true;

            if (forward_finished && backward_finished)
                break;

            if (forward_finished)
                forward_next = false;
            if (backward_finished)
                forward_next = true;

            if (forward_next) {
                forward_settle_node_with_CPD(
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        was_forward_pushed, was_backward_pushed,
                        forward_queue,
                        forward_tentative_distance, backward_tentative_distance,
                        forward_predecessor_node, forward_predecessor_arc,forward_predecessor_type,
                        forward_cpd_nodes, backward_cpd_nodes,
                        number_of_forward_cpd_nodes, number_of_backward_cpd_nodes,
                        target,true
                );
                forward_next = false;
            } else {

                forward_settle_node_with_CPD(
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        was_backward_pushed, was_forward_pushed,
                        backward_queue,
                        backward_tentative_distance, forward_tentative_distance,
                        backward_predecessor_node, backward_predecessor_arc,backward_predecessor_type,
                        backward_cpd_nodes, forward_cpd_nodes,
                        number_of_backward_cpd_nodes, number_of_forward_cpd_nodes,
                        source,false
                );

//                backward_settle_node_with_CPD(
//                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
//                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
//                        was_backward_pushed, was_forward_pushed,
//                        backward_queue,
//                        backward_tentative_distance, forward_tentative_distance,
//                        backward_predecessor_node, backward_predecessor_arc,backward_predecessor_type,
//                        forward_predecessor_node, forward_predecessor_arc,forward_predecessor_type,
//                        backward_cpd_nodes, forward_cpd_nodes,
//                        number_of_backward_cpd_nodes, number_of_forward_cpd_nodes,
//                        source,false
//                );
                forward_next = true;
            }
        }

        state = ContractionHierarchyQuery::InternalState::run;
        return *this;
    }




    ContractionHierarchyQuery& ContractionHierarchyQuery::run_with_direction_cpd(){
        assert(ch && "query object must have an attached CH");
        assert(!forward_queue.empty() && "must add at least one source before calling run");
        assert(!backward_queue.empty() && "must add at least one target before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);
        shortest_path_length = get_min_landmard_upperband(source,target);
        shortest_path_meeting_node = invalid_id;
        number_of_first_move_calls = 0;
        bool forward_next = true;
        number_of_forward_cpd_nodes=0;
        number_of_backward_cpd_nodes = 0;
        number_of_path_extraction = 0;
        number_of_nodes_generated = 0;
        number_of_nodes_expanded = 0;
        cpd_time = 0 ;
        best_cpd_forward = invalid_id;
        best_cpd_backward = invalid_id;
        vector<unsigned> forward_cpd = vector<unsigned>();
        vector<unsigned> backward_cpd = vector<unsigned>();
        for(;;){
            bool forward_finished = false;
            if(forward_queue.empty())
                forward_finished = true;
            else if(forward_queue.peek().key >= shortest_path_length)
                forward_finished = true;

            bool backward_finished = false;
            if(backward_queue.empty())
                backward_finished = true;
            else if(backward_queue.peek().key >= shortest_path_length)
                backward_finished = true;

            if(forward_finished && backward_finished)
                break;

            if(forward_finished)
                forward_next = false;
            if(backward_finished)
                forward_next = true;
            if(forward_next){
                forward_settle_node_with_direction_cpd(
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        was_forward_pushed, was_backward_pushed,
                        forward_queue,
                        forward_tentative_distance, backward_tentative_distance,
                        forward_predecessor_node, forward_predecessor_arc,target
                );
                forward_next = false;
            } else {
                forward_settle_node_with_direction_cpd(
                        ch->backward.first_out, ch->backward.head, ch->backward.weight,
                        ch->forward.first_out, ch->forward.head, ch->forward.weight,
                        was_backward_pushed, was_forward_pushed,
                        backward_queue,
                        backward_tentative_distance, forward_tentative_distance,
                        backward_predecessor_node, backward_predecessor_arc,
                        source
                );
                forward_next = true;
            }
        }
//        for(int i = 0; i < number_of_backward_cpd_nodes; i++){
//            bool exist = false;
//            for(int j = 0; j < number_of_forward_cpd_nodes; j++){
//                if(i !=j){
//                    if(backward_cpd_nodes[i] ==forward_cpd_nodes[j]  ){
//                        exist = true;
//                    }
//                }
//            }
//            if(exist){
//                std::cout<<"eroor"<<std::endl;
//            }
//        }

        state = ContractionHierarchyQuery::InternalState::run;
        return *this;
    }

    unsigned ContractionHierarchyQuery::get_used_source(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        if(shortest_path_meeting_node == invalid_id)
            return invalid_id;
        unsigned x = shortest_path_meeting_node;
        while(forward_predecessor_node[x] != invalid_id)
            x = forward_predecessor_node[x];

        return ch->order[x];
    }

    unsigned ContractionHierarchyQuery::get_used_target(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        if(shortest_path_meeting_node == invalid_id)
            return invalid_id;
        unsigned x = shortest_path_meeting_node;
        while(backward_predecessor_node[x] != invalid_id)
            x = backward_predecessor_node[x];
        return ch->order[x];
    }

    namespace{

        // OnNewInputArc has the signature
        //   on_new_arc(unsigned xy, unsigned y)
        // It is called for every arc xy of the path. y is the head of xy.
        // The source node of the path must be obtained by some other mean

        template<class OnNewInputArc>
        void unpack_forward_arc(const ContractionHierarchy&ch, unsigned arc, const OnNewInputArc&on_new_input_arc);

        template<class OnNewInputArc>
        void unpack_backward_arc(const ContractionHierarchy&ch, unsigned arc, const OnNewInputArc&on_new_input_arc);


        void unpack_forward_arc(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path);

        void unpack_backward_arc(const ContractionHierarchy&ch, unsigned arc,vector<unsigned>& path);

        void reverse_unpack_forward_arc(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path);

        void  reverse_unpack_backward_arc(const ContractionHierarchy&ch, unsigned arc,vector<unsigned>& path);

        void unpack_forward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path, unsigned path_length);

        void unpack_backward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path, unsigned path_length);

        void reverse_unpack_forward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path, unsigned path_length);

        void reverse_unpack_backward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path, unsigned path_length);


        template<class OnNewInputArc>
        void unpack_forward_arc(const ContractionHierarchy&ch, unsigned arc, const OnNewInputArc&on_new_input_arc){
            if(ch.forward.is_shortcut_an_original_arc.is_set(arc)){
                on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
            } else {
                assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
                unpack_backward_arc(ch, ch.forward.shortcut_first_arc[arc], on_new_input_arc);
                unpack_forward_arc(ch, ch.forward.shortcut_second_arc[arc], on_new_input_arc);
            }
        }
//

        template<class OnNewInputArc>
        void unpack_backward_arc(const ContractionHierarchy&ch, unsigned arc, const OnNewInputArc&on_new_input_arc){
            if(ch.backward.is_shortcut_an_original_arc.is_set(arc)){
                on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
            } else {
                assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
                unpack_backward_arc(ch, ch.backward.shortcut_first_arc[arc], on_new_input_arc);
                unpack_forward_arc(ch, ch.backward.shortcut_second_arc[arc], on_new_input_arc);
            }
        }




        void unpack_forward_arc(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path){
            if(ch.forward.is_shortcut_an_original_arc.is_set(arc)){
//                on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
                path.push_back(ch.forward.shortcut_second_arc[arc]);
            } else {
                assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
                unpack_backward_arc(ch, ch.forward.shortcut_first_arc[arc],path);
                unpack_forward_arc(ch, ch.forward.shortcut_second_arc[arc],path);
            }
        }
//

        void unpack_backward_arc(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path){
            if(ch.backward.is_shortcut_an_original_arc.is_set(arc)){
//                on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
                path.push_back(ch.backward.shortcut_second_arc[arc]);
            } else {
                assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
                unpack_backward_arc(ch, ch.backward.shortcut_first_arc[arc],path);
                unpack_forward_arc(ch, ch.backward.shortcut_second_arc[arc], path);
            }
        }


        void reverse_unpack_forward_arc(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path){
            if(ch.forward.is_shortcut_an_original_arc.is_set(arc)){
//                on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
                path.push_back(ch.forward.shortcut_second_arc[arc]);
            } else {
                assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
                reverse_unpack_forward_arc(ch, ch.forward.shortcut_second_arc[arc],path);
                reverse_unpack_backward_arc(ch, ch.forward.shortcut_first_arc[arc],path);
            }
        }
//

        void reverse_unpack_backward_arc(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path){
            if(ch.backward.is_shortcut_an_original_arc.is_set(arc)){
//                on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
                path.push_back(ch.backward.shortcut_second_arc[arc]);
            } else {
                assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
                reverse_unpack_forward_arc(ch, ch.backward.shortcut_second_arc[arc], path);
                reverse_unpack_backward_arc(ch, ch.backward.shortcut_first_arc[arc],path);
            }
        }

        void unpack_forward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path, unsigned path_length){
            if(ch.forward.is_shortcut_an_original_arc.is_set(arc)){
                path.push_back(ch.forward.shortcut_second_arc[arc]);
//            on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
            } else {
                assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
                unpack_backward_arc_prefix(ch, ch.forward.shortcut_first_arc[arc], path,path_length);
                if(path.size()== path_length){
                    return;
                }
                unpack_forward_arc_prefix(ch, ch.forward.shortcut_second_arc[arc], path,path_length);
            }
        }

        void unpack_backward_arc_prefix(const ContractionHierarchy&ch, unsigned arc,vector<unsigned>& path,unsigned path_length){
            if(ch.backward.is_shortcut_an_original_arc.is_set(arc)){
                path.push_back(ch.backward.shortcut_second_arc[arc]);
//            on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
            } else {
                assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
                unpack_backward_arc_prefix(ch, ch.backward.shortcut_first_arc[arc], path,path_length);
                if(path.size()== path_length){
                    return;
                }
                unpack_forward_arc_prefix(ch, ch.backward.shortcut_second_arc[arc], path,path_length);
            }
        }

        void reverse_unpack_forward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path,unsigned path_length){
            if(ch.forward.is_shortcut_an_original_arc.is_set(arc)){
//                on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
                path.push_back(ch.forward.shortcut_second_arc[arc]);
            } else {
                assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
                reverse_unpack_forward_arc_prefix(ch, ch.forward.shortcut_second_arc[arc],path,path_length);
                if(path.size()== path_length){
                    return;
                }
                reverse_unpack_backward_arc_prefix(ch, ch.forward.shortcut_first_arc[arc],path,path_length);
            }
        }
//

        void reverse_unpack_backward_arc_prefix(const ContractionHierarchy&ch, unsigned arc, vector<unsigned>& path,unsigned path_length){
            if(ch.backward.is_shortcut_an_original_arc.is_set(arc)){
//                on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
                path.push_back(ch.backward.shortcut_second_arc[arc]);
            } else {
                assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
                assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
                reverse_unpack_forward_arc_prefix(ch, ch.backward.shortcut_second_arc[arc], path,path_length);
                if(path.size()== path_length){
                    return;
                }
                reverse_unpack_backward_arc_prefix(ch, ch.backward.shortcut_first_arc[arc],path,path_length);
            }
        }

    }

    unsigned ContractionHierarchyQuery::get_distance() {
        assert(state == ContractionHierarchyQuery::InternalState::run);

        if(shortest_path_meeting_node == invalid_id)
            return inf_weight;
        else
            return forward_tentative_distance[shortest_path_meeting_node] + backward_tentative_distance[shortest_path_meeting_node];
    }

    unsigned ContractionHierarchyQuery::get_distance_with_cpd() {
//    assert(state == ContractionHierarchyQuery::InternalState::run);
        return shortest_path_length;
    }

    bool ContractionHierarchyQuery::check_overlapping(unsigned vaiNode){
        if(forward_predecessor_node[vaiNode] == backward_predecessor_node[vaiNode]){
            return true;
            //return false;
        }
        return false;
    }

    std::vector<unsigned>ContractionHierarchyQuery::get_arc_path(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        std::vector<unsigned>path;
        if(shortest_path_meeting_node != invalid_id)
        {
            std::vector<unsigned>up_path;
            {
                unsigned x = shortest_path_meeting_node;
                while(forward_predecessor_node[x] != invalid_id){
                    assert(was_forward_pushed.is_set(x));
                    up_path.push_back(forward_predecessor_arc[x]);
                    //up_path.push_back(find_arc_given_sorted_head(ch->forward.first_out, ch->forward.head, forward_predecessor_node[x], x));
                    x = forward_predecessor_node[x];
                }
            }
            for(unsigned i=up_path.size(); i>0; --i){
                unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(xy);});
            }
            {
                unsigned x = shortest_path_meeting_node;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    //unpack_backward_arc(*ch, find_arc_given_sorted_head(ch->backward.first_out, ch->backward.head, backward_predecessor_node[x], x), [&](unsigned xy, unsigned y){path.push_back(xy);});
                    x = backward_predecessor_node[x];
                }
            }
        }
        return path; // NVRO
    }

    std::vector<unsigned>ContractionHierarchyQuery::get_arc_path_with_cpd(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        std::vector<unsigned>path;
        if(shortest_path_meeting_node != invalid_id)
        {
            std::vector<unsigned>up_path;
            {
                unsigned x = shortest_path_meeting_node;
                while(forward_predecessor_node[x] != invalid_id){
                    assert(was_forward_pushed.is_set(x));
                    up_path.push_back(forward_predecessor_arc[x]);
                    //up_path.push_back(find_arc_given_sorted_head(ch->forward.first_out, ch->forward.head, forward_predecessor_node[x], x));
                    x = forward_predecessor_node[x];
                }
            }
            for(unsigned i=up_path.size(); i>0; --i){
                unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(xy);});
            }
            {
                unsigned x = shortest_path_meeting_node;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    //unpack_backward_arc(*ch, find_arc_given_sorted_head(ch->backward.first_out, ch->backward.head, backward_predecessor_node[x], x), [&](unsigned xy, unsigned y){path.push_back(xy);});
                    x = backward_predecessor_node[x];
                }
            }
        }else{

            if(shortest_path_length == inf_weight){
                return path; //path not found
            }else{
                // construct the path with CPD.
                std::vector<unsigned>up_path;
                {
                    unsigned x = best_cpd_forward;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        //up_path.push_back(find_arc_given_sorted_head(ch->forward.first_out, ch->forward.head, forward_predecessor_node[x], x));
                        x = forward_predecessor_node[x];
                    }
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(xy);});
                }
                vector<unsigned> arc_path;
                vector<unsigned> is_forward_path;
                get_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_backward]-1,ch->top_n_percentage_node_mapper[best_cpd_forward]-1,arc_path,is_forward_path);
                for(int i = 0 ; i < arc_path.size(); i++){
                    if(is_forward_path[i] == 1){
                        unpack_forward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    } else{
                        unpack_backward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    }
                }
                {
                    unsigned x = best_cpd_backward;
                    while (backward_predecessor_node[x] != invalid_id) {
                        assert(was_backward_pushed.is_set(x));
                        unpack_backward_arc(*ch, backward_predecessor_arc[x],
                                            [&](unsigned xy, unsigned y) { path.push_back(xy); });
                        //unpack_backward_arc(*ch, find_arc_given_sorted_head(ch->backward.first_out, ch->backward.head, backward_predecessor_node[x], x), [&](unsigned xy, unsigned y){path.push_back(xy);});
                        x = backward_predecessor_node[x];
                    }
                }
//

            }

        }
        return path; // NVRO
    }

    void ContractionHierarchyQuery::check_node_path(const vector<unsigned>& node_path, unsigned expected_cost,
                                                    const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight){
        unsigned cost = 0;
        for(int i = 0; i < node_path.size()-1; i ++){
            unsigned current =  node_path[i];
            unsigned next = node_path[i+1];
            bool exist = false;
            for(unsigned node_id = first_out[current]; node_id < first_out[current+1]; ++node_id){
                if(head[node_id]==next){
                    exist = true;
                    cost = cost +weight[node_id];
                    break;
                }
            }
            if(!exist){
                std::cout<<"path not exist"<<std::endl;
            }
        }
        if(cost != expected_cost){
            std::cout<<"path cost not match"<<std::endl;
        }
//        std::cout<<"correct"<<std::endl;
    }



    std::vector<unsigned>ContractionHierarchyQuery::get_node_path_with_bi_cpd(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        std::vector<unsigned>path;
        if(shortest_path_meeting_node != invalid_id)
        {
            std::vector<unsigned>up_path;
            {
                unsigned x = shortest_path_meeting_node;
                while(forward_predecessor_node[x] != invalid_id){
                    assert(was_forward_pushed.is_set(x));
                    up_path.push_back(forward_predecessor_arc[x]);
                    //up_path.push_back(find_arc_given_sorted_head(ch->forward.first_out, ch->forward.head, forward_predecessor_node[x], x));
                    x = forward_predecessor_node[x];
                }
            }
            for(unsigned i=up_path.size(); i>0; --i){
                unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(xy);});
            }
            {
                unsigned x = shortest_path_meeting_node;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    //unpack_backward_arc(*ch, find_arc_given_sorted_head(ch->backward.first_out, ch->backward.head, backward_predecessor_node[x], x), [&](unsigned xy, unsigned y){path.push_back(xy);});
                    x = backward_predecessor_node[x];
                }
            }
        }else{

            if(shortest_path_length == inf_weight){
                return path; //path not found
            }else{
                // construct the path with CPD.
                std::vector<unsigned>up_path;
                {
                    unsigned x = best_cpd_forward;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        //up_path.push_back(find_arc_given_sorted_head(ch->forward.first_out, ch->forward.head, forward_predecessor_node[x], x));
                        x = forward_predecessor_node[x];
                    }
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(xy);});
                }


                vector<unsigned> node_path;
                vector<unsigned> arc_path;
                vector<unsigned> is_forward_path;
//                unsigned cost = 0;
                get_bi_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_backward]-1,ch->top_n_percentage_node_mapper[best_cpd_forward]-1,node_path);


//                unsigned first_node = node_path[0];
                bool going_up = node_path[1]>node_path[0];
                for(int i = 1 ; i < node_path.size(); i++){
                    unsigned cur_node = node_path[i];
                    bool is_going_up = cur_node > node_path[i-1];
                    if(!going_up && is_going_up){
                        std::cout<< "error"<< std::endl;
                    }
                    going_up = is_going_up;
                }
//                get_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_forward],ch->top_n_percentage_node_mapper[best_cpd_backward],cpd_path,is_forward_path);
                for(int i = 0 ; i < arc_path.size(); i++){
                    if(is_forward_path[i] == 1){
                        unpack_forward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    } else{
                        unpack_backward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(xy);});
                    }
                }
                {
                    unsigned x = best_cpd_backward;
                    while (backward_predecessor_node[x] != invalid_id) {
                        assert(was_backward_pushed.is_set(x));
                        unpack_backward_arc(*ch, backward_predecessor_arc[x],
                                            [&](unsigned xy, unsigned y) { path.push_back(xy); });
                        //unpack_backward_arc(*ch, find_arc_given_sorted_head(ch->backward.first_out, ch->backward.head, backward_predecessor_node[x], x), [&](unsigned xy, unsigned y){path.push_back(xy);});
                        x = backward_predecessor_node[x];
                    }
                }
//

            }

        }
        return path; // NVRO
    }


    void ContractionHierarchyQuery::get_node_path( std::vector<unsigned>& path){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);
        
//        std::vector<unsigned>path;
        if(shortest_path_meeting_node != invalid_id)
        {
//            unsigned w = 0;
            std::vector<unsigned>up_path;
            {
                unsigned x = shortest_path_meeting_node;
                while(forward_predecessor_node[x] != invalid_id){
                    assert(was_forward_pushed.is_set(x));
                    up_path.push_back(forward_predecessor_arc[x]);
//                    w = w + ch->forward.weight[forward_predecessor_arc[x]];
//                    path.push_back(forward_predecessor_arc[x]);
                    //up_path.push_back(find_arc_given_sorted_head(ch->forward.first_out, ch->forward.head, forward_predecessor_node[x], x));
                    x = forward_predecessor_node[x];
                }
            }
            for(unsigned i=up_path.size(); i>0; --i){
                unpack_forward_arc(*ch, up_path[i-1], path);
            }
            {
                unsigned x = shortest_path_meeting_node;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], path);
//                    w = w + ch->backward.weight[backward_predecessor_arc[x]];
//                    path.push_back(forward_predecessor_arc[x]);
                    //unpack_backward_arc(*ch, find_arc_given_sorted_head(ch->backward.first_out, ch->backward.head, backward_predecessor_node[x], x), [&](unsigned xy, unsigned y){path.push_back(xy);});
                    x = backward_predecessor_node[x];
                }
            }

        }
//        return path; // NVRO
    }


    std::vector<unsigned>ContractionHierarchyQuery::get_node_path(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        std::vector<unsigned>path;
        if(shortest_path_meeting_node != invalid_id)
        {
            std::vector<unsigned>up_path;
            {
                unsigned x = shortest_path_meeting_node;
                while(forward_predecessor_node[x] != invalid_id){
                    assert(was_forward_pushed.is_set(x));
                    up_path.push_back(forward_predecessor_arc[x]);
                    x = forward_predecessor_node[x];
                }
                path.push_back(ch->order[x]);
            }
            for(unsigned i=up_path.size(); i>0; --i){
                unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(y);});
            }
            {
                unsigned x = shortest_path_meeting_node;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], [&](unsigned xy, unsigned y){path.push_back(y);});
                    x = backward_predecessor_node[x];
                }
            }
        }
        return path; // NVRO
    }



    std::vector<std::vector<unsigned>>ContractionHierarchyQuery::get_all_node_paths(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        std::vector<std::vector<unsigned>>paths;
        if(!path_meeting_nodes.empty()){
            std::set<unsigned>unique_meeting_nodes;
            for(unsigned meeting_node : path_meeting_nodes){
                unique_meeting_nodes.insert(meeting_node);
            }


            for(unsigned meeting_node :  unique_meeting_nodes){
                std::vector<unsigned> path;
                {
                    std::vector<unsigned>up_path;
                    unsigned x = meeting_node;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_node[x]);
                        x = forward_predecessor_node[x];
                    }
                    for(int i = up_path.size()-1; i >=0 ;i --){
                        path.push_back(up_path[i]);
                    }
                }

                {
                    std::vector<unsigned>down_path;
                    unsigned x = meeting_node;
                    while(backward_predecessor_node[x] != invalid_id){
                        assert(was_backward_pushed.is_set(x));
                        down_path.push_back(backward_predecessor_node[x]);
                        x = backward_predecessor_node[x];
                    }
                    for(int i = 0; i < down_path.size() ;i ++){
                        path.push_back(down_path[i]);
                    }

                }
                paths.push_back(path);


            }


        }else{
            std::cout<<"no path found"<<std::endl;
        }

        return paths;


    }

    void ContractionHierarchyQuery::get_node_path_with_cpd(std::vector<unsigned>& path){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);
        if(shortest_path_length == inf_weight){
            return;
        }
        if(shortest_path_meeting_node != invalid_id){
            if( ch->top_n_percentage_node_mapper[shortest_path_meeting_node] !=0){
                std::vector<unsigned>up_path;
                std::vector<unsigned>is_forward_up_path;
                {
                    unsigned x = shortest_path_meeting_node;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        is_forward_up_path.push_back(forward_predecessor_type[x]);
                        x = forward_predecessor_node[x];
                    }
                    path.push_back(ch->order[x]);
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    if( is_forward_up_path[i-1]){
                        unpack_forward_arc(*ch, up_path[i-1], path);
                    }else{
                        unpack_backward_arc(*ch, up_path[i-1], path);
                    }
                }

                unsigned x = shortest_path_meeting_node;
                unsigned next = backward_predecessor_node[x];

                if(next != invalid_id && next >= start_id){
                    path.erase(path.end()-1);
                    for(;;){
                        assert(was_backward_pushed.is_set(x));
                        if(next != invalid_id){
                            if( next >= start_id){
                                if(backward_predecessor_type[x]){
                                    reverse_unpack_forward_arc(*ch, backward_predecessor_arc[x], path);
                                }else{
                                    reverse_unpack_backward_arc(*ch, backward_predecessor_arc[x], path);
                                }
                                x =  next;
                            }else{
                                break;
                            }
                        }else{
                            break;
                        }
                        next = backward_predecessor_node[x];

                    }
                    path.push_back(ch->order[x]);
                }
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], path);
                    x = backward_predecessor_node[x];
                }
            }else{

                //normal_path
                std::vector<unsigned>up_path;
                {
                    unsigned x = shortest_path_meeting_node;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        x = forward_predecessor_node[x];
                    }
                    path.push_back(ch->order[x]);
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    unpack_forward_arc(*ch, up_path[i-1], path);
                }
                {
                    unsigned x = shortest_path_meeting_node;
                    while(backward_predecessor_node[x] != invalid_id){
                        assert(was_backward_pushed.is_set(x));
                        unpack_backward_arc(*ch, backward_predecessor_arc[x], path);
                        x = backward_predecessor_node[x];
                    }
                }



            }

        }

    }
//    void ContractionHierarchyQuery::get_node_path_with_cpd(std::vector<unsigned>& path){
//        assert(ch && "query object must have an attached CH");
//        assert(state == ContractionHierarchyQuery::InternalState::run);
//
//        if(shortest_path_meeting_node != invalid_id)
//        {
//            std::vector<unsigned>up_path;
//            {
//                unsigned x = shortest_path_meeting_node;
//                while(forward_predecessor_node[x] != invalid_id){
//                    assert(was_forward_pushed.is_set(x));
//                    up_path.push_back(forward_predecessor_arc[x]);
//                    x = forward_predecessor_node[x];
//                }
//                path.push_back(ch->order[x]);
//            }
//            for(unsigned i=up_path.size(); i>0; --i){
//                unpack_forward_arc(*ch, up_path[i-1], path);
//            }
//            {
//                unsigned x = shortest_path_meeting_node;
//                while(backward_predecessor_node[x] != invalid_id){
//                    assert(was_backward_pushed.is_set(x));
//                    unpack_backward_arc(*ch, backward_predecessor_arc[x], path);
//                    x = backward_predecessor_node[x];
//                }
//            }
//        }
//
//        else{
//
//            if(shortest_path_length == inf_weight){
////                return path; //path not found
//                return;
//            }else{
//                // construct the path with CPD.
//                std::vector<unsigned>up_path;
//                {
//                    unsigned x = best_cpd_forward;
//                    while(forward_predecessor_node[x] != invalid_id){
//                        assert(was_forward_pushed.is_set(x));
//                        up_path.push_back(forward_predecessor_arc[x]);
//                        x = forward_predecessor_node[x];
//                    }
//                    path.push_back(ch->order[x]);
//                }
//                for(unsigned i=up_path.size(); i>0; --i){
//                    unpack_forward_arc(*ch, up_path[i-1], path);
//                }
//                unsigned cur_id = ch->top_n_percentage_node_mapper[best_cpd_forward]-1;
//                unsigned backward_end = ch->top_n_percentage_node_mapper[best_cpd_backward]-1;
//                unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, backward_end );
//                if(next_move == 0xFF){
//                    return ;
//                }
//                unsigned next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
//                for(;;){
//                    if( ch->top_n_percentage_is_forward[next_id] == 1){
//                        unpack_forward_arc(*ch, ch->top_n_percentage_arc_id[next_id], path);
//                    } else{
//                        unpack_backward_arc(*ch, ch->top_n_percentage_arc_id[next_id], path);
//                    }
//                    cur_id = ch->top_n_percentage_head[next_id];
//                    if(cur_id == backward_end){
//                        break;
//                    }
//                    next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, backward_end);
//                    next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
//                }
//
//                unsigned x = best_cpd_backward;
//                while(backward_predecessor_node[x] != invalid_id){
//                    assert(was_backward_pushed.is_set(x));
//                    unpack_backward_arc(*ch, backward_predecessor_arc[x], path);
//                    x = backward_predecessor_node[x];
//                }
//            }
//        }
//    }




    void ContractionHierarchyQuery::get_prefix_node_path_with_cpd(std::vector<unsigned>& path,unsigned path_length){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);
        if(shortest_path_length == inf_weight){
//                return path; //path not found
            return;
        }
        if(shortest_path_meeting_node != invalid_id){
            if( ch->top_n_percentage_node_mapper[shortest_path_meeting_node] !=0){
                std::vector<unsigned>up_path;
                std::vector<unsigned>is_forward_up_path;
                {
                    unsigned x = shortest_path_meeting_node;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        is_forward_up_path.push_back(forward_predecessor_type[x]);
                        x = forward_predecessor_node[x];
                    }
                    path.push_back(ch->order[x]);
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    if( is_forward_up_path[i-1]){
                        unpack_forward_arc_prefix(*ch, up_path[i-1], path,path_length);
                    }else{
                        unpack_backward_arc_prefix(*ch, up_path[i-1], path,path_length);
                    }
                    if(path.size() == path_length){
                        return;
                    }
                }

                unsigned x = shortest_path_meeting_node;
                unsigned next = backward_predecessor_node[x];
                if(next != invalid_id && next >= start_id){
                    path.erase(path.end()-1);
                    for(;;){
                        assert(was_backward_pushed.is_set(x));
                        next = backward_predecessor_node[x];
                        if(next != invalid_id){
                            if( next >= start_id){
                                if(backward_predecessor_type[x]){
                                    reverse_unpack_forward_arc_prefix(*ch, backward_predecessor_arc[x], path,path_length);
                                }else{
                                    reverse_unpack_backward_arc_prefix(*ch, backward_predecessor_arc[x], path,path_length);
                                }
                                if(path.size() == path_length){
                                    return;
                                }
                                x =  next;
                            }else{
                                break;
                            }
                        }else{
                            break;
                        }

                    }
                    path.push_back(ch->order[x]);
                }
                if(path.size() == path_length){
                    return;
                }
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc_prefix(*ch, backward_predecessor_arc[x], path,path_length);
                    if(path.size() == path_length){
                        return;
                    }
                    x = backward_predecessor_node[x];
                }
            }else{

                //normal_path
                std::vector<unsigned>up_path;
                {
                    unsigned x = shortest_path_meeting_node;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        x = forward_predecessor_node[x];
                    }
                    path.push_back(ch->order[x]);
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    unpack_forward_arc_prefix(*ch, up_path[i-1], path, path_length);
                    if(path.size() == path_length){
                        return;
                    }
                }
                {
                    unsigned x = shortest_path_meeting_node;
                    while(backward_predecessor_node[x] != invalid_id){
                        assert(was_backward_pushed.is_set(x));
                        unpack_backward_arc_prefix(*ch, backward_predecessor_arc[x], path, path_length);
                        if(path.size() == path_length){
                            return;
                        }
                        x = backward_predecessor_node[x];
                    }
                }
            }
        }

//        if(shortest_path_meeting_node != invalid_id)
//        {
//            std::vector<unsigned>up_path;
//            {
//                unsigned x = shortest_path_meeting_node;
//                while(forward_predecessor_node[x] != invalid_id){
//                    assert(was_forward_pushed.is_set(x));
//                    up_path.push_back(forward_predecessor_arc[x]);
//                    x = forward_predecessor_node[x];
//                }
//                path.push_back(ch->order[x]);
//            }
//
//
//            for(unsigned i=up_path.size(); i>0; --i){
//                unpack_forward_arc_prefix(*ch, up_path[i-1], path, path_length);
//                if(path.size() == path_length){
//                    return;
//                }
//            }
//            {
//                unsigned x = shortest_path_meeting_node;
//                while(backward_predecessor_node[x] != invalid_id){
//                    assert(was_backward_pushed.is_set(x));
//                    unpack_backward_arc_prefix(*ch, backward_predecessor_arc[x], path, path_length);
//                    if(path.size() == path_length){
//                        return;
//                    }
//                    x = backward_predecessor_node[x];
//                }
//            }
//        }else{
//
//            if(shortest_path_length == inf_weight){
////                return path; //path not found
//                return;
//            }else{
//                // construct the path with CPD.
//                std::vector<unsigned>up_path;
//                {
//                    unsigned x = best_cpd_forward;
//                    while(forward_predecessor_node[x] != invalid_id){
//                        assert(was_forward_pushed.is_set(x));
//                        up_path.push_back(forward_predecessor_arc[x]);
//                        x = forward_predecessor_node[x];
//                    }
//                    path.push_back(ch->order[x]);
//                }
//
//                for(unsigned i=up_path.size(); i>0; --i){
//                    unpack_forward_arc_prefix(*ch, up_path[i-1], path, path_length);
//                    if(path.size() == path_length){
//                        return;
//                    }
//                }
////                vector<unsigned> arc_path;
////                vector<unsigned> is_forward_path;
////                get_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_forward]-1,ch->top_n_percentage_node_mapper[best_cpd_backward]-1,arc_path,is_forward_path);
//////                get_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_forward],ch->top_n_percentage_node_mapper[best_cpd_backward],cpd_path,is_forward_path);
////                for(int i = 0 ; i < arc_path.size(); i++){
////                    if(is_forward_path[i] == 1){
////                        unpack_forward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(y);});
////                    } else{
////                        unpack_backward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(y);});
////                    }
////                }
//
//                unsigned cur_id = ch->top_n_percentage_node_mapper[best_cpd_forward]-1;
//                unsigned backward_end = ch->top_n_percentage_node_mapper[best_cpd_backward]-1;
//                unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, backward_end );
//                if(next_move == 0xFF){
//                    return ;
//                }
//                unsigned next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
//                for(;;){
//                    if( ch->top_n_percentage_is_forward[next_id] == 1){
//                        unpack_forward_arc_prefix(*ch, ch->top_n_percentage_arc_id[next_id], path, path_length);
//                    } else{
//                        unpack_backward_arc_prefix(*ch, ch->top_n_percentage_arc_id[next_id], path, path_length);
//                    }
//                    if(path.size() == path_length){
//                        return;
//                    }
//                    cur_id = ch->top_n_percentage_head[next_id];
//                    if(cur_id == backward_end){
//                        break;
//                    }
//                    next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, backward_end);
//                    next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
//                }
//
//                unsigned x = best_cpd_backward;
//                while(backward_predecessor_node[x] != invalid_id){
//                    assert(was_backward_pushed.is_set(x));
//                    unpack_backward_arc_prefix(*ch, backward_predecessor_arc[x], path, path_length);
//                    if(path.size() == path_length){
//                        return;
//                    }
//                    x = backward_predecessor_node[x];
//
//                }
//            }
//        }
    }





    std::vector<unsigned>ContractionHierarchyQuery::get_node_path_with_cpd(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        std::vector<unsigned>path;
        if(shortest_path_meeting_node != invalid_id)
        {
            std::vector<unsigned>up_path;
            {
                unsigned x = shortest_path_meeting_node;
                while(forward_predecessor_node[x] != invalid_id){
                    assert(was_forward_pushed.is_set(x));
                    up_path.push_back(forward_predecessor_arc[x]);
                    x = forward_predecessor_node[x];
                }
                path.push_back(ch->order[x]);
            }
            for(unsigned i=up_path.size(); i>0; --i){
                unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(y);});
            }
            {
                unsigned x = shortest_path_meeting_node;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], [&](unsigned xy, unsigned y){path.push_back(y);});
                    x = backward_predecessor_node[x];
                }
            }
        }else{

            if(shortest_path_length == inf_weight){
                return path; //path not found
            }else{
                // construct the path with CPD.
                std::vector<unsigned>up_path;
                {
                    unsigned x = best_cpd_forward;
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(was_forward_pushed.is_set(x));
                        up_path.push_back(forward_predecessor_arc[x]);
                        x = forward_predecessor_node[x];
                    }
                    path.push_back(ch->order[x]);
                }
                for(unsigned i=up_path.size(); i>0; --i){
                    unpack_forward_arc(*ch, up_path[i-1], [&](unsigned xy, unsigned y){path.push_back(y);});
                }
                vector<unsigned> arc_path;
                vector<unsigned> is_forward_path;
                get_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_forward]-1,ch->top_n_percentage_node_mapper[best_cpd_backward]-1,arc_path,is_forward_path);
//                get_cpd_path(ch->top_n_percentage_node_mapper[best_cpd_forward],ch->top_n_percentage_node_mapper[best_cpd_backward],cpd_path,is_forward_path);
                for(int i = 0 ; i < arc_path.size(); i++){
                    if(is_forward_path[i] == 1){
                        unpack_forward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(y);});
                    } else{
                        unpack_backward_arc(*ch, arc_path[i], [&](unsigned xy, unsigned y){path.push_back(y);});
                    }
                }
                unsigned x = best_cpd_backward;
                while(backward_predecessor_node[x] != invalid_id){
                    assert(was_backward_pushed.is_set(x));
                    unpack_backward_arc(*ch, backward_predecessor_arc[x], [&](unsigned xy, unsigned y){path.push_back(y);});
                    x = backward_predecessor_node[x];
                }

            }
        }
        return path; // NVRO
    }

    unsigned ContractionHierarchyQuery::get_node_path_with_full_cpd_search(unsigned source, unsigned target,vector<unsigned>&path){
        number_of_first_move_calls = 0;
        target = ch->top_n_percentage_node_mapper[ch->rank[target]]-1 ;
        unsigned cur_id = ch->top_n_percentage_node_mapper[ch->rank[source]]-1;
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
        number_of_first_move_calls++;
        if(next_move == 0xFF){
            return inf_weight;
        }
        path.push_back(source);
//        vector<unsigned> arc_path;
//        vector<unsigned> is_forward;
//        unsigned cost = 0;
        unsigned next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
        for (;;){
//            arc_path.push_back( ch->top_n_percentage_arc_id[next_id]);
//            is_forward.push_back( ch->top_n_percentage_is_forward[next_id]);
//            cost = cost + ch->top_n_percentage_weight[next_id];
            cur_id = ch->top_n_percentage_head[next_id];
            if(ch->top_n_percentage_is_forward[next_id] == 1){
                unpack_forward_arc(*ch, ch->top_n_percentage_arc_id[next_id], path);
            } else{
                unpack_backward_arc(*ch, ch->top_n_percentage_arc_id[next_id], path);
            }
            if(cur_id == target){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
            next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
            number_of_first_move_calls++;
        }

        return 0;
    }



    unsigned ContractionHierarchyQuery::get_prefix_node_path_with_full_cpd_search(unsigned source, unsigned target,vector<unsigned>&path, unsigned path_length){
        number_of_first_move_calls = 0;
        target = ch->top_n_percentage_node_mapper[ch->rank[target]]-1 ;
        unsigned cur_id = ch->top_n_percentage_node_mapper[ch->rank[source]]-1;
        unsigned char next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
        number_of_first_move_calls++;
        if(next_move == 0xFF){
            return inf_weight;
        }
//        unsigned cost = 0;
        path.push_back(source);
        unsigned next_id = ch->top_n_percentage_first_out[cur_id]+next_move;

        for (;;){
//            cost = cost + ch->top_n_percentage_weight[next_id];
            cur_id = ch->top_n_percentage_head[next_id];
            if(ch->top_n_percentage_is_forward[next_id] == 1){
                unpack_forward_arc_prefix(*ch, ch->top_n_percentage_arc_id[next_id], path,path_length);
            } else{
                unpack_backward_arc_prefix(*ch, ch->top_n_percentage_arc_id[next_id], path,path_length);
            }
            if(path.size()== path_length){
                return 0 ;
            }
            if(cur_id == target){
                break;
            }
            next_move = ch->top_n_percentage_cpd.get_first_move(cur_id, target);
            next_id = ch->top_n_percentage_first_out[cur_id]+next_move;
            number_of_first_move_calls++;
        }

        return 0;
//        return cost;
    }





    ContractionHierarchyQuery&ContractionHierarchyQuery::reset_source(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::target_pinned || state == ContractionHierarchyQuery::InternalState::target_run);

        was_forward_pushed.reset_all();
        forward_queue.clear();

        state = ContractionHierarchyQuery::InternalState::target_pinned;
        return *this;
    }

    ContractionHierarchyQuery&ContractionHierarchyQuery::reset_target(){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::source_pinned || state == ContractionHierarchyQuery::InternalState::source_run);

        was_backward_pushed.reset_all();
        backward_queue.clear();

        state = ContractionHierarchyQuery::InternalState::source_pinned;
        return *this;
    }

    namespace{

        // target_list[0] ... target_list[target_count-1] are the target nodes.
        // The IDs are with repect to the CH and not with respect to the input.
        // The nodes are ordered the same way as in the input.

        // select_list[0] ... select_list[select_count-1] are the nodes reachable from a target node in the CH.
        // For a forward search, the nodes in select_list are the ones reachable in the backward CH.
        // For a backward search, the nodes in select_list are the ones reachable in the forward CH.
        // The IDs are with repect to the CH and not with respect to the input.
        // select_list is ordered decreasing by rank.

        void pin(
                const std::vector<unsigned>&external_target_list,
                const std::vector<unsigned>&external_node_to_internal_node,
                std::vector<unsigned>&target_list,
                unsigned&target_count,
                std::vector<unsigned>&select_list,
                unsigned&select_count,
                MinIDQueue&q,
                const std::vector<unsigned>&backward_first_out,
                const std::vector<unsigned>&backward_head,
                const std::vector<unsigned>&backward_weight
        ){
            target_count = external_target_list.size();

            for(unsigned i=0; i<target_count; ++i){
                unsigned t = external_target_list[i];
                t = external_node_to_internal_node[t];
                target_list[i] = t;
                if(!q.contains_id(t))
                    q.push({t,t});
            }

            select_count = 0;
            while(!q.empty()){
                auto x = q.pop().id;
                select_list[select_count++] = x;

                for(unsigned xy=backward_first_out[x]; xy < backward_first_out[x+1]; ++xy){
                    unsigned y = backward_head[xy];
                    assert(x < y);
                    if(!q.contains_id(y))
                        q.push({y,y});
                }
            }

            std::reverse(select_list.begin(), select_list.begin() + select_count);
        }

        //
        // After pinned_run is finished, there are three types of nodes:
        //  1) a source node
        //  2) a node that was reached in the forward CH
        //  3) a node that was reached in the backward CH
        // A node x is categorizes as follows:
        //  1) forward_predecessor_node[x] == invalid_id && has_forward_predecessor.is_set(x)
        //  2) forward_predecessor_node[x] != invalid_id && has_forward_predecessor.is_set(x)
        //  2) !has_forward_predecessor.is_set(x)
        //

        void pinned_run(
                std::vector<unsigned>&select_list,
                unsigned&select_count,

                TimestampFlags&has_forward_predecessor,
                MinIDQueue&forward_queue,
                std::vector<unsigned>&tentative_distance,

                std::vector<unsigned>&forward_predecessor_node,
                std::vector<unsigned>&predecessor_arc,

                const std::vector<unsigned>&forward_first_out,
                const std::vector<unsigned>&forward_head,
                const std::vector<unsigned>&forward_weight,

                const std::vector<unsigned>&backward_first_out,
                const std::vector<unsigned>&backward_head,
                const std::vector<unsigned>&backward_weight
        ){
            full_forward_search(
                    forward_first_out, forward_head, forward_weight,
                    has_forward_predecessor,
                    forward_queue,
                    tentative_distance,
                    forward_predecessor_node, predecessor_arc
            );

            for(unsigned i=0; i<select_count; ++i){
                unsigned
                        x = select_list[i],
                        dist = inf_weight,
                        pred = invalid_id;
                if(has_forward_predecessor.is_set(x))
                    dist = tentative_distance[x];

                for(unsigned xy = backward_first_out[x]; xy < backward_first_out[x+1]; ++xy){
                    unsigned y = backward_head[xy];

                    unsigned new_dist = tentative_distance[y]+backward_weight[xy];
                    if(new_dist < dist){
                        dist = new_dist;
                        pred = xy;
                    }
                }

                if(pred != invalid_id){
                    tentative_distance[x] = dist;
                    predecessor_arc[x] = pred;
                    has_forward_predecessor.reset_one(x);
                }else if(dist == inf_weight){
                    tentative_distance[x] = inf_weight;
                    predecessor_arc[x] = invalid_id;
                }
            }
        }

        void extract_distances_to_targets(
                const std::vector<unsigned>&target_list,
                unsigned target_count,
                const std::vector<unsigned>&forward_tentative_distance,
                unsigned*dist
        ){
            for(unsigned i=0; i<target_count; ++i)
                dist[i] = forward_tentative_distance[target_list[i]];
        }

        std::vector<unsigned> extract_distances_to_targets(
                const std::vector<unsigned>&target_list,
                unsigned target_count,
                const std::vector<unsigned>&forward_tentative_distance
        ){
            std::vector<unsigned>dist(target_count);
            extract_distances_to_targets(target_list, target_count, forward_tentative_distance, &dist[0]);
            return dist; // NVRO
        }
    }

    ContractionHierarchyQuery&ContractionHierarchyQuery::pin_targets(const std::vector<unsigned>&external_target_list){
        assert(ch && "query object must have an attached CH");
        assert((external_target_list.empty() || max_element_of(external_target_list) < ch->node_count()) && "node id out of bounds");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);

        pin(
                external_target_list,
                ch->rank,

                // the following 4 variables happen to be unused and of the
                // required size -> use them to avoid allocating unnecessary
                // memory. Warning: Usage must be consistent over all pinning functions
                backward_predecessor_node, many_to_many_source_or_target_count, backward_tentative_distance, shortest_path_meeting_node,

                backward_queue,

                ch->backward.first_out,
                ch->backward.head,
                ch->backward.weight
        );

        state = ContractionHierarchyQuery::InternalState::target_pinned;
        return *this;
    }

    ContractionHierarchyQuery& ContractionHierarchyQuery::pin_sources(const std::vector<unsigned>&external_source_list){
        assert(ch && "query object must have an attached CH");
        assert((external_source_list.empty() || max_element_of(external_source_list) < ch->node_count()) && "node id out of bounds");
        assert(state == ContractionHierarchyQuery::InternalState::initialized);

        pin(
                external_source_list,
                ch->rank,

                forward_predecessor_node, many_to_many_source_or_target_count, forward_tentative_distance, shortest_path_meeting_node,

                forward_queue,
                ch->forward.first_out,
                ch->forward.head,
                ch->forward.weight
        );

        state = ContractionHierarchyQuery::InternalState::source_pinned;
        return *this;
    }




    ContractionHierarchyQuery& ContractionHierarchyQuery::run_to_pinned_targets(){
        assert(ch && "query object must have an attached CH");
        assert(!forward_queue.empty() && "must add at least one source before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::target_pinned);

        pinned_run(
                backward_tentative_distance, shortest_path_meeting_node,

                was_forward_pushed,
                forward_queue,
                forward_tentative_distance,

                forward_predecessor_node, forward_predecessor_arc,

                ch->forward.first_out,
                ch->forward.head,
                ch->forward.weight,

                ch->backward.first_out,
                ch->backward.head,
                ch->backward.weight
        );

        state = ContractionHierarchyQuery::InternalState::target_run;
        return *this;
    }


    ContractionHierarchyQuery& ContractionHierarchyQuery::run_to_pinned_sources(){
        assert(ch && "query object must have an attached CH");
        assert(!backward_queue.empty() && "must add at least one target before calling run");
        assert(state == ContractionHierarchyQuery::InternalState::source_pinned);

        pinned_run(
                forward_tentative_distance, shortest_path_meeting_node,

                was_backward_pushed,
                backward_queue,
                backward_tentative_distance,

                backward_predecessor_node, backward_predecessor_arc,

                ch->backward.first_out,
                ch->backward.head,
                ch->backward.weight,

                ch->forward.first_out,
                ch->forward.head,
                ch->forward.weight
        );
        state = ContractionHierarchyQuery::InternalState::source_run;
        return *this;
    }


    ContractionHierarchyQuery& ContractionHierarchyQuery::get_distances_to_targets(unsigned*dist){
        assert(state == ContractionHierarchyQuery::InternalState::target_run);
        extract_distances_to_targets(backward_predecessor_node, many_to_many_source_or_target_count, forward_tentative_distance, dist);
        return *this;
    }

    std::vector<unsigned> ContractionHierarchyQuery::get_distances_to_targets(){
        assert(state == ContractionHierarchyQuery::InternalState::target_run);
        return extract_distances_to_targets(backward_predecessor_node, many_to_many_source_or_target_count, forward_tentative_distance);
    }


    ContractionHierarchyQuery& ContractionHierarchyQuery::get_distances_to_sources(unsigned*dist){
        assert(state == ContractionHierarchyQuery::InternalState::source_run);
        extract_distances_to_targets(forward_predecessor_node, many_to_many_source_or_target_count, backward_tentative_distance, dist);
        return *this;
    }

    std::vector<unsigned> ContractionHierarchyQuery::get_distances_to_sources(){
        assert(state == ContractionHierarchyQuery::InternalState::source_run);
        return extract_distances_to_targets(forward_predecessor_node, many_to_many_source_or_target_count, backward_tentative_distance);
    }

    namespace{
        void internal_get_used_sources_to_targets(
                const std::vector<unsigned>&target_list,
                unsigned target_count,

                const TimestampFlags&has_forward_predecessor,
                const std::vector<unsigned>&forward_predecessor_node,
                const std::vector<unsigned>&predecessor_arc,

                const std::vector<unsigned>&backward_head,

                const std::vector<unsigned>&ch_order,

                unsigned*output
        ){
            for(unsigned i=0; i<target_count; ++i){
                unsigned x = target_list[i];
                if(!has_forward_predecessor.is_set(x) && predecessor_arc[x] == invalid_id){
                    output[i] = invalid_id;
                }else{
                    while(!has_forward_predecessor.is_set(x)){
                        unsigned y = backward_head[predecessor_arc[x]];
                        assert(y > x);
                        x = y;
                    }
                    while(forward_predecessor_node[x] != invalid_id){
                        assert(has_forward_predecessor.is_set(x));
                        unsigned y = forward_predecessor_node[x];
                        assert(y < x);
                        x = y;
                    }
                    output[i] = ch_order[x];
                }
            }
        }
    }

    ContractionHierarchyQuery& ContractionHierarchyQuery::get_used_sources_to_targets(unsigned*output){
        assert(state == ContractionHierarchyQuery::InternalState::target_run);

        internal_get_used_sources_to_targets(
                backward_predecessor_node, many_to_many_source_or_target_count,

                was_forward_pushed,
                forward_predecessor_node, forward_predecessor_arc,

                ch->backward.head,
                ch->order,

                output
        );

        return *this;
    }

    std::vector<unsigned> ContractionHierarchyQuery::get_used_sources_to_targets(){
        assert(state == ContractionHierarchyQuery::InternalState::target_run);
        std::vector<unsigned>ret(many_to_many_source_or_target_count);
        get_used_sources_to_targets(&ret[0]);
        return ret; // NVRO
    }

    ContractionHierarchyQuery& ContractionHierarchyQuery::get_used_targets_to_sources(unsigned*output){
        assert(state == ContractionHierarchyQuery::InternalState::source_run);

        internal_get_used_sources_to_targets(
                forward_predecessor_node, many_to_many_source_or_target_count,

                was_backward_pushed,
                backward_predecessor_node, backward_predecessor_arc,

                ch->forward.head,
                ch->order,

                output
        );

        return *this;
    }

    std::vector<unsigned> ContractionHierarchyQuery::get_used_targets_to_sources(){
        assert(state == ContractionHierarchyQuery::InternalState::source_run);
        std::vector<unsigned>ret(many_to_many_source_or_target_count);
        get_used_targets_to_sources(&ret[0]);
        return ret; // NVRO
    }

    template struct ContractionHierarchyExtraWeight<unsigned>;
    template struct ContractionHierarchyExtraWeight<int>;
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<std::vector<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const std::vector<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<std::vector<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const std::vector<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<std::vector<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const std::vector<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<std::vector<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const std::vector<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<ContractionHierarchyExtraWeight<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const ContractionHierarchyExtraWeight<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<ContractionHierarchyExtraWeight<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const ContractionHierarchyExtraWeight<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<ContractionHierarchyExtraWeight<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const ContractionHierarchyExtraWeight<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<ContractionHierarchyExtraWeight<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const ContractionHierarchyExtraWeight<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    template ContractionHierarchyExtraWeight<unsigned>& ContractionHierarchyExtraWeight<unsigned>::reset<std::vector<unsigned>, SaturatedWeightAddition>(const ContractionHierarchy&ch, const std::vector<unsigned>&, const SaturatedWeightAddition&);
    template ContractionHierarchyExtraWeight<int>& ContractionHierarchyExtraWeight<int>::reset<std::vector<int>, SaturatedWeightAddition>(const ContractionHierarchy&ch, const std::vector<int>&, const SaturatedWeightAddition&);
    template unsigned ContractionHierarchyQuery::get_extra_weight_distance<std::vector<unsigned>,SaturatedWeightAddition>(const std::vector<unsigned>&, const SaturatedWeightAddition&);
    template int ContractionHierarchyQuery::get_extra_weight_distance<std::vector<int>,SaturatedWeightAddition>(const std::vector<int>&, const SaturatedWeightAddition&);
    template unsigned ContractionHierarchyQuery::get_extra_weight_distance<ContractionHierarchyExtraWeight<unsigned>,SaturatedWeightAddition>(const ContractionHierarchyExtraWeight<unsigned>&, const SaturatedWeightAddition&);
    template int ContractionHierarchyQuery::get_extra_weight_distance<ContractionHierarchyExtraWeight<int>,SaturatedWeightAddition>(const ContractionHierarchyExtraWeight<int>&, const SaturatedWeightAddition&);

} // namespace RoutingKit

