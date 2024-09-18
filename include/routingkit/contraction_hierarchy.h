#ifndef ROUTING_KIT_CONTRACTION_HIERARCHY_H
#define ROUTING_KIT_CONTRACTION_HIERARCHY_H

#include <routingkit/id_queue.h>
#include <routingkit/timestamp_flag.h>
#include <routingkit/bit_vector.h>
#include <routingkit/permutation.h>

#include <vector>
#include <functional>
#include <cassert>
#include <type_traits>
#include <limits.h>
#include <string>
#include <iostream>
#include <CPD.h>
#include <stack>
#include <landmark.h>
#include "my_timer.h"

namespace RoutingKit{
    class Graph{
    public:
        Graph(){}

        Graph(unsigned node_count, const std::vector<unsigned>&tail, const std::vector<unsigned>&head, const std::vector<unsigned>&weight):
                out_(node_count),
                in_(node_count),
                level_(node_count, 0){

            for(unsigned a=0; a<head.size(); ++a){
                unsigned x = tail[a];
                unsigned y = head[a];
                unsigned w = weight[a];


                if(x != y){
                    out_[x].push_back({y, w, 1, invalid_id});
                    in_[y].push_back({x, w, 1, invalid_id});
                }
            }
        }

        Graph(unsigned node_count):
                out_(node_count),
                in_(node_count),
                level_(node_count, 0){
        }

        void add_out_arc(unsigned node_id, unsigned out_id, unsigned weight,unsigned arc_id, bool is_forward){
            out_[node_id].push_back({out_id, weight, 1, invalid_id, arc_id, is_forward});
        }
        void add_in_arc(unsigned node_id, unsigned in_id, unsigned weight,unsigned arc_id, bool is_forward){
            in_[node_id].push_back({in_id, weight, 1, invalid_id, arc_id, is_forward});
        }


        void add_out_arc(unsigned node_id, unsigned out_id, unsigned weight){
            out_[node_id].push_back({out_id, weight, 1, invalid_id, 0,true});
        }
        void add_in_arc(unsigned node_id, unsigned in_id, unsigned weight){
            in_[node_id].push_back({in_id, weight, 1, invalid_id, 0,true});
        }


        int get_max_out_degree(){
            int max_degree = 0;
            for(int i = 0;i< out_.size()-1; i++){
                int degree  = out_deg(i);
                if(max_degree < degree){
                    max_degree = degree;
                }
            }
            return max_degree;

        }

        int get_max_in_degree(){
            int max_degree = 0;
            for(int i = 0;i< in_.size()-1; i++){
                int degree  = in_deg(i);
                if(max_degree < degree){
                    max_degree = degree;
                }
            }
            return max_degree;

        }

        void add_arc_or_reduce_arc_weight(unsigned x, unsigned mid_node, unsigned y, unsigned weight, unsigned hop_length){
            assert(x != y);

            assert(x < node_count());
            assert(y < node_count());

            auto reduce_arc_if_exists = [weight, hop_length, mid_node](
                    unsigned x, std::vector<Arc>&x_out,
                    unsigned y, std::vector<Arc>&y_in
            ){

                // Does arc exist?
                for(unsigned out_arc = 0; out_arc < x_out.size(); ++out_arc){
                    if(x_out[out_arc].node == y){

                        // Is the existing arc longer?
                        if(x_out[out_arc].weight <= weight)
                            return true;

                        // We need to adjust the weights
                        for(unsigned in_arc = 0; in_arc < y_in.size(); ++in_arc){
                            if(y_in[in_arc].node == x){
                                x_out[out_arc].weight = weight;
                                x_out[out_arc].hop_length = hop_length;
                                x_out[out_arc].mid_node = mid_node;
                                y_in[in_arc].weight = weight;
                                y_in[in_arc].hop_length = hop_length;
                                y_in[in_arc].mid_node = mid_node;
                                return true;
                            }
                        }
                        assert(false && "arc only exists in one direction");
                    }
                }
                return false;
            };

            if(out_[x].size() <= in_[y].size()){
                if(reduce_arc_if_exists(x, out_[x], y, in_[y]))
                    return;
            } else {
                if(reduce_arc_if_exists(y, in_[y], x, out_[x]))
                    return;
            }

            // The edges does not exist -> add the edge
            out_[x].push_back({y,weight,hop_length,mid_node});
            in_[y].push_back({x,weight,hop_length,mid_node});
        }

        void remove_all_incident_arcs(unsigned x){
            assert(x < node_count());

            auto remove_back_arcs = [&](unsigned x, const std::vector<Arc>&x_out, std::vector<std::vector<Arc>>&in){
                for(unsigned out_arc = 0; out_arc < x_out.size(); ++out_arc){
                    unsigned y = x_out[out_arc].node;
                    for(auto in_arc = in[y].begin(); ; ++in_arc){
                        assert(in_arc != in[y].end());
                        if(in_arc->node == x){
                            in[y].erase(in_arc);
                            break;
                        }
                    }
                }
            };

            remove_back_arcs(x, out_[x], in_);
            remove_back_arcs(x, in_[x], out_);

            in_[x].clear();
            out_[x].clear();
            in_[x].shrink_to_fit();
            out_[x].shrink_to_fit();
        }

        void DFSUtil(unsigned v, bool visited[])
        {
            // Mark the current node as visited and
            // print it
            visited[v] = true;

            DFSordering.push_back(v);
//        cout << v << " ";

            // Recur for all the vertices adjacent
            // to this vertex
            for (int i = 0 ; i <out_[v].size(); i++ ){
                if(!visited[out_[v][i].node]){
                    DFSUtil(out_[v][i].node, visited);
                }
            }
        }

        // DFS traversal of the vertices reachable from v.
        // It uses recursive DFSUtil()
        std::vector<unsigned>& DFS(unsigned v)
        {
            // Mark all the vertices as not visited
            bool *visited = new bool[in_.size()];
            for (int i = 0; i < in_.size(); i++)
                visited[i] = false;

            // Call the recursive helper function
            // to print DFS traversal
            DFSUtil(v, visited);
            delete[](visited);
            return DFSordering;
        }


        void DFSUtil_2(unsigned s, vector<bool> &visited)
        {
            // Create a stack for DFS
            stack<int> stack;

            // Push the current source node.
            stack.push(s);

            while (!stack.empty())
            {
                // Pop a vertex from stack and print it
                s = stack.top();
                stack.pop();

                // Stack may contain same vertex twice. So
                // we need to print the popped item only
                // if it is not visited.
                if (!visited[s])
                {
                    DFSordering.push_back(s);
                    visited[s] = true;
                }

                // Get all adjacent vertices of the popped vertex s
                // If a adjacent has not been visited, then push it
                // to the stack.
                for (int i = out_[s].size()-1 ; i >=0; i--){
                    if(!visited[out_[s][i].node]){
                        stack.push(out_[s][i].node);
                    }
                }
            }
        }


//        std::vector<unsigned>& DFS_2(unsigned v)
//        {
//            // Mark all the vertices as not visited
//            bool *visited = new bool[in_.size()];
//            for (int i = 0; i < in_.size(); i++)
//                visited[i] = false;
//
//            // Call the recursive helper function
//            // to print DFS traversal
//            DFSUtil(v, visited);
//            delete[](visited);
//            return DFSordering;
//        }

        std::vector<unsigned>& DFS_2(unsigned start)
        {
            DFSordering.clear();
            // Mark all the vertices as not visited
            vector<bool> visited(out_.size(), false);
//            DFSUtil_2(start, visited);
            for (int i = 0; i < out_.size(); i++)
                if (!visited[i])
                    DFSUtil_2(i, visited);
            return DFSordering;
        }

        void run_DFS(std::vector<unsigned>& DFS_ordering,std::vector<bool>& visited_node){
            std::cout<< "DFS fininished: "<< DFS_ordering.size()/out_.size()<<"%"<<std::endl;
            //get unvisited node;
            std::vector<unsigned> unvisited_node =  std::vector<unsigned>();
            for(unsigned i = 0; i < visited_node.size(); i++){
                if(!visited_node[i]){
                    unvisited_node.push_back(i);
                }
            }
            //select random unvisited node;
            unsigned start =rand() % unvisited_node.size();
            start = unvisited_node[start];
            // run DFS;
//            std::vector<unsigned> current_ordering3 = DFS(start);
            // avoid recursion use new implementation
            std::vector<unsigned> current_ordering = DFS_2(start);

            for(unsigned  i :current_ordering){
                DFS_ordering.push_back(i);
                //set to visit
                visited_node[i] = true;
            }
            if( DFS_ordering.size() != out_.size()){
                run_DFS( DFS_ordering,visited_node);
            }
        }

        std::vector<unsigned> get_DFS_ordering(){
            //set all node unvisited;
//            std::vector<bool> visited_node = std::vector<bool>(out_.size());
//            for(unsigned i = 0; i < out_.size(); i++ ) {
//                visited_node[i] = false;
//            }
//
//
//            std::vector<unsigned> DFS_ordering =  std::vector<unsigned>();
//            run_DFS(DFS_ordering,visited_node);
            std::vector<unsigned> DFS_ordering = DFS_2(0);
            return DFS_ordering;

        }

        struct Arc{
            unsigned node;
            unsigned weight;
            unsigned hop_length;
            unsigned mid_node;
            unsigned arc_id;
            unsigned is_forward;
        };

        unsigned out_deg(unsigned node)const{
            assert(node < node_count());
            return out_[node].size();
        }

        unsigned in_deg(unsigned node)const{
            assert(node < node_count());
            return in_[node].size();
        }

        Arc out(unsigned node, unsigned out_arc)const{
            assert(node < node_count());
            assert(out_arc < out_[node].size());
            return out_[node][out_arc];
        }

        Arc in(unsigned node, unsigned in_arc)const{
            assert(node < node_count());
            assert(in_arc < in_[node].size());
            return in_[node][in_arc];
        }

        unsigned node_count()const{
            assert(in_.size() == out_.size());
            return out_.size();
        }

        unsigned level(unsigned node)const{
            assert(node < node_count());
            return level_[node];
        }

        void raise_level(unsigned node, unsigned level){
            assert(node < node_count());
            if(level >= level_[node]){
                level_[node] = level;
            }
        }

        void graph_to_vector(std::vector<unsigned>& first_out , std::vector<unsigned>& head, std::vector<unsigned>& weight,std::vector<unsigned>& arc_id, std::vector<unsigned>& is_forward) const{
            int arc_count = 0 ;
            for(std::vector<Arc> arcs : out_){
                first_out.push_back(arc_count);
                for(Arc arc : arcs){
                    head.push_back(arc.node);
                    weight.push_back(arc.weight);
                    arc_id.push_back(arc.arc_id);
                    is_forward.push_back(arc.is_forward);
                    arc_count++;
                }
            }
            first_out.push_back(arc_count);
        }

        void vector_to_graph(std::vector<unsigned>& first_out , std::vector<unsigned>& head, std::vector<unsigned>& weight,std::vector<unsigned>& arc_id, std::vector<unsigned>& is_forward){
            out_.clear();
            out_.resize(first_out.size()-1);
            in_.clear();
            in_.resize(first_out.size()-1);
            for(unsigned i = 0 ;  i < first_out.size()-1; i ++){
                for(unsigned j = first_out[i]; j < first_out[i + 1]; ++j){
                    add_out_arc(i,head[j],weight[j],arc_id[j],is_forward[j]);
                    // since we dont need in graph, maybe fix it later
                    add_in_arc(i,head[j],weight[j],arc_id[j],is_forward[j]);
                }

            }
        }

        void vector_to_graph(std::vector<unsigned>& first_out , std::vector<unsigned>& head, std::vector<unsigned>& weight){
            out_.clear();
            out_.resize(first_out.size()-1);
            in_.clear();
            in_.resize(first_out.size()-1);
            for(unsigned i = 0 ;  i < first_out.size()-1; i ++){
                for(unsigned j = first_out[i]; j < first_out[i + 1]; ++j){
                    add_out_arc(i,head[j],weight[j]);
                    // since we dont need in graph, maybe fix it later
                    add_in_arc(i,head[j],weight[j]);
                }

            }
        }
    private:
        std::vector<std::vector<Arc>>out_, in_;
        std::vector<unsigned>level_;
        std::vector<unsigned>DFSordering;
    };

    class ContractionHierarchy{
    public:
        static const unsigned default_max_pop_count = 500;

        static ContractionHierarchy build(
                unsigned node_count, std::vector<unsigned>tail, std::vector<unsigned>head, std::vector<unsigned>weight,
                const std::function<void(std::string)>&log_message = std::function<void(std::string)>(), unsigned max_pop_count = default_max_pop_count
        );

        static ContractionHierarchy build_given_rank(
                std::vector<unsigned>rank,
                std::vector<unsigned>tail, std::vector<unsigned>head, std::vector<unsigned>weight,
                const std::function<void(std::string)>&log_message = std::function<void(std::string)>(), unsigned max_pop_count = default_max_pop_count
        );

        static ContractionHierarchy build_given_order(
                std::vector<unsigned>order,
                std::vector<unsigned>tail, std::vector<unsigned>head, std::vector<unsigned>weight,
                const std::function<void(std::string)>&log_message = std::function<void(std::string)>(), unsigned max_pop_count = default_max_pop_count
        );

        static ContractionHierarchy read(std::function<void(char*, unsigned long long)>data_source);
        static ContractionHierarchy read(std::function<void(char*, unsigned long long)>data_source, unsigned long long file_size);
        static ContractionHierarchy read(std::istream&in);
        static ContractionHierarchy read(std::istream&in, unsigned long long file_size);
        static ContractionHierarchy load_file(const std::string&file_name);

        void write(std::function<void(const char*, unsigned long long)>data_sink) const;
        void write(std::ostream&out) const;
        void save_file(const std::string&file_name) const;

        unsigned node_count()const{
            return rank.size();
        }

        struct Side{
            std::vector<unsigned>first_out;
            std::vector<unsigned>head;
            std::vector<unsigned>weight;

            BitVector is_shortcut_an_original_arc;
            std::vector<unsigned>shortcut_first_arc;  // contains input arc ID if not shortcut
            std::vector<unsigned>shortcut_second_arc; // contains input tail node ID if not shortcut
        };

        std::vector<unsigned>rank, order;
        Side forward, backward;
        std::vector<unsigned> top_n_percentage_node_order,top_n_percentage_node_mapper,top_n_percentage_node_invert_mapper,top_n_percentage_rank_mapper;
        Graph top_n_percentage_graph;
        std::vector<unsigned> top_n_percentage_first_out,top_n_percentage_head,top_n_percentage_weight,top_n_percentage_arc_id;
        std::vector<bool> top_n_percentage_is_forward;
        CPD top_n_percentage_cpd;
        CPD new_cpd;
        std::vector<unsigned> landmark_list;
//        std::vector<unsigned*> landmark_pointer;
        unsigned number_of_landmark;
        void construct_top_n_percentage_graph(double percentage);


        std::vector<unsigned int> reordering_top_n_percentage_graph_based_on_DFS(unsigned start_node_id);
        std::vector<unsigned int> reordering(unsigned start_node_id,const vector<unsigned>& DFS_ordering ,const vector<unsigned>& invert_DFS );
        void load_cpd(FILE *f);


        void load_cpd2(FILE *f);

        void
        convert_to_graph(vector<unsigned int> &first_out, vector<unsigned int> &head, vector<unsigned int> &weight);

        void construct_top_n_percentage_plain_graph(double percentage);

        vector<unsigned int> get_dfs_ordering();

        void construct_top_n_percentage_graph_without_reorder(double percentage);

        ContractionHierarchy load_file(const string &file_name, unsigned long long int &size);

        unsigned long long int get_size();

        void print_graph();
    };

    void check_contraction_hierarchy_for_errors(const ContractionHierarchy&ch);

    template<class Weight>
    struct ContractionHierarchyExtraWeight{

        ContractionHierarchyExtraWeight(){}

        template<class InputWeightContainer, class LinkFunction>
        ContractionHierarchyExtraWeight(const ContractionHierarchy&ch, const InputWeightContainer&extra_weight, const LinkFunction&link){ reset(ch, extra_weight, link); }

        template<class InputWeightContainer, class LinkFunction>
        ContractionHierarchyExtraWeight& reset(const ContractionHierarchy&ch, const InputWeightContainer&extra_weight, const LinkFunction&link);

        std::vector<Weight>forward_weight, backward_weight;
    };

    namespace detail{
        template<class T>
        using ReturnTypeWhenPassedIntOf = typename std::remove_const<typename std::remove_reference<decltype(std::declval<T>()(1))>::type>::type;

        template<class T>
        using ValueTypeOfContainer = typename std::remove_const<typename std::remove_reference<decltype(std::declval<T>()[1])>::type>::type;

        template<class T>
        struct GetExtraWeightTypeHelper{
            typedef ValueTypeOfContainer<T> type;
        };

        template<class T>
        struct GetExtraWeightTypeHelper<ContractionHierarchyExtraWeight<T>>{
            typedef T type;
        };

        template<class T>
        using GetExtraWeightType = typename GetExtraWeightTypeHelper<T>::type;
    }

    class ContractionHierarchyQuery{
    public:
        ContractionHierarchyQuery():ch(0){}
        explicit ContractionHierarchyQuery(const ContractionHierarchy&ch);

        ContractionHierarchyQuery&reset();
        ContractionHierarchyQuery&reset(const ContractionHierarchy&ch);

        ContractionHierarchyQuery&add_source(unsigned s, unsigned dist_to_s = 0);
        ContractionHierarchyQuery&add_target(unsigned t, unsigned dist_to_t = 0);

        ContractionHierarchyQuery&run();

        unsigned get_used_source();
        unsigned get_used_target();

        unsigned get_distance();
        std::vector<unsigned>get_node_path();
        std::vector<unsigned>get_arc_path();

        template<class ExtraWeight, class LinkFunction>
        detail::GetExtraWeightType<ExtraWeight> get_extra_weight_distance(const ExtraWeight&extra_weight, const LinkFunction&link);

        ContractionHierarchyQuery& reset_source();
        ContractionHierarchyQuery& pin_targets(const std::vector<unsigned>&);
        unsigned get_pinned_target_count();
        ContractionHierarchyQuery& run_to_pinned_targets();

        ContractionHierarchyQuery& get_distances_to_targets(unsigned*dist);
        std::vector<unsigned> get_distances_to_targets();

        ContractionHierarchyQuery& reset_target();
        ContractionHierarchyQuery& pin_sources(const std::vector<unsigned>&);
        unsigned get_pinned_source_count();
        ContractionHierarchyQuery& run_to_pinned_sources();

        ContractionHierarchyQuery& get_distances_to_sources(unsigned*dist);
        std::vector<unsigned> get_distances_to_sources();

        // TODO: Mirror these functions in CCH

        ContractionHierarchyQuery& get_used_sources_to_targets(unsigned*dist);
        std::vector<unsigned> get_used_sources_to_targets();

        ContractionHierarchyQuery& get_used_targets_to_sources(unsigned*dist);
        std::vector<unsigned> get_used_targets_to_sources();

        // The get_extra_weight_distances function follow a pattern.
        // The usage pattern is
        //
        //    get_extra_weight_distances_to_[target|source](extra_weight, link, [, tmp, [dist]])
        //
        // There get_extra_weight_distances_to_target is used if targets were pinned, whereas
        // get_extra_weight_distances_to_source is for pinned sources.
        //
        // There are four additional parameters. Their meaning are:
        //
        // * extra_weight: The extra weight according to which the path length should be computed.
        //   This can be a ContractionHierarchyExtraWeight<T> or a container<T>. The former is faster.
        //   container<T> is a placeholder for anything that has an operator[] that provides read-only
        //   access to the arc weight. For every arc a, the expression extra_weight[a] should
        //   give the extra weight. Using ContractionHierarchyExtraWeight<T> is faster.
        //   Typical types for container<T> are const vector<T> and const T*.
        //
        // * The extra weight does not have to be a scalar value. It can be an arbitrary structure.
        //   However, the algorithm needs to know how to concatenate the weights of two arcs.
        //   The link parameter is a functor that tells it how to do it.
        //   For travel time, the link parameter is the addition. If some of your
        //   values are inf_weight, normal addition can run into overflow problems. RoutingKit therefore
        //   provides the SaturatedWeightAddition functor that correctly handles overflows for unsigned and int weights.
        //
        //   The link function will never be provided a default constructed object, i.e.,
        //   it does not exploit that link(T(),foo) == foo.
        //
        // * tmp is a container<T>. Its size must be at least node_count. The content of
        //   the first node_count elements in undefined after this function is completed.
        //   If the parameter is omitted, a temporary buffer is allocated. If the function
        //   is called multiple times, it can be faster to only allocate one buffer.
        //
        // * dist is a container<T>. Its size is the number of pinned sources or targets.
        //   the output is written to it. If it is omited, a vector<T> is allocated and
        //   returned by the function.
        //
        //   If the source and the target nodes are equal, the extra weight length is
        //   default constructed T. If there is no path, the extra weight length is also
        //   default constructed T. For integers default constructed means 0.
        //
        // If both tmp and dist are present, the function is guarenteed to not allocate or free any memory.
        // If tmp and dist are present, and link is guarenteed to not throw, the function is guarenteed to not throw.

        template<class ExtraWeight, class LinkFunction>                                          std::vector<detail::GetExtraWeightType<ExtraWeight>> get_extra_weight_distances_to_targets(const ExtraWeight&extra_weight, const LinkFunction&link);
        template<class ExtraWeight, class LinkFunction, class TmpContainer>                      std::vector<detail::GetExtraWeightType<ExtraWeight>> get_extra_weight_distances_to_targets(const ExtraWeight&extra_weight, const LinkFunction&link, TmpContainer&tmp);
        template<class ExtraWeight, class LinkFunction, class TmpContainer, class DistContainer> ContractionHierarchyQuery&                           get_extra_weight_distances_to_targets(const ExtraWeight&extra_weight, const LinkFunction&link, TmpContainer&tmp, DistContainer&dist);
        template<class ExtraWeight, class LinkFunction>                                          std::vector<detail::GetExtraWeightType<ExtraWeight>> get_extra_weight_distances_to_sources(const ExtraWeight&extra_weight, const LinkFunction&link);
        template<class ExtraWeight, class LinkFunction, class TmpContainer>                      std::vector<detail::GetExtraWeightType<ExtraWeight>> get_extra_weight_distances_to_sources(const ExtraWeight&extra_weight, const LinkFunction&link, TmpContainer&tmp);
        template<class ExtraWeight, class LinkFunction, class TmpContainer, class DistContainer> ContractionHierarchyQuery&                           get_extra_weight_distances_to_sources(const ExtraWeight&extra_weight, const LinkFunction&link, TmpContainer&tmp, DistContainer&dist);

//private:
        const ContractionHierarchy*ch;

        TimestampFlags was_forward_pushed, was_backward_pushed;
        MinIDQueue forward_queue, backward_queue;
        std::vector<unsigned>forward_tentative_distance, backward_tentative_distance;
        std::vector<bool>forward_predecessor_type, backward_predecessor_type;
        std::vector<unsigned>forward_predecessor_node, backward_predecessor_node;
        std::vector<unsigned>forward_predecessor_arc, backward_predecessor_arc;
        std::set<unsigned> path_meeting_nodes;
        unsigned shortest_path_meeting_node;
        unsigned many_to_many_source_or_target_count;
        unsigned shortest_path_length;



        std::vector<vector<unsigned>>forward_caching;
        std::vector<vector<unsigned>>backward_caching;
        std::vector<TimestampFlags>forward_time_flag;
        std::vector<TimestampFlags>backward_time_flag;

        TimestampFlags was_node_cached;
        std::vector<unsigned>distance_cache;
        unsigned best_cpd_forward, best_cpd_backward;
        std::vector<unsigned>cpd_distance_vector,cpd_path_vector;
        vector<unsigned> forward_cpd_nodes, backward_cpd_nodes;
        unsigned number_of_forward_cpd_nodes,number_of_backward_cpd_nodes;
        unsigned number_of_first_move_calls;
        unsigned number_of_path_extraction;
        unsigned number_of_top_n_percentage_nodes;
        unsigned number_of_nodes_generated;
        unsigned number_of_nodes_expanded;
        unsigned start_id;
        double cpd_time;
        CPD top_n_percentage_cpd;
        std::vector<unsigned> top_n_percentage_node_mapper;
        Graph top_n_percentage_graph;
        my_timer timer1 = my_timer();
        std::vector<double> lat,lon;
        unsigned source,target;
        std::vector<landmark> landmark_list;

        enum class InternalState:unsigned{
            initialized,
            run,
            source_pinned,
            source_run,
            target_pinned,
            target_run
        }state;

        unsigned int get_distance_with_cpd();

        bool check_overlapping(unsigned vaiNode);
        void print_bal(unsigned vaiNode);

        vector<unsigned int> get_node_path_with_cpd();
        void get_node_path_with_cpd( vector<unsigned int>& path);
        std::vector<unsigned>get_arc_path_with_cpd();

        std::vector<unsigned>get_node_path_with_bi_cpd();

        unsigned int get_cpd_distance_with_cache(unsigned int forward_start, unsigned int backward_end);

        ContractionHierarchyQuery &run_with_cpd();


        void update_cpd_distance(unsigned int forward_node, const vector<unsigned int> &ch_backward_cpd_nodes,
                                 const unsigned int &ch_number_of_backward_nodes, unsigned int distance_to_popped_node,
                                 const vector<unsigned int> &ch_backward_tentative_distance,
                                 unsigned int &ch_best_cpd_forward,
                                 unsigned int &ch_best_cpd_backward);

        unsigned int get_bi_direction_cpd_distance(unsigned int forward_start, unsigned int backward_end);

        unsigned int bi_cpd(unsigned int &forward_start, unsigned int &backward_end, int &signal);

        unsigned int bi_cpd_path (unsigned& forward_start,unsigned& backward_end,
                                             vector<unsigned>&node_path, int & signal);


        unsigned int
        get_bi_cpd_path(unsigned int forward_start, unsigned int backward_end, vector<unsigned int> &node_path);

        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed,
                                     const TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     const vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                     unsigned int &ch_best_cpd_forward,
                                     unsigned int &ch_best_cpd_backward, unsigned int &target,
                                     TimestampFlags &ch_was_forward_cpd_stored,
                                     TimestampFlags &ch_was_backward_cpd_stored);

        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed,
                                     const TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     const vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                     unsigned int &ch_best_cpd_forward,
                                     unsigned int &ch_best_cpd_backward, unsigned int &target);

        unsigned int get_max_landmard_distance(unsigned int start, unsigned int end);



        void update_cpd_distance(const unsigned int forward_node, const vector<unsigned int> &ch_backward_cpd_nodes,
                                 const unsigned int &ch_number_of_backward_nodes,
                                 const unsigned int distance_to_popped_node,
                                 unsigned int &ch_best_cpd_forward, unsigned int &ch_best_cpd_backward,
                                 vector<unsigned int> &ch_forward_tentative_distance,
                                 vector<unsigned int> &ch_backward_tentative_distance,
                                 TimestampFlags &ch_was_forward_pushed,
                                 TimestampFlags &ch_was_backward_pushed);

        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                     unsigned int &ch_best_cpd_forward,
                                     unsigned int &ch_best_cpd_backward, unsigned int &target);

        void forward_settle_node_until_find_CPD(const vector<unsigned int> &forward_first_out,
                                                const vector<unsigned int> &forward_head,
                                                const vector<unsigned int> &forward_weight,
                                                const vector<unsigned int> &backward_first_out,
                                                const vector<unsigned int> &backward_head,
                                                const vector<unsigned int> &backward_weight,
                                                TimestampFlags &ch_was_forward_pushed,
                                                TimestampFlags &ch_was_backward_pushed,
                                                MinIDQueue &ch_forward_queue,
                                                vector<unsigned int> &ch_forward_tentative_distance,
                                                vector<unsigned int> &ch_backward_tentative_distance,
                                                vector<unsigned int> &ch_forward_predecessor_node,
                                                vector<unsigned int> &ch_forward_predecessor_arc,
                                                vector<unsigned int> &ch_forward_cpd_nodes,
                                                const vector<unsigned int> &ch_backward_cpd_nodes,
                                                unsigned int &ch_number_of_forward_cpd_nodes,
                                                const unsigned int &ch_number_of_backward_cpd_nodes,
                                                unsigned int &ch_best_cpd_forward, unsigned int &ch_best_cpd_backward,
                                                unsigned int &target);

        unsigned int get_cpd_distance_with_cache(unsigned int forward_start, unsigned int backward_end,
                                                 vector<unsigned int> &ch_forward_tentative_distance,
                                                 TimestampFlags &ch_was_forward_pushed,
                                                 const unsigned int &forward_distance);

        unsigned int bi_cpd_with_cache(unsigned int &forward_start, unsigned int &backward_end, int &signal,
                                       vector<unsigned int> &ch_forward_tentative_distance,
                                       TimestampFlags &ch_was_forward_pushed, const unsigned int &forward_distance);

        unsigned int get_bi_cpd_distance_with_cache(unsigned int forward_start, unsigned int backward_end,
                                                    vector<unsigned int> &ch_forward_tentative_distance,
                                                    vector<unsigned int> &ch_backward_tentative_distance,
                                                    TimestampFlags &ch_was_forward_pushed,
                                                    TimestampFlags &ch_was_backward_pushed,
                                                    const unsigned int &forward_distance,
                                                    const unsigned int &backward_distance);

        unsigned int get_cpd_distance_with_cache(unsigned int forward_start, unsigned int backward_end,
                                                 vector<unsigned int> &ch_forward_tentative_distance,
                                                 vector<unsigned int> &ch_backward_tentative_distance,
                                                 TimestampFlags &ch_was_forward_pushed,
                                                 TimestampFlags &ch_was_backward_pushed,
                                                 const unsigned int &forward_distance,
                                                 const unsigned int &backward_distance);

        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, MinIDQueue &ch_backward_queue,
                                     vector<unsigned int> &ch_forward_tentative_distance,
                                     vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                     unsigned int &ch_best_cpd_forward,
                                     unsigned int &ch_best_cpd_backward, unsigned int &target);

        void update_cpd_distance(const unsigned int forward_node, const vector<unsigned int> &ch_backward_cpd_nodes,
                                 const unsigned int &ch_number_of_backward_nodes,
                                 const unsigned int distance_to_popped_node,
                                 unsigned int &ch_best_cpd_forward, unsigned int &ch_best_cpd_backward,
                                 vector<unsigned int> &ch_forward_tentative_distance,
                                 vector<unsigned int> &ch_backward_tentative_distance,
                                 TimestampFlags &ch_was_forward_pushed,
                                 TimestampFlags &ch_was_backward_pushed, MinIDQueue &ch_forward_queue,
                                 MinIDQueue &ch_backward_queue);

        unsigned int get_bi_cpd_distance_with_cache(unsigned int forward_start, unsigned int backward_end,
                                                    vector<unsigned int> &ch_forward_tentative_distance,
                                                    vector<unsigned int> &ch_backward_tentative_distance,
                                                    TimestampFlags &ch_was_forward_pushed,
                                                    TimestampFlags &ch_was_backward_pushed,
                                                    const unsigned int &forward_distance,
                                                    const unsigned int &backward_distance,
                                                    MinIDQueue &ch_forward_queue, MinIDQueue &ch_backward_queue);

        unsigned int bi_cpd_with_cache(unsigned int &forward_start, unsigned int &backward_end, int &signal,
                                       vector<unsigned int> &ch_forward_tentative_distance,
                                       TimestampFlags &ch_was_forward_pushed, const unsigned int &forward_distance,
                                       MinIDQueue &forward_queue);

        void
        forward_settle_node(const vector<unsigned int> &forward_first_out, const vector<unsigned int> &forward_head,
                            const vector<unsigned int> &forward_weight, const vector<unsigned int> &backward_first_out,
                            const vector<unsigned int> &backward_head, const vector<unsigned int> &backward_weight,
                            TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                            MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                            vector<unsigned int> &ch_backward_tentative_distance,
                            vector<unsigned int> &ch_forward_predecessor_node,
                            vector<unsigned int> &ch_forward_predecessor_arc, const unsigned int target);

        void
        forward_settle_node_with_landmark(const vector<unsigned int> &forward_first_out, const vector<unsigned int> &forward_head,
                            const vector<unsigned int> &forward_weight, const vector<unsigned int> &backward_first_out,
                            const vector<unsigned int> &backward_head, const vector<unsigned int> &backward_weight,
                            TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                            MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                            vector<unsigned int> &ch_backward_tentative_distance,
                            vector<unsigned int> &ch_forward_predecessor_node,
                            vector<unsigned int> &ch_forward_predecessor_arc, const unsigned int target);
        unsigned int get_min_landmard_upperband(unsigned int start, unsigned int end);

        ContractionHierarchyQuery &run_with_direction_cpd();

        void forward_settle_node_with_direction_cpd(const vector<unsigned int> &forward_first_out,
                                                    const vector<unsigned int> &forward_head,
                                                    const vector<unsigned int> &forward_weight,
                                                    const vector<unsigned int> &backward_first_out,
                                                    const vector<unsigned int> &backward_head,
                                                    const vector<unsigned int> &backward_weight,
                                                    TimestampFlags &ch_was_forward_pushed,
                                                    TimestampFlags &ch_was_backward_pushed,
                                                    MinIDQueue &ch_forward_queue,
                                                    vector<unsigned int> &ch_forward_tentative_distance,
                                                    vector<unsigned int> &ch_backward_tentative_distance,
                                                    vector<unsigned int> &ch_forward_predecessor_node,
                                                    vector<unsigned int> &ch_forward_predecessor_arc,
                                                    const unsigned int target);

        ContractionHierarchyQuery &run_to_get_all_paths();

        void
        forward_settle_node_all_path(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc, const unsigned int target);

        vector <std::vector<unsigned int>> get_all_node_paths();

        unsigned int run_cpd_search(unsigned int source, unsigned int target, vector<unsigned int> &path);

        ContractionHierarchyQuery &run_with_landmark(const unsigned shortest_path_distance);

        void run_cpd_get_path(unsigned int source, unsigned int target, vector<unsigned int> &path);

        unsigned int run_full_cpd_search(unsigned int source, unsigned int target, vector<unsigned int> &path);

        bool get_cpd_path(unsigned int forward_start, unsigned int backward_end, vector<unsigned int> &arc_path,
                          vector<unsigned int> &is_forward);

        void check_node_path(const vector<unsigned int>& node_path, unsigned int expected_cost,const vector<unsigned int>&first_out,const vector<unsigned int>& head, const vector<unsigned int>& weight);

        unsigned int
        get_node_path_with_full_cpd_search(unsigned int source, unsigned int target, vector<unsigned int> &path);

        unsigned int run_full_cpd_search(unsigned int source, unsigned int target);

        void get_node_path(vector<unsigned int> &path);

        unsigned int
        get_prefix_node_path_with_full_cpd_search(unsigned int source, unsigned int target, vector<unsigned int> &path,
                                                  unsigned int path_lenght);

//        void
//        get_prefix_node_path_with_full_cpd_search(unsigned int source, unsigned int target, vector<unsigned int> &path,
//                                                  unsigned int path_lenght);


        void get_prefix_node_path_with_cpd(vector<unsigned int> &path, unsigned int path_length);

        void forward_settle_node_with_CPD_bi_caching(const vector<unsigned int> &forward_first_out,
                                                     const vector<unsigned int> &forward_head,
                                                     const vector<unsigned int> &forward_weight,
                                                     const vector<unsigned int> &backward_first_out,
                                                     const vector<unsigned int> &backward_head,
                                                     const vector<unsigned int> &backward_weight,
                                                     TimestampFlags &ch_was_forward_pushed,
                                                     TimestampFlags &ch_was_backward_pushed,
                                                     MinIDQueue &ch_forward_queue,
                                                     vector<unsigned int> &ch_forward_tentative_distance,
                                                     vector<unsigned int> &ch_backward_tentative_distance,
                                                     vector<unsigned int> &ch_forward_predecessor_node,
                                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                                     unsigned int &ch_best_cpd_forward,
                                                     unsigned int &ch_best_cpd_backward,
                                                     vector <vector<unsigned int>> &f_cache,
                                                     vector <vector<unsigned int>> &b_cache,
                                                     vector <TimestampFlags> &f_flag,
                                                     vector <TimestampFlags> &b_flag);

        unsigned int
        get_cpd_distance_with_bi_cache(unsigned int forward_start, unsigned int backward_end,
                                       vector<unsigned int> &cache);

        void forward_settle_node_with_CPD_bi_caching(const vector<unsigned int> &forward_first_out,
                                                     const vector<unsigned int> &forward_head,
                                                     const vector<unsigned int> &forward_weight,
                                                     const vector<unsigned int> &backward_first_out,
                                                     const vector<unsigned int> &backward_head,
                                                     const vector<unsigned int> &backward_weight,
                                                     TimestampFlags &ch_was_forward_pushed,
                                                     TimestampFlags &ch_was_backward_pushed,
                                                     MinIDQueue &ch_forward_queue,
                                                     vector<unsigned int> &ch_forward_tentative_distance,
                                                     vector<unsigned int> &ch_backward_tentative_distance,
                                                     vector<unsigned int> &ch_forward_predecessor_node,
                                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                                     unsigned int &ch_best_cpd_forward,
                                                     unsigned int &ch_best_cpd_backward,
                                                     unsigned int target, vector <vector<unsigned int>> &f_cache,
                                                     vector <vector<unsigned int>> &b_cache,
                                                     vector <TimestampFlags> &f_flag,
                                                     vector <TimestampFlags> &b_flag);

        unsigned int
        get_cpd_distance_with_bi_cache(unsigned int forward_start, unsigned int backward_end,
                                       vector<unsigned int> &cache,
                                       TimestampFlags &flag);


        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<unsigned int> &ch_forward_predecessor_type,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                     unsigned int &ch_best_cpd_forward,
                                     unsigned int &ch_best_cpd_backward, unsigned int &target);

        unsigned int get_cpd_distance_with_cache(unsigned int forward_start, unsigned int backward_end,
                                                 vector<unsigned int> &ch_forward_tentative_distance,
                                                 TimestampFlags &ch_was_forward_pushed,
                                                 const unsigned int &forward_distance,
                                                 vector<unsigned int> &forward_node_predecessor,
                                                 vector<unsigned int> &forward_arc_predecessor,
                                                 vector<bool> &forward_type_predecessor);

        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<bool> &ch_forward_predecessor_type,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes,
                                     unsigned int &ch_best_cpd_forward,
                                     unsigned int &ch_best_cpd_backward, unsigned int &target, const bool is_forward);

        void
        forward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                     const vector<unsigned int> &forward_head,
                                     const vector<unsigned int> &forward_weight,
                                     const vector<unsigned int> &backward_first_out,
                                     const vector<unsigned int> &backward_head,
                                     const vector<unsigned int> &backward_weight,
                                     TimestampFlags &ch_was_forward_pushed, TimestampFlags &ch_was_backward_pushed,
                                     MinIDQueue &ch_forward_queue, vector<unsigned int> &ch_forward_tentative_distance,
                                     vector<unsigned int> &ch_backward_tentative_distance,
                                     vector<unsigned int> &ch_forward_predecessor_node,
                                     vector<unsigned int> &ch_forward_predecessor_arc,
                                     vector<bool> &ch_forward_predecessor_type,
                                     vector<unsigned int> &ch_forward_cpd_nodes,
                                     const vector<unsigned int> &ch_backward_cpd_nodes,
                                     unsigned int &ch_number_of_forward_cpd_nodes,
                                     const unsigned int &ch_number_of_backward_cpd_nodes, unsigned int &target,
                                     const bool is_forward);

        void
        backward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                      const vector<unsigned int> &forward_head,
                                      const vector<unsigned int> &forward_weight,
                                      const vector<unsigned int> &backward_first_out,
                                      const vector<unsigned int> &backward_head,
                                      const vector<unsigned int> &backward_weight,
                                      TimestampFlags &ch_was_forward_pushed,
                                      TimestampFlags &ch_was_backward_pushed, MinIDQueue &ch_forward_queue,
                                      vector<unsigned int> &ch_forward_tentative_distance,
                                      vector<unsigned int> &ch_backward_tentative_distance,
                                      vector<unsigned int> &ch_forward_predecessor_node,
                                      vector<unsigned int> &ch_forward_predecessor_arc,
                                      vector<bool> &ch_forward_predecessor_type,
                                      vector<unsigned int> &ch_forward_cpd_nodes,
                                      const vector<unsigned int> &ch_backward_cpd_nodes,
                                      unsigned int &ch_number_of_forward_cpd_nodes,
                                      const unsigned int &ch_number_of_backward_cpd_nodes, unsigned int &target,
                                      const bool is_forward);

        void
        backward_settle_node_with_CPD(const vector<unsigned int> &forward_first_out,
                                      const vector<unsigned int> &forward_head,
                                      const vector<unsigned int> &forward_weight,
                                      const vector<unsigned int> &backward_first_out,
                                      const vector<unsigned int> &backward_head,
                                      const vector<unsigned int> &backward_weight,
                                      TimestampFlags &ch_was_forward_pushed,
                                      TimestampFlags &ch_was_backward_pushed, MinIDQueue &ch_forward_queue,
                                      vector<unsigned int> &ch_forward_tentative_distance,
                                      vector<unsigned int> &ch_backward_tentative_distance,
                                      vector<unsigned int> &ch_forward_predecessor_node,
                                      vector<unsigned int> &ch_forward_predecessor_arc,
                                      vector<bool> &ch_forward_predecessor_type,
                                      vector<unsigned int> &ch_backward_predecessor_node,
                                      vector<unsigned int> &ch_backward_predecessor_arc,
                                      vector<bool> &ch_backward_predecessor_type,
                                      vector<unsigned int> &ch_forward_cpd_nodes,
                                      const vector<unsigned int> &ch_backward_cpd_nodes,
                                      unsigned int &ch_number_of_forward_cpd_nodes,
                                      const unsigned int &ch_number_of_backward_cpd_nodes, unsigned int &target,
                                      const bool is_forward);
    };

    struct SaturatedWeightAddition{
        unsigned operator()(unsigned l, unsigned r)const;
        int operator()(int l, int r)const;
    };

// ------ Template & inline implementations; no more interface descriptions beyond this line -------

    inline
    unsigned ContractionHierarchyQuery::get_pinned_target_count(){
        assert(state == ContractionHierarchyQuery::InternalState::target_run || state == ContractionHierarchyQuery::InternalState::target_pinned);
        return many_to_many_source_or_target_count;
    }

    inline
    unsigned ContractionHierarchyQuery::get_pinned_source_count(){
        assert(state == ContractionHierarchyQuery::InternalState::source_run || state == ContractionHierarchyQuery::InternalState::source_pinned);
        return many_to_many_source_or_target_count;
    }

    template<class ExtraWeight, class LinkFunction>
    std::vector<detail::GetExtraWeightType<ExtraWeight>> ContractionHierarchyQuery::get_extra_weight_distances_to_targets(
            const ExtraWeight&extra_weight,
            const LinkFunction&link
    ){
        std::vector<detail::GetExtraWeightType<ExtraWeight>>tmp(ch->node_count()), dist(get_pinned_target_count());
        get_extra_weight_distances_to_targets(extra_weight, link, tmp, dist);
        return dist; // NVRO
    }

    template<class ExtraWeight, class LinkFunction, class TmpContainer>
    std::vector<detail::GetExtraWeightType<ExtraWeight>> ContractionHierarchyQuery::get_extra_weight_distances_to_targets(
            const ExtraWeight&extra_weight,
            const LinkFunction&link,
            TmpContainer&tmp
    ){
        std::vector<detail::GetExtraWeightType<ExtraWeight>>dist(get_pinned_target_count());
        get_extra_weight_distances_to_targets(extra_weight, link, tmp, dist);
        return dist; // NVRO
    }


    template<class ExtraWeight, class LinkFunction>
    std::vector<detail::GetExtraWeightType<ExtraWeight>> ContractionHierarchyQuery::get_extra_weight_distances_to_sources(
            const ExtraWeight&extra_weight,
            const LinkFunction&link
    ){
        std::vector<detail::GetExtraWeightType<ExtraWeight>>tmp(ch->node_count()), dist(get_pinned_source_count());
        get_extra_weight_distances_to_sources(extra_weight, link, tmp, dist);
        return dist; // NVRO
    }

    template<class ExtraWeight, class LinkFunction, class TmpContainer>
    std::vector<detail::GetExtraWeightType<ExtraWeight>> ContractionHierarchyQuery::get_extra_weight_distances_to_sources(
            const ExtraWeight&extra_weight,
            const LinkFunction&link,
            TmpContainer&tmp
    ){
        std::vector<detail::GetExtraWeightType<ExtraWeight>>dist(get_pinned_source_count());
        get_extra_weight_distances_to_sources(extra_weight, link, tmp, dist);
        return dist; // NVRO
    }

    inline
    unsigned SaturatedWeightAddition::operator()(unsigned l, unsigned r)const{
        assert(l <= inf_weight && "unsigned weight must not be larger than inf_weight");
        assert(r <= inf_weight && "unsigned weight must not be larger than inf_weight");
        if(l >= inf_weight-r)
            return inf_weight;
        else
            return l+r;
    }

    inline
    int SaturatedWeightAddition::operator()(int l, int r)const{
        static_assert(inf_weight == INT_MAX, "this function assumes that inf_weight is INT_MAX");
        if(l > 0){
            if (r > INT_MAX - l){
                return INT_MAX;
            }
        }else if(r < INT_MIN - l){
            return INT_MIN;
        }

        return l + r;
    }

    namespace detail{
        template<class LinkFunction>
        struct InverseLinkFunction{
            explicit InverseLinkFunction(const LinkFunction&link):link(link){}

            template<class L, class R>
            auto operator()(L&&l, R&&r)const
            -> decltype(std::declval<LinkFunction>()(std::forward<R>(r), std::forward<L>(l)))
            {
                return link(std::forward<R>(r), std::forward<L>(l));
            }

            const LinkFunction&link;
        };

        template<class LinkFunction>
        InverseLinkFunction<LinkFunction>inverse_link_function(const LinkFunction&link){
            return InverseLinkFunction<LinkFunction>(link);
        }

        template<class InputWeightContainer, class LinkFunction>
        struct ShortcutWeights{
            typedef typename std::remove_reference<decltype(std::declval<InputWeightContainer>()[0])>::type Weight;

            const InputWeightContainer&input_weight;
            const LinkFunction&link;

            const ContractionHierarchy&ch;

            ShortcutWeights(const InputWeightContainer&input_weight, const LinkFunction&link, const ContractionHierarchy&ch):
                    input_weight(input_weight), link(link), ch(ch){}

            Weight get_forward_weight(unsigned a)const{
                assert(a < ch.forward.is_shortcut_an_original_arc.size());
                if(ch.forward.is_shortcut_an_original_arc.is_set(a)){
                    return input_weight[ch.forward.shortcut_first_arc[a]];
                }else{
                    return link(get_backward_weight(ch.forward.shortcut_first_arc[a]), get_forward_weight(ch.forward.shortcut_second_arc[a]));
                }
            }

            Weight get_backward_weight(unsigned a)const{
                assert(a < ch.backward.is_shortcut_an_original_arc.size());
                if(ch.backward.is_shortcut_an_original_arc.is_set(a)){
                    return input_weight[ch.backward.shortcut_first_arc[a]];
                }else{
                    return link(get_backward_weight(ch.backward.shortcut_first_arc[a]), get_forward_weight(ch.backward.shortcut_second_arc[a]));
                }
            }
        };

        template<class WeightT, class LinkFunction>
        struct ShortcutWeights<ContractionHierarchyExtraWeight<WeightT>, LinkFunction>{

            typedef WeightT Weight;

            const ContractionHierarchyExtraWeight<WeightT>&extra_weight;

            explicit ShortcutWeights(const ContractionHierarchyExtraWeight<WeightT>&extra_weight, const LinkFunction&, const ContractionHierarchy&):
                    extra_weight(extra_weight){}

            const Weight&get_forward_weight(unsigned a)const{
                return extra_weight.forward_weight[a];
            }

            const Weight&get_backward_weight(unsigned a)const{
                return extra_weight.backward_weight[a];
            }
        };

        template<class ExtraWeight, class LinkFunction>
        ShortcutWeights<ExtraWeight, LinkFunction> make_shortcut_weights(const ExtraWeight&extra_weight, const LinkFunction&link, const ContractionHierarchy&ch){
            return ShortcutWeights<ExtraWeight, LinkFunction>{extra_weight, link, ch};
        }

        template<class ShortcutWeights>
        struct InvertShorcutWeights{
            typedef typename ShortcutWeights::Weight Weight;

            const ShortcutWeights&shortcut_weights;

            explicit InvertShorcutWeights(const ShortcutWeights&shortcut_weights):shortcut_weights(shortcut_weights){}

            auto get_forward_weight(unsigned a)const->decltype(std::declval<ShortcutWeights>().get_backward_weight(a)){
                return shortcut_weights.get_backward_weight(a);
            }

            auto get_backward_weight(unsigned a)const->decltype(std::declval<ShortcutWeights>().get_forward_weight(a)){
                return shortcut_weights.get_forward_weight(a);
            }
        };

        template<class ShortcutWeights>
        InvertShorcutWeights<ShortcutWeights>inverse_shortcut_weights(const ShortcutWeights&shortcut_weights){
            return InvertShorcutWeights<ShortcutWeights>(shortcut_weights);
        }

        template<class GetForwardWeight, class LinkFunction>
        ReturnTypeWhenPassedIntOf<GetForwardWeight> get_extra_weight_up_distance(
                unsigned shortest_path_meeting_node,
                const std::vector<unsigned>&forward_predecessor_node,
                const std::vector<unsigned>&forward_predecessor_arc,
                const GetForwardWeight&get_forward_extra_weight,
                const LinkFunction&link
        ){
            using Weight = ReturnTypeWhenPassedIntOf<GetForwardWeight>;

            unsigned x = shortest_path_meeting_node;
            assert(forward_predecessor_node[x] != invalid_id);
            Weight ret = get_forward_extra_weight(forward_predecessor_arc[x]);
            x = forward_predecessor_node[x];

            while(forward_predecessor_node[x] != invalid_id){
                ret = link(get_forward_extra_weight(forward_predecessor_arc[x]), ret);
                x = forward_predecessor_node[x];
            }
            return ret;
        }


        template<class ShortcutWeights, class LinkFunction>
        typename ShortcutWeights::Weight internal_get_extra_weight_distance(
                const ShortcutWeights&shortcut_weights,
                const LinkFunction&link,
                unsigned shortest_path_meeting_node,
                const std::vector<unsigned>&forward_predecessor_node,
                const std::vector<unsigned>&forward_predecessor_arc,
                const std::vector<unsigned>&backward_predecessor_node,
                const std::vector<unsigned>&backward_predecessor_arc
        ){
            using Weight = typename ShortcutWeights::Weight;

            if(shortest_path_meeting_node == invalid_id)
                return Weight{};


            auto inverted_link = detail::inverse_link_function(link);

            bool has_up_part = forward_predecessor_node[shortest_path_meeting_node] != invalid_id;
            bool has_down_part = backward_predecessor_node[shortest_path_meeting_node] != invalid_id;

            // This if-then-else chain avoids the need for a neutral element with respect to link

            if(!has_up_part && !has_down_part) {
                return Weight{};
            } else if(has_up_part && has_down_part) {
                return link(
                        detail::get_extra_weight_up_distance(
                                shortest_path_meeting_node,
                                forward_predecessor_node,
                                forward_predecessor_arc,
                                [&](unsigned a)->decltype(shortcut_weights.get_forward_weight(a)){return shortcut_weights.get_forward_weight(a);},
                                link
                        ),
                        detail::get_extra_weight_up_distance(
                                shortest_path_meeting_node,
                                backward_predecessor_node,
                                backward_predecessor_arc,
                                [&](unsigned a)->decltype(shortcut_weights.get_backward_weight(a)){return shortcut_weights.get_backward_weight(a);},
                                inverted_link
                        )
                );
            } else if(has_up_part) {
                return detail::get_extra_weight_up_distance(
                        shortest_path_meeting_node,
                        forward_predecessor_node,
                        forward_predecessor_arc,
                        [&](unsigned a)->decltype(shortcut_weights.get_forward_weight(a)){return shortcut_weights.get_forward_weight(a);},
                        link
                );
            } else {
                return detail::get_extra_weight_up_distance(
                        shortest_path_meeting_node,
                        backward_predecessor_node,
                        backward_predecessor_arc,
                        [&](unsigned a)->decltype(shortcut_weights.get_backward_weight(a)){return shortcut_weights.get_backward_weight(a);},
                        inverted_link
                );
            }
        }

    }

    template<class ExtraWeight, class LinkFunction>
    detail::GetExtraWeightType<ExtraWeight> ContractionHierarchyQuery::get_extra_weight_distance(
            const ExtraWeight&extra_weight,
            const LinkFunction&link){
        assert(ch && "query object must have an attached CH");
        assert(state == ContractionHierarchyQuery::InternalState::run);

        auto shortcut_weight = detail::make_shortcut_weights(extra_weight, link, *ch);

        return detail::internal_get_extra_weight_distance(
                shortcut_weight, link,
                shortest_path_meeting_node,
                forward_predecessor_node, forward_predecessor_arc,
                backward_predecessor_node, backward_predecessor_arc
        );
    }

    template<class Weight> template<class InputWeightContainer, class LinkFunction>
    ContractionHierarchyExtraWeight<Weight>& ContractionHierarchyExtraWeight<Weight>::reset(const ContractionHierarchy&ch, const InputWeightContainer&input_extra_weight, const LinkFunction&link){
        const unsigned node_count = ch.node_count();

        forward_weight.resize(ch.forward.weight.size());
        backward_weight.resize(ch.backward.weight.size());

        for(unsigned x=0; x<node_count; ++x){
            for(unsigned xy=ch.forward.first_out[x]; xy<ch.forward.first_out[x+1]; ++xy){
                if(ch.forward.is_shortcut_an_original_arc.is_set(xy)){
                    forward_weight[xy] = input_extra_weight[ch.forward.shortcut_first_arc[xy]];
                } else {
                    forward_weight[xy] = link(
                            backward_weight[ch.forward.shortcut_first_arc[xy]],
                            forward_weight[ch.forward.shortcut_second_arc[xy]]
                    );
                }
            }
            for(unsigned xy=ch.backward.first_out[x]; xy<ch.backward.first_out[x+1]; ++xy){
                if(ch.backward.is_shortcut_an_original_arc.is_set(xy)){
                    backward_weight[xy] = input_extra_weight[ch.backward.shortcut_first_arc[xy]];
                } else {
                    backward_weight[xy] = link(
                            backward_weight[ch.backward.shortcut_first_arc[xy]],
                            forward_weight[ch.backward.shortcut_second_arc[xy]]
                    );
                }
            }
        }
        return *this;
    }


    namespace detail{
        template<class LinkFunction, class ExtraWeight, class TmpContainer, class DistContainer>
        void extract_distances_to_targets(
                const std::vector<unsigned>&target_list,
                unsigned target_count,

                const TimestampFlags&has_forward_predecessor,
                const std::vector<unsigned>&forward_predecessor_node,
                const std::vector<unsigned>&predecessor_arc,

                const ExtraWeight&extra_weight,
                const std::vector<unsigned>&forward_first_out,
                const std::vector<unsigned>&forward_head,
                const std::vector<unsigned>&backward_first_out,
                const std::vector<unsigned>&backward_head,

                TmpContainer&source_to_node_distance,
                TimestampFlags&has_source_to_node_distance,

                DistContainer&output,

                std::vector<unsigned>&stack,

                const LinkFunction&link
        ){
            using Weight = typename ExtraWeight::Weight;

            has_source_to_node_distance.reset_all();

            unsigned nodes_on_stack_with_forward_precessor_count = 0;
            unsigned stack_size = 0;

            auto push = [&](unsigned x){
                assert(stack_size < stack.size());
                stack[stack_size++] = x;
            };

            auto pop = [&]{
                assert(stack_size != 0);
                return stack[--stack_size];
            };

            auto single_forward_step_expand_distance_to_node = [&](unsigned x){
                assert(!has_source_to_node_distance.is_set(x));
                assert(has_forward_predecessor.is_set(x));
                assert(forward_predecessor_node[x] != invalid_id);

                unsigned a = predecessor_arc[x];
                unsigned p = forward_predecessor_node[x];

                if(has_source_to_node_distance.is_set(p))
                    source_to_node_distance[x] = link(source_to_node_distance[p], extra_weight.get_forward_weight(a));
                else
                    source_to_node_distance[x] = extra_weight.get_forward_weight(a); // p is source
                has_source_to_node_distance.set(x);
            };

            auto single_backward_step_expand_distance_to_node = [&](unsigned x){
                assert(!has_source_to_node_distance.is_set(x));
                assert(!has_forward_predecessor.is_set(x));

                unsigned a = predecessor_arc[x];
                unsigned p = backward_head[a];

                if(has_source_to_node_distance.is_set(p))
                    source_to_node_distance[x] = link(source_to_node_distance[p], extra_weight.get_backward_weight(a));
                else
                    source_to_node_distance[x] = extra_weight.get_backward_weight(a); // p is source

                has_source_to_node_distance.set(x);
            };

            auto push_non_reached_nodes = [&](unsigned x){
                assert(!has_source_to_node_distance.is_set(x));
                while(!has_forward_predecessor.is_set(x)){
                    push(x);
                    unsigned y = backward_head[predecessor_arc[x]];
                    assert(y > x);
                    x = y;
                    if(has_source_to_node_distance.is_set(x))
                        return;
                }
                while(forward_predecessor_node[x] != invalid_id){
                    assert(has_forward_predecessor.is_set(x));
                    push(x);
                    ++nodes_on_stack_with_forward_precessor_count;
                    unsigned y = forward_predecessor_node[x];
                    assert(y < x);
                    x = y;
                    if(has_source_to_node_distance.is_set(x))
                        return;
                }
                source_to_node_distance[x] = Weight{}; // x is source node
            };

            auto expand_distances_from_source_to_node = [&](unsigned x){
                if(!has_source_to_node_distance.is_set(x)){
                    push_non_reached_nodes(x);
                    while(stack_size != 0){
                        unsigned y = pop();
                        if(nodes_on_stack_with_forward_precessor_count != 0){
                            --nodes_on_stack_with_forward_precessor_count;
                            single_forward_step_expand_distance_to_node(y);
                        }else{
                            single_backward_step_expand_distance_to_node(y);
                        }
                    }
                }
            };

            for(unsigned i=0; i<target_count; ++i){
                unsigned t = target_list[i];
                assert(t != invalid_id);

                if(!has_forward_predecessor.is_set(t) && predecessor_arc[t] == invalid_id){
                    output[i] = Weight{}; // t is not reachable
                }else{
                    expand_distances_from_source_to_node(t);
                    output[i] = source_to_node_distance[t];
                }
            }
        }
    }

    template<class ExtraWeight, class LinkFunction, class TmpContainer, class DistContainer>
    ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets(
            const ExtraWeight&extra_weight,
            const LinkFunction&link,
            TmpContainer&tmp,
            DistContainer&dist
    ){
        assert(state == ContractionHierarchyQuery::InternalState::target_run);

        auto shortcut_weight = detail::make_shortcut_weights(extra_weight, link, *ch);

        detail::extract_distances_to_targets(
                backward_predecessor_node, many_to_many_source_or_target_count,

                was_forward_pushed,
                forward_predecessor_node, forward_predecessor_arc,

                shortcut_weight,
                ch->forward.first_out,
                ch->forward.head,
                ch->backward.first_out,
                ch->backward.head,

                tmp,
                was_backward_pushed,

                dist,

                backward_predecessor_arc,

                link
        );

        return *this;
    }

    template<class ExtraWeight, class LinkFunction, class TmpContainer, class DistContainer>
    ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources(
            const ExtraWeight&extra_weight,
            const LinkFunction&link,
            TmpContainer&tmp,
            DistContainer&dist
    ){
        assert(state == ContractionHierarchyQuery::InternalState::source_run);

        auto inverted_link = detail::inverse_link_function(link);

        auto shortcut_weight = detail::make_shortcut_weights(extra_weight, link, *ch);
        auto inverted_shortcut_weight = detail::inverse_shortcut_weights(shortcut_weight);

        detail::extract_distances_to_targets(
                forward_predecessor_node, many_to_many_source_or_target_count,

                was_backward_pushed,
                backward_predecessor_node, backward_predecessor_arc,

                inverted_shortcut_weight,
                ch->backward.first_out,
                ch->backward.head,
                ch->forward.first_out,
                ch->forward.head,

                tmp,
                was_forward_pushed,

                dist,

                forward_predecessor_arc,

                inverted_link
        );

        return *this;
    }

    extern template struct ContractionHierarchyExtraWeight<unsigned>;
    extern template struct ContractionHierarchyExtraWeight<int>;
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<std::vector<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const std::vector<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<std::vector<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const std::vector<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<std::vector<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const std::vector<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<std::vector<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const std::vector<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<ContractionHierarchyExtraWeight<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const ContractionHierarchyExtraWeight<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<ContractionHierarchyExtraWeight<int>, SaturatedWeightAddition, std::vector<int>, std::vector<int>>(const ContractionHierarchyExtraWeight<int>&, const SaturatedWeightAddition&, std::vector<int>&, std::vector<int>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_targets<ContractionHierarchyExtraWeight<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const ContractionHierarchyExtraWeight<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    extern template ContractionHierarchyQuery& ContractionHierarchyQuery::get_extra_weight_distances_to_sources<ContractionHierarchyExtraWeight<unsigned>, SaturatedWeightAddition, std::vector<unsigned>, std::vector<unsigned>>(const ContractionHierarchyExtraWeight<unsigned>&, const SaturatedWeightAddition&, std::vector<unsigned>&, std::vector<unsigned>&);
    extern template ContractionHierarchyExtraWeight<unsigned>& ContractionHierarchyExtraWeight<unsigned>::reset<std::vector<unsigned>, SaturatedWeightAddition>(const ContractionHierarchy&ch, const std::vector<unsigned>&, const SaturatedWeightAddition&);
    extern template ContractionHierarchyExtraWeight<int>& ContractionHierarchyExtraWeight<int>::reset<std::vector<int>, SaturatedWeightAddition>(const ContractionHierarchy&ch, const std::vector<int>&, const SaturatedWeightAddition&);
    extern template unsigned ContractionHierarchyQuery::get_extra_weight_distance<std::vector<unsigned>,SaturatedWeightAddition>(const std::vector<unsigned>&, const SaturatedWeightAddition&);
    extern template int ContractionHierarchyQuery::get_extra_weight_distance<std::vector<int>,SaturatedWeightAddition>(const std::vector<int>&, const SaturatedWeightAddition&);
    extern template unsigned ContractionHierarchyQuery::get_extra_weight_distance<ContractionHierarchyExtraWeight<unsigned>,SaturatedWeightAddition>(const ContractionHierarchyExtraWeight<unsigned>&, const SaturatedWeightAddition&);
    extern template int ContractionHierarchyQuery::get_extra_weight_distance<ContractionHierarchyExtraWeight<int>,SaturatedWeightAddition>(const ContractionHierarchyExtraWeight<int>&, const SaturatedWeightAddition&);

} // namespace RoutingKit

#endif