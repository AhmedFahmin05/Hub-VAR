#include <routingkit/vector_io.h>
#include <routingkit/my_timer.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include "CPD.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <routingkit/geo_dist.h>
#include <iomanip>
#include <random>
#include "../src/graph.h"
#include "../src/coverage_ordering_path.h"
#include "../src/performance_metrics.h"
#include "verify.h"
#include <cassert>

using namespace RoutingKit;
using namespace std;

bool comp(tuple<int, int> &a, tuple<int, int> &b )
{
    return (get<0>(a) == get<0>(b)) && (get<1>(a) == get<1>(b));
}


void normalise_vector(vector<double> &val){
    double max = *max_element(val.begin(), val.end());
    double min = *min_element(val.begin(), val.end());
    //cout << max << " " << min << endl;
    for(unsigned i = 0; i < val.size(); i++){
        //cout << val[i] << " ";
        val[i] = (val[i] - min) / (max - min);
        //cout << val[i] << endl;
    }

}

int main(int argc, char*argv[]){

    try{
        string file_location;
        string first_out_file;
        string head_file;
        string weight_file;
        string query_file;
        int num_of_alternative_paths;
        int qcount;

        if(argc != 4){
            cerr << argv[0] << "file_location query_count  " << endl;
            return 1;
        }else{
            file_location = argv[1];
            first_out_file = file_location + ".first";
            head_file = file_location + ".head";
            weight_file = file_location + ".weight";
            query_file = argv[2];
            num_of_alternative_paths = strtol(argv[3], NULL, 10);
        }


        cout << "Loading graph ... " << flush;

        vector<unsigned>first_out = load_vector<unsigned>(first_out_file);
        vector<unsigned>head = load_vector<unsigned>(head_file);
        vector<unsigned>weight = load_vector<unsigned>(weight_file);

        cout << "done" << endl;


        cout << "Validity tests ... " << flush;
        check_if_graph_is_valid(first_out, head);
        cout << "done" << endl;


        std::vector<tuple <int, int>> query_points;
        std::ifstream ifs(query_file);
        unsigned a , b;
        ifs >> std::ws;
        while (ifs.good()){
            ifs >> a >> b ;
            query_points.push_back(std::make_tuple(a , b));
            ifs >> std::ws;
        }
        ifs.close();
        cout << "Number of queries : " << query_points.size() << endl;


        cout << "Loading Hub Labels" << endl;
        string graph_name = file_location + ".txt";
        string label_name = file_location + ".label";
        string order_name = file_location + ".order";
        const char* graphFileName=  graph_name.c_str();
        const char* labelFileName = label_name.c_str();
        const char* orderFileName = order_name.c_str();

        unsigned node_count = first_out.size()-1;

        PLabel lab;
        lab.load_labels(labelFileName);
        ifstream order_ifs(orderFileName);
        vector<NodeID> rank(node_count);
        vector<NodeID> inv(node_count);
        for(int i = 0; i < node_count; ++i){
            NodeID tv;
            order_ifs >> tv;
            rank[tv] = i;
            inv[i] = tv;
        }
        cout << "Done" << endl;
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(50); // seed the generator
        std::uniform_int_distribution<> distr(0, node_count); // define the range

        unsigned from,to;
        vector<pair<unsigned , unsigned >> via_nodes;
//        unsigned dis = lab.generating_via_nodes(from, to, rank, inv, via_nodes);
//        cout << dis << " " << INF_WEIGHT <<  endl;
//        cout << via_nodes.size() << endl;


        /*
        int j = 0;
        // The loop while queue is empty
        cout << "Total via nodes : " << via_nodes.size() << endl;
        while(!via_nodes.empty())
        {
            // Prints and pops the element in the top of the queue
            pair<int,int> top=via_nodes.top();
            via_nodes.pop();
            cout<<"("<<top.first<<","<<top.second<<")"<<endl;
            if(top.second < double(dis) * 1.5){
                j++;
            }
        }
        vector<int> path = vector<int>();
        int  cost = lab.query_path(from, to, rank, inv, path);
        cout << cost << endl;
        cout << "The path is";
        for(auto x:path)
            cout << " " << x;
        cout << endl;
        cout << "Total via nodes after distance puring : " << j << endl;*/

        std::vector<int> via_numbers;
        std::vector<int> via_numbers_after_distance_puring;
        std::vector<int> via_numbers_after_detour_filtering;
        std::vector<int> via_numbers_after_similarity_filtering;

        std::vector<double> s_score;
        std::vector<double> b_s_score;
        std::vector<double> lo_scores;
        vector<vector< tuple < vector<unsigned> , unsigned, double >>> alternative_Paths;


        int total = 0;
        int c = 0;
        for (auto&& point: query_points)
        {
            std::tie(from, to) = point;
            vector<unsigned> meeting_points;
            meeting_points.push_back(0);
            c++;
            cout << "Query Count : " << c << endl;
            cout << "From: " << from << endl;
            cout << "To: " << to << endl;
            //cout << from << "  " << to << endl;
            if(lab.query(from, to) == INF_WEIGHT){
                cout << "Unreachable" << endl;
                continue;
            }
            vector< tuple < vector<unsigned> , unsigned, double >> a_paths;
            vector<int> s_path(0);
            unsigned shortest_distance = lab.query_path(from, to, rank, inv, s_path);
            a_paths.push_back(make_tuple(vector<unsigned>(s_path.begin(),s_path.end()), shortest_distance, 0));


            vector< tuple<unsigned, unsigned, unsigned> > via_nodes_union;
            unsigned dis1 = lab.generating_via_nodes_union(from, to, rank, inv, via_nodes_union, shortest_distance);
            vector<unsigned> nodes;
            for(auto x:via_nodes_union) {
                unsigned v,d, s;
                std::tie(v, d, s) = x;
                if( s == 0 ){
                    if(lab.query(v, to) + d  < double(shortest_distance) * 1.5){
                        nodes.push_back(v);
                    }
                } else if (s == 1){
                    if(lab.query(from, v) + d  < double(shortest_distance) * 1.5){
                        nodes.push_back(v);
                    }
                } else{
                   nodes.push_back(v);
                }
            }

            int j = 0;
            int k = 0;

            for(auto x:nodes)
            {

                j++;
                unsigned v_node = x;
                int distance1, distance2;

                vector<int> path1(0);
                distance1 = lab.query_path(from, v_node, rank, inv, path1);

                vector<int> path2(0);
                distance2 = lab.query_path(v_node, to, rank, inv, path2);

                if(path1[path1.size() - 2] == path2[1]){
                    continue;
                }

                vector<unsigned> path(path1.begin(), path1.end());
                path.insert(path.end(), path2.begin() + 1, path2.end());

                // Only checking with the optimal
                /*int flag = 1;
                double max = 0;
                double sim_score = similarity_check(get<0>(a_paths[0]), path, get<1>(a_paths[0]), distance1 + distance2, first_out, head, weight);
                if(sim_score > 0.50){
                    flag = 0;
                }

                if(flag == 1){
                    a_paths.push_back(make_tuple(path, distance1 + distance2, sim_score));
                    meeting_points.push_back(path1.size());
                    // cout << "******" << endl;
                    // cout << path1.size() << endl;
                }*/

                int flag = 1;
                double max = 0;
                for(auto p:a_paths){
                    double sim_score = similarity_check(get<0>(p), path, get<1>(p), distance1 + distance2, first_out, head, weight);
                    if(sim_score > 0.50){
                        flag = 0;
                        break;
                    }
                    if(sim_score > max){
                        max = sim_score;
                    }
                }
                if(flag == 1){
                    a_paths.push_back(make_tuple(path, distance1 + distance2, max));
                    meeting_points.push_back(path1.size());
                }
            }
            if(a_paths.size() >= num_of_alternative_paths){
                std::sort(begin(a_paths) , end(a_paths), [](const tuple < vector<unsigned > , unsigned, double >& a,
                                                            const tuple < vector<unsigned > , unsigned, double >& b) -> bool
                {
                    return std::get<1>(a) < std::get<1>(b);
                });

                cout << "Length of " << num_of_alternative_paths << " alternative paths are :" << endl;
                for(int i = 0; i < num_of_alternative_paths; i++){
                    cout << get<1>(a_paths[i]) << endl;
                }

//                double max = 0;
//                vector<double> similarity_scores(a_paths.size());
//                vector<double> bounded_stretch_scores(a_paths.size());
//                vector<double> local_optimality_scores(a_paths.size());
//                vector<double> length_scores(a_paths.size());
//                similarity_scores[0] = 0;
//                bounded_stretch_scores[0] = 1;
//                local_optimality_scores[0] = 1;
//                length_scores[0] = 0;
//                double bs_score, lo_score;
//
//                for (int i = 1; i < a_paths.size(); i++) {
//                    similarity_scores[i] = get<2>(a_paths[i]);
//                    std::tie(bs_score, lo_score) = get_Local_Optimality_and_Bounded_Stretch_Optimal(get<0>(a_paths[i]), meeting_points[i], lab, first_out, head, weight, shortest_distance, 10);
//
//                    bounded_stretch_scores[i] = bs_score;
//                    local_optimality_scores[i] = lo_score;
//                    length_scores[i] = (double) (get<1>(a_paths[i]) - get<1>(a_paths[0])) / get<1>(a_paths[0]);
//
//                }
//
//                normalise_vector(similarity_scores);
//                normalise_vector(bounded_stretch_scores);
//                normalise_vector(local_optimality_scores);
//                normalise_vector(length_scores);
/*
                vector< tuple < vector<unsigned> , unsigned, double >> result_paths;
                result_paths.push_back(make_tuple(get<0>(a_paths[0]), get<1>(a_paths[0]), get<2>(a_paths[0])));
                a_paths.erase(a_paths.begin() + 0);
                do {
                    for (int i = 0; i < a_paths.size(); i++) {
                        double normalised_score =
                                -bounded_stretch_scores[i] - similarity_scores[i] + local_optimality_scores[i] -
                                length_scores[i];
                        // cout << similarity_scores[i] << "  " << bounded_stretch_scores[i] << "  " << local_optimality_scores[i] << "  " << length_scores[i] << " " << normalised_score << endl;
                        get<2>(a_paths[i]) = normalised_score;
                    }

                    std::sort(begin(a_paths) , end(a_paths), [](const tuple < vector<unsigned > , unsigned, double >& a,
                                                                const tuple < vector<unsigned > , unsigned, double >& b) -> bool
                    {
                        return std::get<2>(a) > std::get<2>(b);
                    });

                    result_paths.push_back(make_tuple(get<0>(a_paths[0]), get<1>(a_paths[0]), get<2>(a_paths[0])));
                    a_paths.erase(a_paths.begin() + 0);

                    if(result_paths.size() == num_of_alternative_paths) break;

                    for (int i = 0; i < a_paths.size(); i++) {
                        double max = 0;
                        for(int j = 0 < j < result_paths.size(); j++){
                            double score = similarity_check(get<0>(a_paths[i]), get<0>(result_paths[j]), get<1>(a_paths[i]), get<1>(result_paths[j]), first_out, head, weight);
                            // cout << score << endl;
                            if (score > max){
                                max = score;
                            }
                        }

                        double normalised_score =
                                -bounded_stretch_scores[i] - similarity_scores[i] + local_optimality_scores[i] -
                                length_scores[i];
                        // cout << similarity_scores[i] << "  " << bounded_stretch_scores[i] << "  " << local_optimality_scores[i] << "  " << length_scores[i] << " " << normalised_score << endl;
                        get<2>(a_paths[i]) = normalised_score;
                    }

                } while(result_paths.size() <= num_of_alternative_paths);*/


//                for (int i = 0; i < a_paths.size(); i++) {
//                    double normalised_score =
//                            -bounded_stretch_scores[i] - similarity_scores[i] + local_optimality_scores[i] -
//                            length_scores[i];
//                        // cout << similarity_scores[i] << "  " << bounded_stretch_scores[i] << "  " << local_optimality_scores[i] << "  " << length_scores[i] << " " << normalised_score << endl;
//                        get<2>(a_paths[i]) = normalised_score;
//                }
//                alternative_Paths.push_back(a_paths);
            }
            else{
                cout << "The algorithm couldn't return " << num_of_alternative_paths << " paths" << endl;
            }

        }

//        cout << "Number of queries returned k path : " << total << "  out of  " << c << endl;
//        cout << "*********** Via Nodes Summary ***********" << endl;
//        cout << "Minimum : " << *std::min_element(via_numbers.begin(), via_numbers.end()) << endl;
//        cout << "Maximum : " << *std::max_element(via_numbers.begin(), via_numbers.end()) << endl;
//        int count = 0 ;
//        for(auto x:via_numbers)
//            count += x;
//        cout << "Average : " << count / total << endl;
//        cout << endl;
//
//        cout << "*********** Via Nodes after distance puring Summary ***********" << endl;
//        cout << "Minimum : " << *std::min_element(via_numbers_after_distance_puring.begin(), via_numbers_after_distance_puring.end()) << endl;
//        cout << "Maximum : " << *std::max_element(via_numbers_after_distance_puring.begin(), via_numbers_after_distance_puring.end()) << endl;
//        count = 0 ;
//        for(auto x:via_numbers_after_distance_puring)
//            count += x;
//        cout << "Average : " << count / total << endl;
//        cout << endl;
//
//        cout << "*********** Via Nodes after detour filtering Summary ***********" << endl;
//        cout << "Minimum : " << *std::min_element(via_numbers_after_detour_filtering.begin(), via_numbers_after_detour_filtering.end()) << endl;
//        cout << "Maximum : " << *std::max_element(via_numbers_after_detour_filtering.begin(), via_numbers_after_detour_filtering.end()) << endl;
//        count = 0 ;
//        for(auto x:via_numbers_after_detour_filtering)
//            count += x;
//        cout << "Average : " << count / total << endl;
//        cout << endl;
//
//        cout << "*********** Via Nodes after similarity filtering Summary ***********" << endl;
//        cout << "Minimum : " << *std::min_element(via_numbers_after_similarity_filtering.begin(), via_numbers_after_similarity_filtering.end()) << endl;
//        cout << "Maximum : " << *std::max_element(via_numbers_after_similarity_filtering.begin(), via_numbers_after_similarity_filtering.end()) << endl;
//        count = 0 ;
//        for(auto x:via_numbers_after_similarity_filtering)
//            count += x;
//        cout << "Average : " << count / total << endl;
//        cout << endl;
//
//
//        cout << "Start measuring " << endl;
//        double avg_length = 0;
//        for(auto x:alternative_Paths){
//            // cout << endl;
//            std::sort(begin(x) , end(x), [](const tuple < vector<unsigned > , unsigned, double >& a,
//                                            const tuple < vector<unsigned > , unsigned, double >& b) -> bool
//            {
//                return std::get<2>(a) > std::get<2>(b);
//            });
//            /*for(int i = 0; i < x.size() ; i++){
//                cout << get<1>(x[i]) << "  " << get<2>(x[i]) << endl;
//            }*/
//
//            double max = 0;
//            for(int i = 0; i < num_of_alternative_paths ; i++){
//                for(int j = 0; j < num_of_alternative_paths ; j ++){
//                    if(i > j){
//                        double score = similarity_check(get<0>(x[i]), get<0>(x[j]), get<1>(x[i]), get<1>(x[j]), first_out, head, weight);
//                        // cout << score << endl;
//                        if (score > max){
//                            max = score;
//                        }
//                    }
//                }
//            }
//            s_score.push_back(max);
//
//            double bs_max = 0;
//            double lo_min = 1;
//            double bs_score, lo_score;
//            double total_lengths = 0;
//            for(int i = 0; i < num_of_alternative_paths ; i++){
//                // cout << score << endl;
//                std::tie(bs_score, lo_score) = get_Local_Optimality_and_Bounded_Stretch(get<0>(x[i]), lab, first_out, head, weight);;
//                // cout << bs_score << "   " << lo_score << endl;
//                if (bs_score > bs_max){
//                    bs_max = bs_score;
//                }
//                if (lo_score < lo_min){
//                    lo_min = lo_score;
//                }
//                total_lengths += get<1>(x[i]);
//            }
//            // cout << "Results "  << bs_max << "  " << lo_min << endl;
//            lo_scores.push_back(lo_min);
//            b_s_score.push_back(bs_max);
//            avg_length += (total_lengths / 3);
//        }
//
//        cout << "***********  Similarity score for 3 alternative paths per query ***********" << endl;
//        cout << "Minimum : " << *std::min_element(s_score.begin(), s_score.end()) << endl;
//        cout << "Maximum : " << *std::max_element(s_score.begin(), s_score.end()) << endl;
//        double s_count = 0 ;
//        for(auto x:s_score)
//            s_count += x;
//        cout << "Average : " <<  s_count / s_score.size() << endl;
//        cout << endl;
//
//        cout << "***********  Bounded Stretch score for 3 alternative paths per query ***********" << endl;
//        cout << "Minimum : " << *std::min_element(b_s_score.begin(), b_s_score.end()) << endl;
//        cout << "Maximum : " << *std::max_element(b_s_score.begin(), b_s_score.end()) << endl;
//        double b_s_count = 0 ;
//        for(auto x:b_s_score)
//            b_s_count += x;
//        cout << "Average : " <<  b_s_count / b_s_score.size() << endl;
//        cout << endl;
//
//        cout << "***********  Local Optimality score for 3 alternative paths per query ***********" << endl;
//        cout << "Minimum : " << *std::min_element(lo_scores.begin(), lo_scores.end()) << endl;
//        cout << "Maximum : " << *std::max_element(lo_scores.begin(), lo_scores.end()) << endl;
//        double lo_count = 0 ;
//        for(auto x:lo_scores)
//            lo_count += x;
//        cout << "Average : " <<  lo_count / lo_scores.size() << endl;
//        cout << endl;

        return 0;
    }catch(exception&err){
        cerr << "Stopped on exception : " << err.what() << endl;
    }
}