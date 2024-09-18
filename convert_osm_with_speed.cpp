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
#include <routingkit/dijkstra.h>
using namespace RoutingKit;
using namespace std;

void build_shp_undirect_graph(const char *graphFileName,const char * labelFileName,const char * output){
    RoutingKit::my_timer t1 = RoutingKit::my_timer();
    t1.start();
    WGraph wgraph;
    wgraph.load_graph(graphFileName);
    Coverage_Ordering_Path coverage_ordering(wgraph);

    string orderFileName(labelFileName);
    orderFileName.append(".order");
    coverage_ordering.save_rank(orderFileName.c_str());

    double average_label  = 0 ;
    unsigned long long memory = 0;
    unsigned long long original_memory = 0;
    string labelFile(labelFileName);
    labelFile.append(".label");
    coverage_ordering.plabels.save_labels(labelFile.c_str(),average_label,memory,original_memory);
    t1.stop();
    double time = t1.elapsed_time_mins();

//    string output_file = argv[3];
    string output_file = output;
    std::ofstream myFile(output_file + "index/shp_construction.csv");
//    std::cout<<output_file + "index/shp_construction.csv"<<std::endl;
    myFile<<"time,average_label,memory,original_memory\n";
    myFile<<std::fixed<<setprecision(8)<<time<<","<<average_label
          <<","<<memory/(double)1000000
          <<","<<original_memory/(double)1000000
          <<"\n";
    myFile.close();




}

struct myComp {
    constexpr bool operator()(
            tuple<unsigned , unsigned > const& a,
            tuple<unsigned , unsigned > const& b)
    const noexcept
    {
        return get<0>(a) > get<0>(b);
    }
};




void path_planning(unsigned s, unsigned t, const unordered_map<unsigned, unsigned>& idx,
                   std::vector<RoutingKit::tree_node>& col,
                   std::vector<unsigned>& path) {
    if (col[idx.at(t) - 1].get_g() != INT32_MAX) {
        unsigned v = t;
        path.push_back(t);
        while (v != s) {
            col[idx.at(v) - 1].set_expanded(true);
            v = col[idx.at(v) - 1].get_parent();
            path.push_back(v);
        }
        col[idx.at(v) - 1].set_expanded(true);
        std::reverse(path.begin(), path.end());
        assert(path.front() == s && path.back() == t);
    }
}

void get_path_using_vai(unsigned s, unsigned t, unsigned u, const unordered_map<unsigned, unsigned>&  f_idx,
                        const std::vector<RoutingKit::tree_node>& ft, const unordered_map<unsigned, unsigned>&  b_idx,
                        const std::vector<RoutingKit::tree_node>& bt,
                        std::vector<unsigned>& ap) {
    // get path s-> .. -> u
    unsigned v = u;
    while (v != s) {
        ap.push_back(v);
        v = ft[f_idx.at(v) - 1].get_parent();
    }
    ap.push_back(v);

    std::reverse(ap.begin(), ap.end());

    assert(ap.front() == s && ap.back() == u);
    // get path u -> .. -> t
    v = u;
    while (v != t) {
        v = bt[b_idx.at(v) - 1].get_parent();
        ap.push_back(v);
    }
    assert(ap.front() == s && ap.back() == t);
}




bool comp(tuple<int, int> &a, tuple<int, int> &b )
{
    return (get<0>(a) == get<0>(b)) && (get<1>(a) == get<1>(b));
}



int main(int argc, char*argv[]){

    try{
        string file_location;
        string query_location;
        string ch_file;
        string first_out_file;
        string head_file;
        string weight_file;
        int query_count;
        int num_of_alternative_paths;
//        string result_file;
//        string percentage;
//        string number_of_landmark;



       if(argc != 2){
            cerr << argv[0] << "ch_file query_count  " << endl;
            return 1;
        }else{
            file_location = argv[1];
            //query_location = argv[2];
            //ch_file = argv[1];
            //first_out_file = argv[2];
            //head_file = argv[3];
            //weight_file = argv[4];
            //query_count = strtol(argv[3], NULL, 10);
            //num_of_alternative_paths = strtol(argv[4], NULL, 10);
//            result_file = argv[2];
//            percentage = argv[3];
//            number_of_landmark = argv[4];
//            order_file = argv[3];
//            mapper_file = argv[4];
        }


        std::ifstream gr(file_location);
        ofstream outfile;
        string file_path = "./dataset/MEL/AUS-road-t.MEL_speed.gr";
        outfile.open(file_path, ios::out);










        // V: number of nodes
        // E: number of edges
        unsigned V, E;
        if (!(gr >> V >> E)) {
            std::cout << "Error to get V and E \n";
            return 1;
        }
        std::cout << V << " " << E << endl;
        outfile << "p sp "<< V << " " << E * 2 << endl;
        unsigned x, y, speed, j;
        double wt;
        for (unsigned i = 0; i < E; i++) {
            if (!(gr >> j >> x >> y >> wt >> speed)) {
                //cout << j << " " << x << " " << y << " " << wt << " " << speed;
                std::cout << "Error to read graph file \n";
                return 1;
            }
            // From Km/hr to m/s
            double s = speed / 3.6;
            double val = ( wt / s ) * 1000;
            unsigned v = round(val);
            //g.addEdge(x , y , test, speed);
            // cout << x + 1 << " " << y + 1 << " " << val << " " << v <<  endl;
            outfile << "a " << x + 1 << " " << y + 1 << " "  << v <<  " " << speed <<  endl;
            outfile << "a " << y + 1 << " " << x + 1 << " "  << v << " " << speed <<  endl;

        }
        gr.close();
        outfile.close();
        return 0;

    }catch(exception&err){
        cerr << "Stopped on exception : " << err.what() << endl;
    }
}