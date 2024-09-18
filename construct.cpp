//
// Created by Bojie Shen on 22/11/20.
//


#include "../src/graph.h"
#include "../src/coverage_ordering_path.h"
#include "../src/my_timer.h"
#include <stdexcept>
#include <iomanip>

void build_shp_undirect_graph(const char *graphFileName,const char * labelFileName,const char * output){
    //RoutingKit::my_timer t1 = RoutingKit::my_timer();
    //t1.start();
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
    //t1.stop();
    //double time = t1.elapsed_time_mins();

//    string output_file = argv[3];
    string output_file = output;
    std::ofstream myFile(output_file + "index/shp_construction.csv");
//    std::cout<<output_file + "index/shp_construction.csv"<<std::endl;
//    myFile<<"time,average_label,memory,original_memory\n";
//    myFile<<std::fixed<<setprecision(8)<<time<<","<<average_label
//    <<","<<memory/(double)1000000
//    <<","<<original_memory/(double)1000000
//    <<"\n";
//    myFile.close();




}


int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "graph_file label_file output_file\n");
        return 0;
    }
    build_shp_undirect_graph(argv[1],argv[2],argv[3]);

    return 0;
}