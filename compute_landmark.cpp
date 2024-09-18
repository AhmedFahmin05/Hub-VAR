//
// Created by Bojie Shen on 31/8/20.
//
#include <routingkit/vector_io.h>
#include <routingkit/my_timer.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/min_max.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/dijkstra.h>
#include <modified_dijkstra.h>
#include "verify.h"
#include "CPD.h"
#include <routingkit/timer.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <omp.h>
#include <unordered_map>
#include "landmark.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace RoutingKit;
using namespace std;


vector<landmark> build_landmark(const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight, int number){
    std::vector<landmark> l_list;

    int random =0;

    auto tail = invert_inverse_vector(first_out);
    Dijkstra dij(first_out,tail,head);
    std::vector<unsigned> d_list =dij.get_all_distance(random,weight);
    int max_id =0;
    unsigned  max = 0;
    int cur_id = 0;
    for(unsigned d: d_list){
        if(d > max){
            max = d;
            max_id = cur_id;
        }
        cur_id ++;
    }

    for(int i = 0 ; i < number; i++){
        d_list =dij.get_all_distance(max_id,weight);
        landmark l  = landmark(d_list,max_id);
        l_list.push_back(l);

        max_id = 0;
        max = 0;
        for(int j = 0; j < first_out.size()-1; j ++){
            unsigned cur_distance =0;
            for( const landmark& ll : l_list){
                if(ll.id == j){
                    cur_distance = 0;
                    break;
                }
                cur_distance = cur_distance + ll.get_distance(j);
            }
            if(cur_distance > max) {
                max= cur_distance;
                max_id = j;
            }
        }
    }
    return l_list;
}

void export_landmark(std::ostream& outfile,std::vector<landmark> l_list){
    for(const landmark& ll : l_list){
        outfile << ll.id << " "<< ll.get_distance_list().size() ;
        for (unsigned d : ll.get_distance_list()){
            outfile << " "<< d ;
        }
        outfile << "\n";
    }
}


void build_landmark_for_top_n(string ch_file, string percentage, string number_of_landmark){
    cout << "Load top "<< percentage <<" percentage graph ... " << endl;


    ch_file = ch_file +"_"+ percentage;

    vector<unsigned> first_out = load_vector<unsigned>(ch_file+".first");
    vector<unsigned> head = load_vector<unsigned>(ch_file+".head");
    vector<unsigned> weight = load_vector<unsigned>(ch_file+".weight");
    cout << "Build landmark for top "<< percentage <<" percentage graph ... " << endl;

//        std::filebuf fb;
//        std::string file = ch_file+".landmark";
//        fb.open (file,std::ios::out);
//        std::ostream os(&fb);
    auto l_list = build_landmark(first_out,head,weight,stoi(number_of_landmark));
//        export_landmark(os, build_landmark(first_out,head,weight,stoi(number_of_landmark)));
//        fb.close();

    vector<unsigned > merged_l_list =vector<unsigned>(stoi(number_of_landmark)*(first_out.size()-1));

    for(int i = 0; i< l_list.size(); i ++){
        vector<unsigned> current = l_list[i].get_distance_list();
        for(int j = 0; j < current.size();j++){
            merged_l_list[i*(first_out.size()-1)+j] = current[j];
        }
    }
    save_vector(ch_file+".landmark"+number_of_landmark,merged_l_list);

    cout << "done" << endl;
}


void build_landmark_for_entire_graph(string ch_file, string number_of_landmark){
    cout << "Load entire graph ... " << endl;

//    vector<unsigned> first_out = load_vector<unsigned>(ch_file+".first");
//    vector<unsigned> head = load_vector<unsigned>(ch_file+".head");
//    vector<unsigned> weight = load_vector<unsigned>(ch_file+".weight");

    vector<unsigned> first_out ;
    vector<unsigned> head ;
    vector<unsigned> weight;

    ContractionHierarchy ch = ContractionHierarchy::load_file(ch_file+ ".ch");
    cout << "Build landmark for entire graph ... " << endl;
    ch.convert_to_graph(first_out,head,weight);
    auto l_list = build_landmark(first_out,head,weight,stoi(number_of_landmark));

    vector<unsigned > merged_l_list =vector<unsigned>(stoi(number_of_landmark)*(first_out.size()-1));

    unsigned landmarks = stoi(number_of_landmark);
    for(int i = 0; i< first_out.size()-1; i ++){
        for(int j = 0; j < l_list.size();j++){
            merged_l_list[i*landmarks+j] =  l_list[j].get_distance_list()[i];
        }
    }

    save_vector(ch_file+".landmark"+number_of_landmark,merged_l_list);

    cout << "done" << endl;
}


int main(int argc, char*argv[]){

    try{
        string ch_file;
        string percentage;
        string number_of_landmark;

        if(argc != 3){
            cerr << argv[0] << "ch_file number_of_landmark" << endl;
            return 1;
        }else{
            ch_file = argv[1];
            // percentage = argv[2];
            number_of_landmark = argv[2];
//            order_file = argv[3];
//            mapper_file = argv[4];
        }

        build_landmark_for_entire_graph(ch_file,number_of_landmark);

//
//        cout << "Load top "<< percentage <<" percentage graph ... " << endl;
//
//
//        ch_file = ch_file +"_"+ percentage;
//
//        vector<unsigned> first_out = load_vector<unsigned>(ch_file+".first");
//        vector<unsigned> head = load_vector<unsigned>(ch_file+".head");
//        vector<unsigned> weight = load_vector<unsigned>(ch_file+".weight");
//        cout << "Build landmark for top "<< percentage <<" percentage graph ... " << endl;
//
////        std::filebuf fb;
////        std::string file = ch_file+".landmark";
////        fb.open (file,std::ios::out);
////        std::ostream os(&fb);
//        auto l_list = build_landmark(first_out,head,weight,stoi(number_of_landmark));
////        export_landmark(os, build_landmark(first_out,head,weight,stoi(number_of_landmark)));
////        fb.close();
//
//        vector<unsigned > merged_l_list =vector<unsigned>(stoi(number_of_landmark)*(first_out.size()-1));
//
//        for(int i = 0; i< l_list.size(); i ++){
//            vector<unsigned> current = l_list[i].get_distance_list();
//            for(int j = 0; j < current.size();j++){
//                merged_l_list[i*(first_out.size()-1)+j] = current[j];
//            }
//        }
//        save_vector(ch_file+".landmark"+number_of_landmark,merged_l_list);
//
//        cout << "done" << endl;
    }catch(exception&err){
        cerr << "Stopped on exception : " << err.what() << endl;
    }
}
