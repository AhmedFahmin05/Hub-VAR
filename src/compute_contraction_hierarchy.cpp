#include <routingkit/vector_io.h>
#include <routingkit/timer.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/min_max.h>
#include <routingkit/inverse_vector.h>

#include "verify.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <routingkit/my_timer.h>
#include <iomanip>

using namespace RoutingKit;
using namespace std;

//int main(int argc, char*argv[]){
//
//	try{
//		string graph_first_out;
//		string graph_head;
//		string graph_weight;
//
//		string ch_file;
//
//		if(argc != 2){
//			cerr << argv[0] << " ch_file" << endl;
//			return 1;
//		}else{
//			ch_file = argv[1];
//		}
//        size_t lastindex = ch_file.find_last_of(".");
//        string rawname = ch_file.substr(0, lastindex);
//        graph_first_out = rawname +".first";
//        graph_head = rawname +".head";
//        graph_weight = rawname +".weight";
//
//		cout << "Loading graph ... " << flush;
//
//		vector<unsigned>first_out = load_vector<unsigned>(graph_first_out);
//		vector<unsigned>head = load_vector<unsigned>(graph_head);
//		vector<unsigned>weight = load_vector<unsigned>(graph_weight);
//
//		cout << "done" << endl;
//
//		cout << "Validity tests ... " << flush;
//		check_if_graph_is_valid(first_out, head);
//		cout << "done" << endl;
//
//
//		const unsigned node_count = first_out.size()-1;
//		const unsigned arc_count = head.size();
//
//		if(first_out.front() != 0)
//			throw runtime_error("The first element of first out must be 0.");
//		if(first_out.back() != arc_count)
//			throw runtime_error("The last element of first out must be the arc count.");
//		if(head.empty())
//			throw runtime_error("The head vector must not be empty.");
//		if(max_element_of(head) >= node_count)
//			throw runtime_error("The head vector contains an out-of-bounds node id.");
//		if(weight.size() != arc_count)
//			throw runtime_error("The weight vector must be as long as the number of arcs");
//
//		my_timer t1 = my_timer();
//		t1.start();
//		vector<unsigned> rank = {0,1,2,3,4,5,6,7,8,9};
////		auto ch = ContractionHierarchy::build(node_count, invert_inverse_vector(first_out), head, weight, [](string msg){cout << msg << endl;});
//        auto ch = ContractionHierarchy::build_given_rank
//                (rank, invert_inverse_vector(first_out), head, weight, [](string msg){cout << msg << endl;});
//
//		t1.stop();
//        double building_time = t1.elapsed_time_mins();
//
//		check_contraction_hierarchy_for_errors(ch);
//		unsigned long long ch_size = 0 ;
//		ch.save_file(ch_file);
//		ch.print_graph();
////        ch_size = ch.get_size()/1000000;
//
//
//        lastindex = ch_file.find_last_of("/");
//        rawname = ch_file.substr(0, lastindex);
//
//        string filename = ch_file.substr(lastindex);
//        lastindex = ch_file.find_last_of(".");
//        filename = ch_file.substr(0,lastindex);
//        lastindex =filename.find_last_of("/");
//        filename = filename.substr(lastindex+1);
//
//
//        string outputname = rawname+"/results/index/" +filename;
//
//        std::ofstream myFile(outputname+"_ch.csv");
//        myFile<<"time,size\n";
//        myFile<<std::fixed<<setprecision(8)<<building_time
//                  <<","<<ch.get_size()/(double)1000000
//                  <<"\n";
//        myFile.close();
//        std::cout<<ch.get_size()/(double)1000000<<std::endl;
//
//	}catch(exception&err){
//		cerr << "Stopped on exception : " << err.what() << endl;
//	}
//}


#include <routingkit/vector_io.h>
#include <routingkit/timer.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/min_max.h>
#include <routingkit/inverse_vector.h>

#include "verify.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <routingkit/my_timer.h>
#include <iomanip>

using namespace RoutingKit;
using namespace std;

int main(int argc, char*argv[]){

    try{
        string graph_first_out;
        string graph_head;
        string graph_weight;

        string ch_file;

        if(argc != 2){
            cerr << argv[0] << " ch_file" << endl;
            return 1;
        }else{
            ch_file = argv[1];
        }
        size_t lastindex = ch_file.find_last_of(".");
        string rawname = ch_file.substr(0, lastindex);
        graph_first_out = rawname +".first";
        graph_head = rawname +".head";
        graph_weight = rawname +".weight";
        my_timer t1 = my_timer();

        cout << "Loading graph ... " << flush;

//        vector<unsigned>first_out = load_vector<unsigned>(graph_first_out);
//        vector<unsigned>head = load_vector<unsigned>(graph_head);
//        vector<unsigned>weight = load_vector<unsigned>(graph_weight);

//        cout << "done" << endl;
//
//        cout << "Validity tests ... " << flush;
//        check_if_graph_is_valid(first_out, head);
//        cout << "done" << endl;
//
//
//        const unsigned node_count = first_out.size()-1;
//        const unsigned arc_count = head.size();
//
//        if(first_out.front() != 0)
//            throw runtime_error("The first element of first out must be 0.");
//        if(first_out.back() != arc_count)
//            throw runtime_error("The last element of first out must be the arc count.");
//        if(head.empty())
//            throw runtime_error("The head vector must not be empty.");
//        if(max_element_of(head) >= node_count)
//            throw runtime_error("The head vector contains an out-of-bounds node id.");
//        if(weight.size() != arc_count)
//            throw runtime_error("The weight vector must be as long as the number of arcs");

        t1.start();


        vector<unsigned>first_out = load_vector<unsigned>(graph_first_out);
        vector<unsigned>head = load_vector<unsigned>(graph_head);
        vector<unsigned>weight = load_vector<unsigned>(graph_weight);
//        cout<<first_out.size()<<endl;
//        cout<<head.size()<<endl;
//        cout<<head[733843]<<endl;
//        cout<<head[2909]<< " " << weight[2909] << endl;

        const unsigned node_count = first_out.size()-1;
        auto ch = ContractionHierarchy::build(node_count, invert_inverse_vector(first_out), head, weight, [](string msg){cout << msg << endl;});
        ch.save_file(ch_file);


        t1.stop();
        double building_time = t1.elapsed_time_mins();

        check_contraction_hierarchy_for_errors(ch);

        ch.save_file(ch_file);
//        ch_size = ch.get_size()/1000000;
        cout<<ch.forward.head.size()<<endl;
        cout<<ch.forward.is_shortcut_an_original_arc.size()<<endl;
        cout<<ch.backward.is_shortcut_an_original_arc.size()<<endl;
        cout<<ch.forward.shortcut_first_arc.size()<<endl;
        cout<<ch.forward.shortcut_second_arc.size()<<endl;
        cout<<ch.backward.shortcut_first_arc.size()<<endl;
        cout<<ch.backward.shortcut_second_arc.size()<<endl;


        int j = 0; int k=0;
        set<int> original_arc;
        for(int i = 0 ; i< ch.forward.head.size();i ++){
            if(ch.forward.is_shortcut_an_original_arc.is_set(i)){
                original_arc.insert(ch.forward.shortcut_first_arc[i]);
                j++;
            }else{
                k++;
            }
        }
        for(int i = 0 ; i< ch.backward.head.size();i ++){
            if(ch.backward.is_shortcut_an_original_arc.is_set(i)){
                original_arc.insert(ch.backward.shortcut_first_arc[i]);
                j++;
            }else{
                k++;
            }
        }
        cout<<original_arc.size()<<endl;
        cout<< " original arc:" << j << "shortcust: "<< k <<endl;
        lastindex = ch_file.find_last_of("/");
        rawname = ch_file.substr(0, lastindex);

        string filename = ch_file.substr(lastindex);
        lastindex = ch_file.find_last_of(".");
        filename = ch_file.substr(0,lastindex);
        lastindex =filename.find_last_of("/");
        filename = filename.substr(lastindex+1);


        string outputname = rawname+"/results/index/" +filename;

        std::ofstream myFile(outputname+"_ch.csv");
        myFile<<"time,size\n";
        myFile<<std::fixed<<setprecision(8)<<building_time
              <<","<<ch.get_size()/(double)1000000
              <<"\n";
        myFile.close();

    }catch(exception&err){
        cerr << "Stopped on exception : " << err.what() << endl;
    }
}

