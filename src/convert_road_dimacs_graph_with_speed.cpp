#include <routingkit/vector_io.h>
#include <routingkit/permutation.h>
#include <routingkit/sort.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/graph_util.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RoutingKit;
using namespace std;

void sort_arcs_and_remove_multi_and_loop_arcs(
        unsigned node_count, std::vector<unsigned>&tail, std::vector<unsigned>&head, std::vector<unsigned>&weight, std::vector<unsigned>&speed, std::vector<unsigned>&input_arc_id
){

    long long timer = 0;  // initialize to avoid warning, not needed
    {
        auto p = compute_inverse_sort_permutation_first_by_tail_then_by_head_and_apply_sort_to_tail(node_count, tail, head);
        head = apply_inverse_permutation(p, std::move(head));
        weight = apply_inverse_permutation(p, std::move(weight));
        speed = apply_inverse_permutation(p, std::move(speed));
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
                speed[out] = speed[in];
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
                speed[out] = speed[in];
                input_arc_id[out] = input_arc_id[in];
                ++out;
            }else{
                if(weight[in] < weight[out-1]){
                    weight[out-1] = weight[in];
                    speed[out-1] = speed[in];
                    input_arc_id[out-1] = input_arc_id[in];
                }
            }
        }
        arc_count = out;
    }

    tail.erase(tail.begin()+arc_count, tail.end());
    head.erase(head.begin()+arc_count, head.end());
    weight.erase(weight.begin()+arc_count, weight.end());
    speed.erase(speed.begin()+arc_count, speed.end());
    input_arc_id.erase(input_arc_id.begin()+arc_count, input_arc_id.end());
}



int main(int argc, char*argv[]){

	try{
		string dimacs_file;
		string output_file;
		string first_out_file;
		string head_file;
		string weight_file;
        string speed_file;

		if(argc != 3){
			cerr << argv[0] << " dimacs_file" << endl;
			return 1;
		}else{
			dimacs_file = argv[1];
			output_file = argv[2];
//			first_out_file = argv[2];
//			head_file = argv[3];
//			weight_file = argv[4];
		}

        first_out_file = output_file +".first";
        head_file = output_file +".head";
        weight_file = output_file +".weight";
        speed_file = output_file +".speed";

		cout << "Loading data ... " << flush;

		ifstream in(dimacs_file);
		if(!in)
			throw runtime_error("Can not open \""+dimacs_file+"\"");

		string line;
		unsigned line_num = 0;
		unsigned next_arc = 0;

		vector<unsigned>tail, head, weight, speed;
		unsigned node_count, arc_count;

		bool was_header_read = false;
		while(std::getline(in, line)){
			++line_num;
			if(line.empty() || line[0] == 'c')
				continue;

			std::istringstream lin(line);
			if(!was_header_read){
				was_header_read = true;
				std::string p, sp;
				if(!(lin >> p >> sp >> node_count >> arc_count))
					throw std::runtime_error("Can not parse header in dimacs file.");
				if(p != "p" || sp != "sp")
					throw std::runtime_error("Invalid header in dimacs file.");

				tail.resize(arc_count);
				head.resize(arc_count);
				weight.resize(arc_count);
                speed.resize(arc_count);

			}else{
				std::string a;
				unsigned h, t, w, s;
				if(!(lin >> a >> t >> h >> w >> s))
					throw std::runtime_error("Can not parse line num "+std::to_string(line_num)+" \""+line+"\" in dimacs file.");
				--h;
				--t;
				if(a != "a" || h >= node_count || t >= node_count)
					throw std::runtime_error("Invalid arc in line num "+std::to_string(line_num)+" \""+line+"\" in dimacs file.");
				if(next_arc < arc_count){
					head[next_arc] = h;
					tail[next_arc] = t;
					weight[next_arc] = w;
                    speed[next_arc] = s;
				}
				++next_arc;
			}
		}

		if(next_arc != arc_count)
			throw std::runtime_error("The arc count in the header ("+to_string(arc_count)+") does not correspond with the actual number of arcs ("+to_string(next_arc)+").");

		cout << "done" << endl;

		cout << "Ordering arcs ... " << flush;

		vector<unsigned>first_out;
        std::vector<unsigned>input_arc_id = identity_permutation(head.size());


        {
            sort_arcs_and_remove_multi_and_loop_arcs(node_count, tail, head, weight, speed, input_arc_id);
        }

		{
			auto p = compute_inverse_stable_sort_permutation_using_key(head, node_count, [](unsigned x){return x;});
			tail = apply_inverse_permutation(p, tail);
			auto q = compute_inverse_stable_sort_permutation_using_key(tail, node_count, [](unsigned x){return x;});
			tail = apply_inverse_permutation(q, tail);
			auto r = chain_permutation_first_left_then_right(q, p);
			head = apply_inverse_permutation(r, head);
			weight = apply_inverse_permutation(r, weight);
            speed = apply_inverse_permutation(r, speed);
			first_out = invert_vector(tail, node_count);
		}

		cout << "done" << endl;

		cout << "Saving file ... " << flush;
		save_vector(first_out_file, first_out);
		save_vector(head_file, head);
		save_vector(weight_file, weight);
        save_vector(speed_file, speed);
		cout << "done" << endl;

	}catch(exception&err){
		cerr << "Stopped on exception : " << err.what() << endl;
	}
}
