//
// Created by Ahmed Fahmin on 14/12/2022.
//

#ifndef ROUTINGKIT_PERFORMANCE_METRICS_H
#define ROUTINGKIT_PERFORMANCE_METRICS_H
#include "../src/coverage_ordering_path.h"

int getRandom(int a,int b){
    int num = (rand() % (b - a + 1)) + a;
    return num;
}

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

unsigned get_edge_weight(unsigned s, unsigned t , const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight){
    unsigned edge_weight = INT32_MAX;
    for (unsigned j = first_out[s]; j < first_out[s + 1]; j++) {
        if ( head[j] == t && edge_weight > weight[j]) {
            edge_weight = weight[j];
        }
    }
    return edge_weight;
}

unsigned get_edge_speed(unsigned s, unsigned t , const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& speed){
    unsigned edge_weight = INT32_MAX;
    for (unsigned j = first_out[s]; j < first_out[s + 1]; j++) {
        if ( head[j] == t && edge_weight > speed[j]) {
            edge_weight = speed[j];
        }
    }
    return edge_weight;
}

unsigned get_edge_distance(unsigned s, unsigned t , const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned >& distance){
    unsigned edge_weight = INT32_MAX;
    for (unsigned j = first_out[s]; j < first_out[s + 1]; j++) {
        if ( head[j] == t && edge_weight > distance[j]) {
            edge_weight = distance[j];
        }
    }
    return edge_weight;
}

unsigned get_path_weight(const vector<unsigned>& p, const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight){
    unsigned total_weight = 0;
    for (unsigned i = 0; i < p.size() - 1; i++) {
        unsigned s = p[i];
        unsigned t = p[i + 1];
        total_weight += get_edge_weight(s, t, first_out, head, weight);
    }
    return total_weight;
}

unsigned get_total_path_distance(const vector<unsigned>& p, const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& distance){
    unsigned total_weight = 0;
    for (unsigned i = 0; i < p.size() - 1; i++) {
        unsigned s = p[i];
        unsigned t = p[i + 1];
        total_weight += get_edge_distance(s, t, first_out, head, distance);
    }
    return total_weight;
}
double similarity_check(const vector<unsigned>& p1, const vector<unsigned>& p2, unsigned distance1, unsigned distance2
        ,const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight){

    if(p1.size() > p2.size()) {
        return similarity_check(p2, p1, distance2, distance1, first_out, head, weight);
    }

    std::unordered_map<unsigned , unsigned> lookup;
    unsigned pos = 1;
    for (const auto& i : p1) {
        lookup[i] = pos;
        pos++;
    }
    unsigned common = 0;
    for (unsigned i = 0; i < p2.size() - 1; i++) {
        if(lookup.find(p2[i]) == lookup.end() || lookup.find(p2[i + 1]) == lookup.end() ){
            // not existed
            continue;
        }
        if ((lookup.at(p2[i]) + 1) == lookup.at(p2[i + 1])) {
            common += get_edge_weight(p2[i], p2[i + 1], first_out, head, weight);
        }
    }
    unsigned union_ = distance1 + distance2 - common;
    assert(union_ >= common);
    return (double) common/union_;
}

double similarity_check(const vector<unsigned>& p1, const vector<unsigned>& p2, unsigned distance1, unsigned distance2
        ,const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight, const double threshold){

    if(p1.size() > p2.size()) {
        return similarity_check(p2, p1, distance2, distance1, first_out, head, weight, threshold);
    }

    std::unordered_map<unsigned , unsigned> lookup;
    double upper_limit = (threshold / (1 + threshold)) * ((distance1 + distance2)) ;
    unsigned pos = 1;
    for (const auto& i : p1) {
        lookup[i] = pos;
        pos++;
    }
    unsigned common = 0;
    for (unsigned i = 0; i < p2.size() - 1; i++) {
        if(lookup.find(p2[i]) == lookup.end() || lookup.find(p2[i + 1]) == lookup.end() ){
            // not existed
            continue;
        }
        if ((lookup.at(p2[i]) + 1) == lookup.at(p2[i + 1])) {
            common += get_edge_weight(p2[i], p2[i + 1], first_out, head, weight);
        }
        if(common > upper_limit) {
            return -1;
        }
    }
    unsigned union_ = distance1 + distance2 - common;
    assert(union_ >= common);
    return (double) common/union_;
}

tuple<double, double> get_Local_Optimality_and_Bounded_Stretch_Optimal(const vector<unsigned >& p, unsigned length, const PLabel &lab, const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight,  const unsigned shortest_distance, const int step){
    unsigned val;
    unsigned s_distance;
    double min_lo_score = 1;
    double max_bs_score = 1;
    unsigned dist = shortest_distance;
    // unsigned dist = lab.query(p[0], p[p.size() - 1]);
    unsigned cost = 0;
    vector<unsigned> cumulative_distance(p.size(), 0);
    for(unsigned i = 1 ; i < p.size() ; i++){
        // cout << lab.query(p[i-1], p[i]) << "  " << get_edge_weight(p[i-1], p[i], first_out, head, weight) << " " <<endl;
        cost +=  get_edge_weight(p[i-1], p[i], first_out, head, weight);
        cumulative_distance[i] = cost;
    }

    assert(p.size() == cumulative_distance.size());

    double lo_score;
    double bs_score;
    //cout << "Start " << p.size() << "  "  << length << endl;
    for (int m = 0; m < length ; m = m + step) {
        //cout << p.size() - 1 << " " << length << " " << step << endl;
        for(int n = p.size() - 1 ; n >= length ; n  -= step){

            s_distance = lab.query(p[m], p[n]);
            //cout << n << " " << length <<  " " << " " << s_distance << endl;
            val = cumulative_distance[n] - cumulative_distance[m];
            if(s_distance == val){
                break;
            }
            // cout << s_distance << "  " << val << endl;
            bs_score = (double) val/s_distance;
            if(bs_score > max_bs_score){
                max_bs_score = bs_score;
            }
            if(val > s_distance){
                lo_score = (double) val/dist;
                if(lo_score < min_lo_score){
                    min_lo_score = lo_score;
                }
            }
        }
    }
    //cout << max_bs_score << " " << min_lo_score << endl;
    return make_tuple(max_bs_score, min_lo_score);
}

tuple<double, double> get_Local_Optimality_and_Bounded_Stretch_Optimal_Caching(const vector<unsigned >& p, unsigned length, const PLabel &lab, const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight, std::unordered_map<std::pair<unsigned, unsigned>, unsigned, pair_hash>& hub_lookup, const unsigned shortest_distance, unsigned step){
    unsigned val;
    unsigned s_distance;
    double min_lo_score = 1;
    double max_bs_score = 1;
    //unsigned dist = lab.query(p[0], p[p.size() - 1]);
    unsigned dist = shortest_distance;
    unsigned cost = 0;
    vector<unsigned> cumulative_distance(p.size(), 0);
    for(unsigned i = 1 ; i < p.size() ; i++){
        // cout << lab.query(p[i-1], p[i]) << "  " << get_edge_weight(p[i-1], p[i], first_out, head, weight) << " " <<endl;
        cost +=  get_edge_weight(p[i-1], p[i], first_out, head, weight);
        cumulative_distance[i] = cost;
    }

    assert(p.size() == cumulative_distance.size());

    double lo_score;
    double bs_score;
    for (unsigned m = 0; m < length ; m = m + step) {
        for(unsigned n = length ; n < p.size() ; n = n + step){

            if(hub_lookup.find(make_pair(p[m], p[n])) == hub_lookup.end()){
                s_distance = lab.query(p[m], p[n]);
                hub_lookup.insert({make_pair(p[m], p[n]), s_distance});

            }
            else{
                s_distance = hub_lookup.at(make_pair(p[m], p[n]));
            }
            //s_distance = lab.query(p[m], p[n]);

            val = cumulative_distance[n] - cumulative_distance[m];
            // cout << s_distance << "  " << val << endl;
            bs_score = (double) val/s_distance;
            if(bs_score > max_bs_score){
                max_bs_score = bs_score;
            }
            if(val > s_distance){
                lo_score = (double) val/dist;
                if(lo_score < min_lo_score){
                    min_lo_score = lo_score;
                }
            }
        }
    }
    // cout << max_bs_score << " " << min_lo_score << endl;
    return make_tuple(max_bs_score, min_lo_score);
}

tuple<double, double> get_Local_Optimality_and_Bounded_Stretch(const vector<unsigned >& p, const PLabel &lab, const vector<unsigned>& first_out, const vector<unsigned>& head, const vector<unsigned>& weight){
    unsigned val;
    unsigned s_distance;
    double min_lo_score = 1;
    double max_bs_score = 0;
    unsigned dist = lab.query(p[0], p[p.size() - 1]);
    //cout << dist << endl;

    unsigned cost = 0;
    vector<unsigned> cumulative_distance(p.size(), 0);
    for(unsigned i = 1 ; i < p.size() ; i++){
        //cout << lab.query(p[i-1], p[i]) << "  " << get_edge_weight(p[i-1], p[i], first_out, head, weight) << " " <<endl;
        cost +=  get_edge_weight(p[i-1], p[i], first_out, head, weight);
        cumulative_distance[i] = cost;
    }
    assert(p.size() == cumulative_distance.size());
    /*cout << "The Values : " ;
    for(unsigned i = 0 ; i < cumulative_distance.size(); i++){
        cout << cumulative_distance[i] << " " ;
    }
    cout << endl;*/
    double lo_score;
    double bs_score;
    for (unsigned m = 0; m < p.size() ; m++) {
        for(unsigned n = m + 1 ; n < p.size(); n ++){
            s_distance = lab.query(p[m], p[n]);
            val = cumulative_distance[n] - cumulative_distance[m];
            // cout << s_distance << "  " << val << endl;
            bs_score = (double) val/s_distance;
            if(bs_score > max_bs_score){
                max_bs_score = bs_score;
            }
            if(val > s_distance){
                lo_score = (double) val/dist;
                if(lo_score < min_lo_score){
                    min_lo_score = lo_score;
                }
            }
        }
    }
    // cout << max_bs_score << " " << min_lo_score << endl;
    return make_tuple(max_bs_score, min_lo_score);
}

double Bounded_Stretch(vector<int> p, PLabel lab, vector<unsigned> first_out, vector<unsigned> head, vector<unsigned > weight){
    int val;
    int s_distance;
    double score = -INFINITY;

    double cost = 0;
    vector<int> cumulative_distance;
    cumulative_distance.emplace_back(cost);
    for(int i = 1; i < p.size() ; i++){
        cost += get_edge_weight(p[i-1], p[i], first_out, head, weight);
        cumulative_distance.emplace_back(cost);
    }
    assert(p.size() == cumulative_distance.size());

    double t_score;
    for (int m = 0; m < p.size() - 2; m++) {
        for(int n = m + 2 ; n < p.size(); n ++){
            s_distance = lab.query(p[m], p[n]);
            val = cumulative_distance[n] - cumulative_distance[m];
            t_score = (double) val/s_distance;
            if(t_score > score){
                score = t_score;
            }
            //cout << p[m] << "  " << p[n] << "   " << s_distance << "   " << val << endl;
        }
        //cout << endl;
        /*
        int s = p[m];
        int t = p[m + 1];
        val += get_edge_weight(s, t, first_out, head, weight);
         */
    }
    return score;
}



double Local_Optimality(vector<int> p, PLabel lab, vector<unsigned> first_out, vector<unsigned> head, vector<unsigned > weight){
    int val;
    int s_distance;
    double score = 1;
    int dist = lab.query(p[0], p[p.size() - 1]);
    double cost = 0;
    vector<int> cumulative_distance;
    cumulative_distance.emplace_back(cost);
    for(int i = 1; i < p.size() ; i++){
        cost += get_edge_weight(p[i-1], p[i], first_out, head, weight);
        cumulative_distance.emplace_back(cost);
    }

    // assert(p.size() == cumulative_distance.size());

    double t_score;
    for (int m = 0; m < p.size() - 2; m++) {
        for(int n = m + 2 ; n < p.size(); n ++){
            s_distance = lab.query(p[m], p[n]);
            val = cumulative_distance[n] - cumulative_distance[m];
            if(val > s_distance){
                t_score = (double) val/dist;
                if(t_score < score){
                    score = t_score;
                }
            }
        }
    }
    return score;
}


#endif //ROUTINGKIT_PERFORMANCE_METRICS_H
