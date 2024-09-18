//
// Created by Bojie Shen on 31/8/20.
//

#pragma once
#include <vector>
#include <cmath>

namespace RoutingKit{
    class landmark{
    private:
        std::vector<unsigned> distance_list;



    public:
        int id;

        landmark(std::vector<unsigned> d_list,int l_id) : distance_list(d_list), id(l_id){ }


        unsigned get_lower_bound(const int& start_id, const int& end_id) const{
            const unsigned& start_distance  =distance_list[start_id];
            const unsigned& end_distance  =distance_list[end_id];
//            if(start_distance == -1 || end_distance == -1){
//                return -1;
//            } else {
////                unsigned result1 =  abs((long long int)start_distance - end_distance );
////                unsigned result2 = (start_distance > end_distance ) ? (start_distance - end_distance ) : (end_distance - start_distance);
//                return abs((long long int)start_distance - end_distance );;
//            }
            return abs((long long int)start_distance - end_distance );;
        }

//        unsigned get_upper_bound(const int& start_id, const int& end_id) const{
//            const unsigned& start_distance  =distance_list[start_id];
//            const unsigned& end_distance  =distance_list[end_id];
//            if(start_distance == -1 || end_distance == -1){
//                return -1;
//            } else {
//                return abs(start_distance+end_distance);
//            }
//        }

        unsigned get_distance (const int& centroid_id) const {
            if(centroid_id == id ){
                return 0;
            }else{
                return distance_list[centroid_id];
            }
        }

        const std::vector<unsigned>& get_distance_list () const {
            return distance_list;
        }
    };


}
