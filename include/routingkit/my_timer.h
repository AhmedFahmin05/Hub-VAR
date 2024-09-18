//
// Created by Bojie Shen on 10/8/20.
//

#ifndef ROUTINGKIT_MY_TIMER_H
#define ROUTINGKIT_MY_TIMER_H
#include <chrono>
namespace RoutingKit
{

    class my_timer
    {
        std::chrono::time_point<std::chrono::steady_clock> stop_time;
        std::chrono::time_point<std::chrono::steady_clock> start_time;
//#ifdef OS_MAC
//    uint64_t start_time;
//    uint64_t stop_time;
//    mach_timebase_info_data_t timebase;
//#else
//    static std::chrono::time_point<std::chrono::steady_clock> stop_time;
//    static std::chrono::time_point<std::chrono::steady_clock> start_time;
////    timespec stop_time;
////    timespec start_time;
//#endif

    public:
        my_timer();
        void reset();
        void start();
        void stop();

        double current_time_nano(std::chrono::time_point<std::chrono::steady_clock> current_time);
        double elapsed_time_nano();
        double elapsed_time_micro();
        double get_time_nano();

        double elapsed_time_mins();
    };

}
#endif //ROUTINGKIT_MY_TIMER_H
