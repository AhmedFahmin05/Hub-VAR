//
// Created by Bojie Shen on 10/8/20.
//


#include <chrono>
#include "my_timer.h"

RoutingKit::my_timer::my_timer() = default;


double RoutingKit::my_timer::get_time_nano()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();
//#ifdef OS_MAC
//    uint64_t raw_time = mach_absolute_time();
//    return (double)(raw_time * timebase.numer / timebase.denom);
//#else
//    timespec raw_time;
//    clock_gettime(CLOCK_MONOTONIC , &raw_time);
//    return (double)(raw_time.tv_nsec);
//#endif
}

void RoutingKit::my_timer::start()
{
    start_time = std::chrono::steady_clock::now();
    stop_time = start_time;
//#ifdef OS_MAC
//    start_time = mach_absolute_time();
//    stop_time = start_time;
//#else
////    clock_gettime(CLOCK_MONOTONIC , &start_time);
//#endif
}

void RoutingKit::my_timer::stop()
{
    stop_time =  std::chrono::steady_clock::now();
//#ifdef OS_MAC
//    stop_time = mach_absolute_time();
//#else
////    clock_gettime(CLOCK_MONOTONIC , &stop_time);
////    auto stime = std::chrono::steady_clock::now();
////    auto etime = std::chrono::steady_clock::now();
//#endif
}

double RoutingKit::my_timer::current_time_nano(std::chrono::time_point<std::chrono::steady_clock> current_time){
    return std::chrono::duration_cast<std::chrono::nanoseconds>(current_time - start_time).count();
}
double RoutingKit::my_timer::elapsed_time_nano()
{

    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();

//#ifdef OS_MAC
//    uint64_t elapsed_time = stop_time - start_time;
//    return (double)(elapsed_time * timebase.numer / timebase.denom);
////    Nanoseconds nanosecs = AbsoluteToNanoseconds(*(AbsoluteTime*)&elapsed_time);
////    return (double) UnsignedWideToUInt64(nanosecs) ;
//
//#else
////    if ((stop_time.tv_nsec-start_time.tv_nsec)<0)
////        return (double)(1000000000+stop_time.tv_nsec-start_time.tv_nsec);
////    else
////        return (double)(stop_time.tv_nsec - start_time.tv_nsec);
//#endif
}

void RoutingKit::my_timer::reset()
{
//#ifdef OS_MAC
//    start_time = stop_time = 0;
//#else
////    start_time.tv_sec = 0;
////    start_time.tv_nsec = 0;
////    stop_time.tv_sec = 0;
////    stop_time.tv_nsec = 0;
//#endif
}

double RoutingKit::my_timer::elapsed_time_micro()
{
    return elapsed_time_nano() / 1000.0;
}
double RoutingKit::my_timer::elapsed_time_mins()
{
    return elapsed_time_nano() / 1000.0/ 1000000.0 / 60.0;;
}
