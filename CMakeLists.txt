cmake_minimum_required(VERSION 3.13)
project(RoutingKit)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_CXX_FLAGS "-pedantic -Wall -Wno-strict-aliasing -Wno-long-long -Wno-deprecated -Wno-deprecated-declarations")
#set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")
#set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fPIC -Iinclude")

#set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG  -mcpu=apple-m1")
#set(DEV_CXX_FLAGS "-g -ggdb -O0 -fno-omit-frame-pointer")
#set(RELEASE_CXX_FLAGS "-O3 -DNDEBUG")

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -mcpu=apple-m1")
set(DEV_CXX_FLAGS "-g -ggdb -O0 -fno-omit-frame-pointer")
set(RELEASE_CXX_FLAGS "-O3 -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${DEV_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${RELEASE_CXX_FLAGS}")


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
set(ZLIB_INCLUDE_DIR "/opt/homebrew/opt/zlib/include")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include -fopenmp")



#include_directories("/usr/local/Cellar/libomp/10.0.0/include")
include_directories("/opt/homebrew/Cellar/libomp/14.0.6/include")
include_directories("/opt/homebrew/Cellar/google-sparsehash/2.0.4/include")

include_directories(include)
include_directories(src)

set(SRC
        include/routingkit/my_timer.h
        include/routingkit/bit_vector.h
        include/routingkit/constants.h
        include/routingkit/contraction_hierarchy.h
        include/routingkit/dijkstra.h
        include/routingkit/filter.h
        include/routingkit/geo_dist.h
        include/routingkit/geo_position_to_node.h
        include/routingkit/graph_util.h
        include/routingkit/id_mapper.h
        include/routingkit/id_queue.h
        include/routingkit/id_set_queue.h
        include/routingkit/inverse_vector.h
        include/routingkit/min_max.h
        include/routingkit/nested_dissection.h
        include/routingkit/osm_decoder.h
        include/routingkit/osm_graph_builder.h
        include/routingkit/osm_profile.h
        include/routingkit/osm_simple.h
        include/routingkit/permutation.h
        include/routingkit/sort.h
        include/routingkit/strongly_connected_component.h
        include/routingkit/tag_map.h
        include/routingkit/timer.h
        include/routingkit/timestamp_flag.h
        include/routingkit/vector_io.h
        include/routingkit/shortest_dijkstra_single.h
        include/routingkit/ch_dijkstra_single.h
        include/routingkit/ch_dijkstra_hl.h
        include/routingkit/ori_dijkstra.h

        src/searchnode.h
        src/modified_dijkstra.h
        src/cpool.h
        src/landmark.h
        src/my_timer.cpp
        src/bit_select.cpp
        src/bit_vector.cpp
        src/buffered_asynchronous_reader.cpp
        src/tree_node.h
        #    src/compare_vector.cpp
        #    src/compute_contraction_hierarchy.cpp
        #    src/compute_geographic_distance_weights.cpp
        #    src/compute_nested_dissection_order.cpp
        src/contraction_hierarchy.cpp
        #    src/convert_road_dimacs_coordinates.cpp
        #    src/convert_road_dimacs_graph.cpp
        #    src/decode_vector.cpp
        #    src/encode_vector.cpp
        #    src/examine_ch.cpp
        src/expect.cpp
        #    src/export_road_dimacs_graph.cpp
        src/file_data_source.cpp
        #    src/generate_constant_vector.cpp
        #    src/generate_dijkstra_rank_test_queries.cpp
        #    src/generate_random_node_list.cpp
        #    src/generate_random_source_times.cpp
        src/geo_position_to_node.cpp
        #    src/graph_to_dot.cpp
        #    src/graph_to_svg.cpp
        src/graph_util.cpp
        src/id_mapper.cpp
        src/nested_dissection.cpp
        src/osm_decoder.cpp
        #   src/osm_extract.cpp
        src/osm_graph_builder.cpp
        src/osm_profile.cpp
        src/osm_simple.cpp
        src/protobuf.cpp
        #    src/randomly_permute_nodes.cpp
        #    src/run_contraction_hierarchy_query.cpp
        #    src/run_dijkstra.cpp
        #    src/show_path.cpp
        src/strongly_connected_component.cpp
        src/timer.cpp
        src/vector_io.cpp
        src/verify.cpp


        src/bit_select.h
        src/buffered_asynchronous_reader.h
        src/emulate_gcc_builtin.h
        src/expect.h
        src/file_data_source.h
        src/protobuf.h
        src/verify.h
        src/CPD.h
        src/CPD.cpp
        src/DCPD.h
        src/DCPD.cpp
        src/my_timer.h
        src/my_timer.cpp

        src/construction.h
        src/coverage_ordering.h
        src/coverage_ordering_bp.h
        src/coverage_ordering_compress.h
        src/coverage_ordering_path.h
        src/graph.h
        src/graph_search.h
        src/heap.h
        src/labels.h
        src/ordering.h
        src/paras.h
        src/sse2neon.h
        src/time_util.h

        )

add_executable(convert_graph ${SRC} src/convert_road_dimacs_graph.cpp)
#add_executable(construct ${SRC} construct.cpp)
add_executable(convert_graph_with_speed ${SRC} src/convert_road_dimacs_graph_with_speed.cpp)
add_executable(convert_coordinates ${SRC} src/convert_road_dimacs_coordinates.cpp)
#add_executable(run_dijkstra ${SRC} src/run_dijkstra.cpp)
add_executable(compute_contraction_hierarchy ${SRC} src/compute_contraction_hierarchy.cpp)
add_executable(convert_osm ${SRC} convert_osm.cpp)
add_executable(benchmark_altervative_using_hublabeling_admissible ${SRC} benchmark_altervative_using_hublabeling_admissible.cpp)
add_executable(compute_landmark ${SRC} compute_landmark.cpp)

find_library(omp_lib omp "/opt/homebrew/opt")
if (omp_lib)
    target_link_libraries(benchmark_altervative_using_hublabeling_admissible  ${omp_lib})
endif()


find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
target_link_libraries(convert_graph ZLIB::ZLIB Threads::Threads)
target_link_libraries(convert_graph_with_speed ZLIB::ZLIB Threads::Threads)
target_link_libraries(convert_coordinates ZLIB::ZLIB Threads::Threads)
target_link_libraries(compute_contraction_hierarchy ZLIB::ZLIB Threads::Threads)
target_link_libraries(convert_osm  ZLIB::ZLIB Threads::Threads)
target_link_libraries(benchmark_altervative_using_hublabeling_admissible  ZLIB::ZLIB Threads::Threads)
target_link_libraries(compute_landmark ZLIB::ZLIB Threads::Threads)