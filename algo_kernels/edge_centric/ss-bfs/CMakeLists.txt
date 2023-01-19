set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -Wall -Wextra -Wl,--no-as-needed -funroll-loops -march=native -flto -pedantic -std=c++14")
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}")

file(GLOB_RECURSE SRC_BFS ${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/common_harness/*.cpp)
include_directories(${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/edge_centric/ss-bfs)

add_executable(bfs-ss ${SRC_BFS})

set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)

find_package(Boost  COMPONENTS system thread serialization)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR}) 
  #message( "Found Boost at ${Boost_LIBRARIES}")
else()
   message(FATAL_ERROR "Failed to find boost at " ${Boost_Lib_Dirs} " Need system thread")
endif()

target_compile_definitions(bfs-ss PRIVATE -DSSHEBACKEND)

target_include_directories(bfs-ss PRIVATE "${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")

#############################################
#            Link our libraries             #
#############################################
target_link_libraries(bfs-ss cryptoTools TaskHandler)
