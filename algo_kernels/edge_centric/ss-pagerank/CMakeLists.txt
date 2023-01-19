set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -Wall -Wextra -Wl,--no-as-needed -funroll-loops -march=native -flto -pedantic -std=c++14")
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}")

# add_executable(pagerank-sgx ${SRC_PAGERANK})
# target_compile_definitions(pagerank-sgx PRIVATE -DSGXBACKEND)

# target_include_directories(pagerank-sgx PRIVATE "${CMAKE_SOURCE_DIR}/SGX-Worker/Common/")

# # add_untrusted_executable(pagerank SRCS ${SRC_PAGERANK} EDL ${CMAKE_SOURCE_DIR}/SGX-Worker/Enclave/Enclave.edl EDL_SEARCH_PATHS ${CMAKE_SOURCE_DIR}/SGX-Worker/Enclave/)

# #############################################
# #            Link our libraries             #
# #############################################
# target_link_libraries(pagerank-sgx cryptoTools)
# target_link_libraries(pagerank-sgx TaskHandler)

# # get_target_property(OUT pagerank LINK_LIBRARIES)
# # message(STATUS ${OUT})


file(GLOB_RECURSE SRC_PAGERANK ${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/common_harness/*.cpp)
add_executable(pagerank-ss ${SRC_PAGERANK})

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

target_compile_definitions(pagerank-ss PRIVATE -DSSHEBACKEND)
target_include_directories(pagerank-ss PRIVATE "${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")
target_include_directories(pagerank-ss PRIVATE ${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/edge_centric/ss-pagerank)
target_link_libraries(pagerank-ss PRIVATE cryptoTools TaskHandler)

# find_package (Boost COMPONENTS thread system REQUIRED)

# FILE(GLOB_RECURSE TEST_COMM_SRCS ${CMAKE_SOURCE_DIR}/SS-Worker/test/*.cpp)
# add_executable(test-comm ${TEST_COMM_SRCS})
# target_compile_definitions(test-comm PRIVATE -DSSBACKEND)
# target_include_directories(test-comm PRIVATE "${CMAKE_SOURCE_DIR}/SS-Worker/include/")
# target_include_directories(test-comm PRIVATE "${Boost_INCLUDE_DIRS}")
# target_link_libraries(test-comm PRIVATE ${Boost_LIBRARIES} cryptoTools TaskHandler ${EMP-OT_LIBRARIES} OPHELib::ophelib)