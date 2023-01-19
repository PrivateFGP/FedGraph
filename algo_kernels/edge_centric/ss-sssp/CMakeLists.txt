set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -Wall -Wextra -Wl,--no-as-needed -funroll-loops -march=native -flto -pedantic -std=c++14")
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}")

file(GLOB_RECURSE SRC_SSSP ${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/common_harness/*.cpp)
include_directories(${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/edge_centric/ss-sssp)

# add_executable(sssp-sgx ${SRC_SSSP})
# target_compile_definitions(sssp-sgx PRIVATE -DSGXBACKEND)

# target_include_directories(sssp-sgx PRIVATE "${CMAKE_SOURCE_DIR}/SGX-Worker/Common/")

# # add_untrusted_executable(sssp SRCS ${SRC_SSSP} EDL ${CMAKE_SOURCE_DIR}/SGX-Worker/Enclave/Enclave.edl EDL_SEARCH_PATHS ${CMAKE_SOURCE_DIR}/SGX-Worker/Enclave/)

# #############################################
# #            Link our libraries             #
# #############################################
# target_link_libraries(sssp-sgx cryptoTools TaskHandler)

# # add_dependencies(sssp libenclave-sign)

add_executable(sssp-ss ${SRC_SSSP})

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

target_compile_definitions(sssp-ss PRIVATE -DSSHEBACKEND)

target_include_directories(sssp-ss PRIVATE "${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")

#############################################
#            Link our libraries             #
#############################################
target_link_libraries(sssp-ss cryptoTools TaskHandler)
